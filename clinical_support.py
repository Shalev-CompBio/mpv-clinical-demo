"""
Clinical Support Engine - Main API for the Clinical Decision Support Framework.

Provides a unified interface for:
- Phenotype-based queries
- Gene-based queries
- Interactive decision support
"""

from typing import List, Set, Optional, Union
from pathlib import Path

from data_loader import DataLoader, get_loader
from scoring_engine import ScoringEngine, ScoringConfig
from prediction_engine import PredictionEngine, PredictionConfig
from output_models import (
    QueryResult, GeneQueryResult, ModuleMatch, GeneCandidate,
    PhenotypePrediction, PhenotypeInfo, ExplainabilityItem
)


class ClinicalSupportEngine:
    """
    Main entry point for the Clinical Decision Support Framework.
    
    Usage:
        engine = ClinicalSupportEngine()
        
        # Query with phenotypes
        result = engine.query(
            observed=["HP:0001513", "HP:0010442"],  # Obesity, Polydactyly
            excluded=["HP:0000365"]  # Hearing impairment
        )
        
        # Query with gene
        result = engine.query_gene("RPGR")
    """
    
    def __init__(self, 
                 data_dir: Optional[Path] = None,
                 scoring_config: Optional[ScoringConfig] = None,
                 prediction_config: Optional[PredictionConfig] = None):
        """
        Initialize the clinical support engine.
        
        Args:
            data_dir: Optional path to data directory containing Excel and CSV files
            scoring_config: Optional scoring configuration
            prediction_config: Optional prediction configuration
        """
        self.loader = DataLoader(data_dir)
        self.loader.load()
        
        self.scorer = ScoringEngine(self.loader, scoring_config)
        self.predictor = PredictionEngine(self.loader, prediction_config)
    
    def _resolve_phenotypes(self, 
                            phenotypes: List[str]) -> tuple[Set[str], List[str]]:
        """
        Resolve phenotype identifiers to HPO IDs.
        
        Returns:
            Tuple of (resolved HPO IDs, unmatched identifiers)
        """
        resolved = set()
        unmatched = []
        
        for identifier in phenotypes:
            hpo_id = self.loader.resolve_phenotype(identifier)
            if hpo_id:
                resolved.add(hpo_id)
            else:
                unmatched.append(identifier)
        
        return resolved, unmatched
    
    def query(self,
              observed: Optional[List[str]] = None,
              excluded: Optional[List[str]] = None,
              top_genes: int = 20,
              top_predictions: int = 10) -> QueryResult:
        """
        Query with observed and/or excluded phenotypes.
        
        Args:
            observed: List of phenotypes that are present (HPO IDs or names)
            excluded: List of phenotypes that are absent (HPO IDs or names)
            top_genes: Maximum number of candidate genes to return
            top_predictions: Maximum number of predicted phenotypes to return
        
        Returns:
            QueryResult with matched modules, candidate genes, and predictions
        """
        observed = observed or []
        excluded = excluded or []
        
        # Resolve phenotypes
        observed_ids, unmatched_obs = self._resolve_phenotypes(observed)
        excluded_ids, unmatched_exc = self._resolve_phenotypes(excluded)
        
        all_unmatched = unmatched_obs + unmatched_exc
        
        # Score and rank modules
        ranked_modules = self.scorer.rank_modules(observed_ids, excluded_ids)
        
        # Get best module
        best_module = ranked_modules[0] if ranked_modules else None
        
        # Rank genes within best module
        candidate_genes = []
        if best_module and best_module.score > 0:
            candidate_genes = self.scorer.rank_genes(
                best_module.module_id, observed_ids
            )[:top_genes]
        
        # Predict missing phenotypes
        predicted = []
        if best_module and best_module.score > 0:
            predicted = self.predictor.predict_missing_phenotypes(
                best_module.module_id, observed_ids, excluded_ids
            )[:top_predictions]
            
        # --- Alternative Candidates Logic (Safeguard) ---
        alternative_genes = []
        if len(ranked_modules) > 1:
            # Look at next 4 best modules (ranks 2-5)
            next_best_modules = ranked_modules[1:5]
            
            all_alternatives = []
            for mod in next_best_modules:
                if mod.score > 0:
                    # Get top genes from this alternative module
                    mod_genes = self.scorer.rank_genes(mod.module_id, observed_ids)
                    all_alternatives.extend(mod_genes)
            
            # Sort all alternatives by support score
            # We want the highest scoring genes globally, outside the best module
            all_alternatives.sort(key=lambda x: (x.support_score, x.stability_score), reverse=True)
            
            # Take top 3
            alternative_genes = all_alternatives[:3]
        # ------------------------------------------------
        
        # --- Discriminative Questions (Safeguard) ---
        discriminative = self.predictor.get_discriminative_questions(
            ranked_modules, observed_ids, excluded_ids, top_n=20
        )
        # --------------------------------------------
            
        # Build explanation
        explanation = []
        if best_module and best_module.score > 0:
            explanation = self.predictor.explain_scoring(
                best_module.module_id, observed_ids, excluded_ids
            )
        
        return QueryResult(
            matched_modules=ranked_modules,
            best_module=best_module,
            candidate_genes=candidate_genes,
            alternative_genes=alternative_genes,
            predicted_phenotypes=predicted,
            discriminative_questions=discriminative,
            explanation=explanation,
            observed_phenotypes=list(observed_ids),
            excluded_phenotypes=list(excluded_ids),
            unmatched_inputs=all_unmatched
        )
    
    def query_gene(self, gene: str, top_phenotypes: int = 20) -> Optional[GeneQueryResult]:
        """
        Query with a single gene to understand its phenotype context.
        
        Args:
            gene: Gene symbol
            top_phenotypes: Number of characteristic phenotypes to return
        
        Returns:
            GeneQueryResult or None if gene not found
        """
        gene_info = self.loader.get_gene_info(gene)
        if gene_info is None:
            return None
        
        module_id = gene_info.module_id
        
        # Get other genes in module
        module_genes = []
        for g in self.loader.get_module_genes(module_id):
            g_info = self.loader.get_gene_info(g)
            if g_info and g != gene:
                module_genes.append(GeneCandidate(
                    gene=g,
                    module_id=module_id,
                    support_score=0.0,
                    stability_score=g_info.stability_score,
                    classification=g_info.classification,
                    supporting_phenotypes=[]
                ))
        
        # Sort by stability
        module_genes.sort(key=lambda x: x.stability_score, reverse=True)
        
        # Get characteristic phenotypes
        char_phenotypes = self.predictor.get_expected_phenotypes(module_id, top_phenotypes)
        
        # Convert to PhenotypeInfo
        pheno_infos = []
        profile = self.loader.get_module_profile(module_id)
        if profile:
            for pred in char_phenotypes:
                pheno = profile.phenotypes.get(pred.hpo_id)
                if pheno:
                    pheno_infos.append(PhenotypeInfo(
                        name=pheno.name,
                        hpo_id=pred.hpo_id,
                        prevalence=pheno.prevalence,
                        specificity=pheno.specificity,
                        genes=pheno.genes_with
                    ))
        
        return GeneQueryResult(
            gene=gene,
            module_id=module_id,
            stability_score=gene_info.stability_score,
            classification=gene_info.classification,
            module_genes=module_genes,
            characteristic_phenotypes=pheno_infos
        )
    
    def suggest_next_phenotype(self,
                               observed: Optional[List[str]] = None,
                               excluded: Optional[List[str]] = None) -> Optional[PhenotypePrediction]:
        """
        Suggest the most informative phenotype to ask about next.
        
        Args:
            observed: Currently observed phenotypes
            excluded: Currently excluded phenotypes
        
        Returns:
            PhenotypePrediction for the most discriminative phenotype
        """
        observed = observed or []
        excluded = excluded or []
        
        observed_ids, _ = self._resolve_phenotypes(observed)
        excluded_ids, _ = self._resolve_phenotypes(excluded)
        
        ranked = self.scorer.rank_modules(observed_ids, excluded_ids)
        
        return self.predictor.suggest_next_question(ranked, observed_ids, excluded_ids)
    
    def get_module_summary(self, module_id: int) -> dict:
        """Get summary information about a module."""
        profile = self.loader.get_module_profile(module_id)
        if profile is None:
            return {}
        
        genes = list(self.loader.get_module_genes(module_id))
        gene_infos = [self.loader.get_gene_info(g) for g in genes]
        core_genes = [g.gene for g in gene_infos if g and g.classification == 'core']
        
        top_phenotypes = self.predictor.get_expected_phenotypes(module_id, top_n=10)
        
        return {
            "module_id": module_id,
            "total_genes": len(genes),
            "core_genes": core_genes,
            "top_phenotypes": [
                {"name": p.name, "prevalence": p.prevalence, "specificity": p.specificity}
                for p in top_phenotypes
            ]
        }


# Convenience function
def create_engine(data_dir: Optional[Path] = None) -> ClinicalSupportEngine:
    """Create a configured clinical support engine."""
    return ClinicalSupportEngine(data_dir)


if __name__ == "__main__":
    # Demo usage
    print("=" * 70)
    print("CLINICAL DECISION SUPPORT FRAMEWORK - Demo")
    print("=" * 70)
    
    engine = ClinicalSupportEngine()
    
    # Demo 1: BBS-like phenotype query
    print("\n[DEMO 1] Phenotype Query: Obesity + Polydactyly + Rod-cone dystrophy")
    print("-" * 70)
    
    result = engine.query(
        observed=["Obesity", "Polydactyly", "Rod-cone dystrophy"],
        excluded=[]
    )
    
    print(result.summary())
    
    # Demo 2: Gene query
    print("\n" + "=" * 70)
    print("[DEMO 2] Gene Query: RPGR")
    print("-" * 70)
    
    gene_result = engine.query_gene("RPGR")
    if gene_result:
        print(gene_result.summary())
    
    # Demo 3: Usher-like query
    print("\n" + "=" * 70)
    print("[DEMO 3] Phenotype Query: Sensorineural hearing impairment + Rod-cone dystrophy")
    print("-" * 70)
    
    result2 = engine.query(
        observed=["Sensorineural hearing impairment", "Rod-cone dystrophy"],
        excluded=[]
    )
    
    print(result2.summary())
    
    # Demo 4: Next question suggestion
    print("\n" + "=" * 70)
    print("[DEMO 4] Next Question Suggestion")
    print("-" * 70)
    
    next_q = engine.suggest_next_phenotype(
        observed=["Rod-cone dystrophy"]
    )
    if next_q:
        print(f"Suggested next question: {next_q.name}")
        print(f"  HPO ID: {next_q.hpo_id}")
        print(f"  Reason: {next_q.reason}")
