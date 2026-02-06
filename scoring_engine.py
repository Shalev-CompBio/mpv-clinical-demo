"""
Scoring Engine for the Clinical Decision Support Framework.

Implements:
- Module scoring based on phenotype prevalence and specificity
- Gene ranking within modules
- Confidence calculation based on score separation
"""

from typing import List, Dict, Set, Optional, Tuple
from dataclasses import dataclass

from data_loader import DataLoader, PhenotypeProfile, ModuleProfile, GeneInfo
from output_models import (
    ModuleMatch, GeneCandidate, PhenotypeInfo, ExplainabilityItem
)


@dataclass
class ScoringConfig:
    """Configuration for scoring weights and penalties."""
    prevalence_weight: float = 1.0  # Weight for phenotype prevalence
    specificity_weight: float = 1.0  # Weight for phenotype specificity
    exclusion_penalty: float = 0.5  # Penalty multiplier for excluded phenotypes
    stability_bonus: float = 0.1  # Bonus for core genes over unstable
    min_confidence_threshold: float = 0.1  # Minimum score separation for confidence


class ScoringEngine:
    """Scores modules and genes based on phenotype input."""
    
    def __init__(self, loader: DataLoader, config: Optional[ScoringConfig] = None):
        """Initialize with data loader and optional config."""
        self.loader = loader
        self.config = config or ScoringConfig()
    
    def score_phenotype(self, pheno: PhenotypeProfile) -> float:
        """
        Calculate phenotype contribution to module score.
        
        Score = prevalence * specificity (both normalized to 0-1)
        Higher scores for phenotypes that are common in the module AND rare elsewhere.
        """
        # Normalize percentages to 0-1 range
        prevalence = pheno.prevalence / 100.0
        specificity = pheno.specificity / 100.0
        
        return (self.config.prevalence_weight * prevalence + 
                self.config.specificity_weight * specificity)
    
    def score_module(self, 
                     module_id: int,
                     observed_hpo_ids: Set[str],
                     excluded_hpo_ids: Set[str]) -> Tuple[float, List[ExplainabilityItem]]:
        """
        Score a single module against observed and excluded phenotypes.
        
        Returns:
            Tuple of (score, explanation items)
        """
        profile = self.loader.get_module_profile(module_id)
        if profile is None:
            return 0.0, []
        
        score = 0.0
        explanations = []
        
        # Score observed phenotypes present in module
        for hpo_id in observed_hpo_ids:
            pheno = profile.phenotypes.get(hpo_id)
            if pheno is not None:
                contribution = self.score_phenotype(pheno)
                score += contribution
                explanations.append(ExplainabilityItem(
                    phenotype_name=pheno.name,
                    hpo_id=hpo_id,
                    contribution=contribution,
                    explanation=f"Present in module (prevalence: {pheno.prevalence:.1f}%, "
                               f"specificity: {pheno.specificity:.1f}%)"
                ))
        
        # Penalize excluded phenotypes present in module
        for hpo_id in excluded_hpo_ids:
            pheno = profile.phenotypes.get(hpo_id)
            if pheno is not None and pheno.prevalence > 0:
                # Only penalize if phenotype is actually observed in module genes
                penalty = self.score_phenotype(pheno) * self.config.exclusion_penalty
                score -= penalty
                explanations.append(ExplainabilityItem(
                    phenotype_name=pheno.name,
                    hpo_id=hpo_id,
                    contribution=-penalty,
                    explanation=f"Excluded but present in module (penalty)"
                ))
        
        return max(0.0, score), explanations
    
    def rank_modules(self,
                     observed_hpo_ids: Set[str],
                     excluded_hpo_ids: Set[str]) -> List[ModuleMatch]:
        """
        Score and rank all modules against the input phenotypes.
        
        Returns:
            List of ModuleMatch objects sorted by score (descending)
        """
        results = []
        
        for module_id in self.loader.get_all_module_ids():
            score, explanations = self.score_module(
                module_id, observed_hpo_ids, excluded_hpo_ids
            )
            
            profile = self.loader.get_module_profile(module_id)
            gene_count = len(profile.all_genes) if profile else 0
            
            # Collect contributing phenotypes
            contributing = []
            penalized = []
            for exp in explanations:
                hpo_id = exp.hpo_id
                pheno = profile.phenotypes.get(hpo_id) if profile else None
                if pheno:
                    info = PhenotypeInfo(
                        name=pheno.name,
                        hpo_id=hpo_id,
                        prevalence=pheno.prevalence,
                        specificity=pheno.specificity,
                        genes=pheno.genes_with
                    )
                    if exp.contribution > 0:
                        contributing.append(info)
                    else:
                        penalized.append(info)
            
            results.append(ModuleMatch(
                module_id=module_id,
                score=score,
                confidence=0.0,  # Calculated after ranking
                gene_count=gene_count,
                contributing_phenotypes=contributing,
                penalized_phenotypes=penalized
            ))
        
        # Sort by score descending
        results.sort(key=lambda x: x.score, reverse=True)
        
        # Calculate confidence based on score separation
        if len(results) >= 2 and results[0].score > 0:
            separation = (results[0].score - results[1].score) / results[0].score
            results[0].confidence = min(1.0, max(0.0, separation))
            
            # Lower confidence for other matches
            for i, match in enumerate(results[1:], 1):
                if results[0].score > 0:
                    relative_score = match.score / results[0].score
                    match.confidence = relative_score * (1 - results[0].confidence)
        elif len(results) == 1 and results[0].score > 0:
            results[0].confidence = 1.0
        
        return results
    
    def rank_genes(self,
                   module_id: int,
                   observed_hpo_ids: Set[str]) -> List[GeneCandidate]:
        """
        Rank genes within a module based on phenotype support.
        
        A gene's support score = sum of (prevalence * specificity) for 
        observed phenotypes that the gene is annotated with.
        """
        profile = self.loader.get_module_profile(module_id)
        if profile is None:
            return []
        
        gene_scores: Dict[str, Tuple[float, List[str]]] = {}  # gene -> (score, pheno list)
        
        # Calculate phenotype support for each gene
        for hpo_id in observed_hpo_ids:
            pheno = profile.phenotypes.get(hpo_id)
            if pheno is None:
                continue
            
            contribution = self.score_phenotype(pheno)
            
            for gene in pheno.genes_with:
                if gene not in gene_scores:
                    gene_scores[gene] = (0.0, [])
                current_score, pheno_list = gene_scores[gene]
                gene_scores[gene] = (current_score + contribution, 
                                     pheno_list + [pheno.name])
        
        # Add all module genes even if no phenotype support
        module_genes = self.loader.get_module_genes(module_id)
        for gene in module_genes:
            if gene not in gene_scores:
                gene_scores[gene] = (0.0, [])
        
        # Build candidate list with gene info
        candidates = []
        for gene, (support_score, supporting_phenos) in gene_scores.items():
            gene_info = self.loader.get_gene_info(gene)
            if gene_info is None:
                continue
            
            # Apply stability bonus
            final_score = support_score
            if gene_info.classification == 'core':
                final_score += self.config.stability_bonus
            elif gene_info.classification == 'unstable':
                final_score -= self.config.stability_bonus * 0.5
            
            candidates.append(GeneCandidate(
                gene=gene,
                module_id=module_id,
                support_score=final_score,
                stability_score=gene_info.stability_score,
                classification=gene_info.classification,
                supporting_phenotypes=supporting_phenos
            ))
        
        # Sort by support score descending, then by stability
        candidates.sort(key=lambda x: (x.support_score, x.stability_score), reverse=True)
        
        return candidates
    
    def get_module_explanation(self,
                               module_id: int,
                               observed_hpo_ids: Set[str],
                               excluded_hpo_ids: Set[str]) -> List[ExplainabilityItem]:
        """Get detailed explanation for module scoring."""
        _, explanations = self.score_module(module_id, observed_hpo_ids, excluded_hpo_ids)
        # Sort by absolute contribution
        explanations.sort(key=lambda x: abs(x.contribution), reverse=True)
        return explanations


if __name__ == "__main__":
    # Test scoring
    from data_loader import DataLoader
    
    loader = DataLoader()
    loader.load()
    
    engine = ScoringEngine(loader)
    
    # Test with some phenotypes that might be characteristic of BBS (module 7)
    # Using phenotype names that should be resolved
    test_observed = set()
    test_excluded = set()
    
    # Try to find obesity-related and polydactyly phenotypes
    for name in ["Obesity", "Polydactyly", "Rod-cone dystrophy", "Retinitis pigmentosa"]:
        hpo_id = loader.resolve_phenotype(name)
        if hpo_id:
            test_observed.add(hpo_id)
            print(f"Resolved '{name}' -> {hpo_id}")
    
    print(f"\nScoring with {len(test_observed)} observed phenotypes...")
    
    matches = engine.rank_modules(test_observed, test_excluded)
    
    print("\nTop 5 module matches:")
    for match in matches[:5]:
        print(f"  Module {match.module_id}: score={match.score:.3f}, "
              f"confidence={match.confidence:.2%}, genes={match.gene_count}")
    
    if matches and matches[0].score > 0:
        best = matches[0]
        print(f"\nTop genes in module {best.module_id}:")
        genes = engine.rank_genes(best.module_id, test_observed)
        for gene in genes[:10]:
            print(f"  {gene.gene}: {gene.classification}, "
                  f"support={gene.support_score:.3f}")
