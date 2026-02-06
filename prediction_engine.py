"""
Prediction Engine for the Clinical Decision Support Framework.

Implements:
- Missing phenotype prediction
- Expected phenotype suggestions
- Score explainability
"""

from typing import List, Set, Optional
from dataclasses import dataclass

from data_loader import DataLoader, PhenotypeProfile, ModuleProfile
from output_models import PhenotypePrediction, ExplainabilityItem


@dataclass
class PredictionConfig:
    """Configuration for prediction thresholds."""
    min_prevalence: float = 20.0  # Minimum prevalence % to predict
    min_specificity: float = 5.0  # Minimum specificity % to consider
    max_predictions: int = 10  # Maximum number of predictions to return


class PredictionEngine:
    """Predicts missing phenotypes and provides explanations."""
    
    def __init__(self, loader: DataLoader, config: Optional[PredictionConfig] = None):
        """Initialize with data loader and optional config."""
        self.loader = loader
        self.config = config or PredictionConfig()
    
    def predict_missing_phenotypes(self,
                                   module_id: int,
                                   observed_hpo_ids: Set[str],
                                   excluded_hpo_ids: Set[str]) -> List[PhenotypePrediction]:
        """
        Predict phenotypes that are likely missing from the input.
        
        Returns high-prevalence phenotypes in the matched module that are:
        - Not in observed set
        - Not in excluded set
        - Have high prevalence and/or specificity
        """
        profile = self.loader.get_module_profile(module_id)
        if profile is None:
            return []
        
        candidates = []
        already_mentioned = observed_hpo_ids | excluded_hpo_ids
        
        for hpo_id, pheno in profile.phenotypes.items():
            # Skip if already mentioned
            if hpo_id in already_mentioned:
                continue
            
            # Filter by thresholds
            if pheno.prevalence < self.config.min_prevalence:
                continue
            if pheno.specificity < self.config.min_specificity:
                continue
            
            # Calculate prediction score (higher = more important to ask about)
            prediction_score = pheno.prevalence * pheno.specificity / 100.0
            
            reason = self._build_reason(pheno)
            
            candidates.append((prediction_score, PhenotypePrediction(
                name=pheno.name,
                hpo_id=hpo_id,
                prevalence=pheno.prevalence,
                specificity=pheno.specificity,
                reason=reason
            )))
        
        # Sort by prediction score descending
        candidates.sort(key=lambda x: x[0], reverse=True)
        
        # Return top N predictions
        return [pred for _, pred in candidates[:self.config.max_predictions]]
    
    def _build_reason(self, pheno: PhenotypeProfile) -> str:
        """Build human-readable reason for prediction."""
        if pheno.prevalence >= 80:
            return f"Very common in module ({pheno.prevalence:.0f}% of genes)"
        elif pheno.prevalence >= 50:
            return f"Common in module ({pheno.prevalence:.0f}% of genes)"
        elif pheno.specificity >= 50:
            return f"Highly specific to module ({pheno.specificity:.0f}% of genes with phenotype are in module)"
        else:
            return f"Characteristic (prevalence: {pheno.prevalence:.0f}%, specificity: {pheno.specificity:.0f}%)"
    
    def get_expected_phenotypes(self,
                                module_id: int,
                                top_n: int = 20) -> List[PhenotypePrediction]:
        """
        Get the most characteristic phenotypes for a module.
        
        Returns phenotypes sorted by (prevalence * specificity).
        """
        profile = self.loader.get_module_profile(module_id)
        if profile is None:
            return []
        
        candidates = []
        for hpo_id, pheno in profile.phenotypes.items():
            # Score by prevalence * specificity
            score = pheno.prevalence * pheno.specificity / 100.0
            
            candidates.append((score, PhenotypePrediction(
                name=pheno.name,
                hpo_id=hpo_id,
                prevalence=pheno.prevalence,
                specificity=pheno.specificity,
                reason=self._build_reason(pheno)
            )))
        
        candidates.sort(key=lambda x: x[0], reverse=True)
        return [pred for _, pred in candidates[:top_n]]
    
    def suggest_next_question(self,
                              ranked_modules: List,  # List of ModuleMatch
                              observed_hpo_ids: Set[str],
                              excluded_hpo_ids: Set[str]) -> Optional[PhenotypePrediction]:
        """
        Suggest the most informative phenotype to ask about next.
        
        This finds phenotypes that best discriminate between top candidates.
        """
        if len(ranked_modules) < 2:
            return None
        
        top_module = ranked_modules[0]
        second_module = ranked_modules[1]
        
        top_profile = self.loader.get_module_profile(top_module.module_id)
        second_profile = self.loader.get_module_profile(second_module.module_id)
        
        if top_profile is None or second_profile is None:
            return None
        
        already_asked = observed_hpo_ids | excluded_hpo_ids
        
        # Find phenotypes that discriminate between top two modules
        best_discriminator = None
        best_discrimination_score = 0.0
        
        for hpo_id, pheno in top_profile.phenotypes.items():
            if hpo_id in already_asked:
                continue
            
            # Get prevalence in top module
            top_prevalence = pheno.prevalence
            
            # Get prevalence in second module
            second_pheno = second_profile.phenotypes.get(hpo_id)
            second_prevalence = second_pheno.prevalence if second_pheno else 0.0
            
            # Discrimination score = difference in prevalence
            discrimination = abs(top_prevalence - second_prevalence)
            
            # Weight by average prevalence (prefer common phenotypes)
            avg_prevalence = (top_prevalence + second_prevalence) / 2.0
            if avg_prevalence < 10:  # Ignore very rare phenotypes
                continue
            
            weighted_score = discrimination * (avg_prevalence / 100.0)
            
            if weighted_score > best_discrimination_score:
                best_discrimination_score = weighted_score
                direction = "more common in best match" if top_prevalence > second_prevalence else "less common in best match"
                best_discriminator = PhenotypePrediction(
                    name=pheno.name,
                    hpo_id=hpo_id,
                    prevalence=top_prevalence,
                    specificity=pheno.specificity,
                    reason=f"Discriminating phenotype ({direction})"
                )
        
        return best_discriminator
    
    def get_discriminative_questions(self,
                                   ranked_modules: List,
                                   observed_hpo_ids: Set[str],
                                   excluded_hpo_ids: Set[str],
                                   top_n: int = 5) -> List[PhenotypePrediction]:
        """
        Get the top phenotypes to discriminate between the top 2 candidates.
        
        Returns:
            List of PhenotypePrediction objects sorted by discriminatory power.
        """
        if len(ranked_modules) < 2:
            return []
        
        top_module = ranked_modules[0]
        second_module = ranked_modules[1]
        
        # Only suggest if the second best is somewhat competitive
        # e.g., score > 0
        if second_module.score <= 0:
            return []
        
        top_profile = self.loader.get_module_profile(top_module.module_id)
        second_profile = self.loader.get_module_profile(second_module.module_id)
        
        if top_profile is None or second_profile is None:
            return []
        
        already_asked = observed_hpo_ids | excluded_hpo_ids
        
        candidates = []
        
        for hpo_id, pheno in top_profile.phenotypes.items():
            if hpo_id in already_asked:
                continue
            
            # Get prevalence in top module
            top_prevalence = pheno.prevalence
            
            # Get prevalence in second module
            second_pheno = second_profile.phenotypes.get(hpo_id)
            second_prevalence = second_pheno.prevalence if second_pheno else 0.0
            
            # Discrimination score = difference in prevalence
            discrimination = abs(top_prevalence - second_prevalence)
            
            # Weight by average prevalence (prefer common phenotypes)
            avg_prevalence = (top_prevalence + second_prevalence) / 2.0
            if avg_prevalence < 10:  # Ignore very rare phenotypes
                continue
            
            weighted_score = discrimination * (avg_prevalence / 100.0)
            
            direction = "more" if top_prevalence > second_prevalence else "less"
            reason = f"{top_prevalence:.0f}% of Module {top_module.module_id} VS {second_prevalence:.0f}% in Module {second_module.module_id}"
            
            candidates.append((weighted_score, PhenotypePrediction(
                name=pheno.name,
                hpo_id=hpo_id,
                prevalence=top_prevalence,
                specificity=pheno.specificity,
                reason=reason
            )))
        
        # Sort by weighted score descending
        candidates.sort(key=lambda x: x[0], reverse=True)
        
        # Return top N
        return [pred for _, pred in candidates[:top_n]]
    
    def explain_scoring(self,
                        module_id: int,
                        observed_hpo_ids: Set[str],
                        excluded_hpo_ids: Set[str]) -> List[ExplainabilityItem]:
        """
        Provide detailed explanation for why a module scored as it did.
        """
        profile = self.loader.get_module_profile(module_id)
        if profile is None:
            return []
        
        explanations = []
        
        # Explain observed phenotypes
        for hpo_id in observed_hpo_ids:
            pheno = profile.phenotypes.get(hpo_id)
            if pheno:
                score = (pheno.prevalence + pheno.specificity) / 200.0  # Normalized
                explanations.append(ExplainabilityItem(
                    phenotype_name=pheno.name,
                    hpo_id=hpo_id,
                    contribution=score,
                    explanation=f"Supports module (prevalence: {pheno.prevalence:.1f}%, "
                               f"specificity: {pheno.specificity:.1f}%)"
                ))
            else:
                explanations.append(ExplainabilityItem(
                    phenotype_name="Unknown",
                    hpo_id=hpo_id,
                    contribution=0.0,
                    explanation="Phenotype not present in module profile"
                ))
        
        # Explain excluded phenotypes
        for hpo_id in excluded_hpo_ids:
            pheno = profile.phenotypes.get(hpo_id)
            if pheno and pheno.prevalence > 0:
                penalty = (pheno.prevalence + pheno.specificity) / 400.0  # Half weight
                explanations.append(ExplainabilityItem(
                    phenotype_name=pheno.name,
                    hpo_id=hpo_id,
                    contribution=-penalty,
                    explanation=f"Penalizes module (excluded but {pheno.prevalence:.1f}% prevalence)"
                ))
        
        # Sort by absolute contribution
        explanations.sort(key=lambda x: abs(x.contribution), reverse=True)
        
        return explanations


if __name__ == "__main__":
    from data_loader import DataLoader
    
    loader = DataLoader()
    loader.load()
    
    engine = PredictionEngine(loader)
    
    # Test with module 7 (BBS module)
    print("Expected phenotypes for Module 7 (BBS):")
    expected = engine.get_expected_phenotypes(7, top_n=10)
    for pheno in expected:
        print(f"  {pheno.name}: {pheno.prevalence:.1f}% prevalence, {pheno.specificity:.1f}% specificity")
    
    print("\n" + "="*60)
    
    # Test missing prediction
    print("\nPredicted missing phenotypes for Module 7 (assuming some observed):")
    observed = set()
    # Add some phenotypes we "observed"
    for pheno in expected[:3]:
        observed.add(pheno.hpo_id)
    
    missing = engine.predict_missing_phenotypes(7, observed, set())
    for pheno in missing[:5]:
        print(f"  {pheno.name}: {pheno.reason}")
