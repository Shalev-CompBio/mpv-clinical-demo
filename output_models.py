"""
Output data models for the Clinical Decision Support Framework.

Defines structured result types for module matching, gene ranking,
and phenotype prediction.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict


@dataclass
class PhenotypeInfo:
    """Information about a single phenotype."""
    name: str
    hpo_id: str
    prevalence: float  # Percentage of genes in module with this phenotype
    specificity: float  # Percentage of all genes with this phenotype that are in module
    genes: List[str] = field(default_factory=list)  # Genes in module with this phenotype


@dataclass
class ModuleMatch:
    """Result of module scoring."""
    module_id: int
    score: float
    confidence: float  # 0-1, based on score separation from next best
    gene_count: int
    contributing_phenotypes: List[PhenotypeInfo] = field(default_factory=list)
    penalized_phenotypes: List[PhenotypeInfo] = field(default_factory=list)


@dataclass
class GeneCandidate:
    """A prioritized gene candidate within a module."""
    gene: str
    module_id: int
    support_score: float  # Based on phenotype overlap
    stability_score: float
    classification: str  # core, peripheral, unstable
    supporting_phenotypes: List[str] = field(default_factory=list)


@dataclass
class PhenotypePrediction:
    """A predicted missing or expected phenotype."""
    name: str
    hpo_id: str
    prevalence: float
    specificity: float
    reason: str  # e.g., "High prevalence in matched module"


@dataclass
class ExplainabilityItem:
    """Explanation for a scoring decision."""
    phenotype_name: str
    hpo_id: str
    contribution: float  # Positive or negative contribution to score
    explanation: str


@dataclass
class QueryResult:
    """Complete result from a phenotype query."""
    # Matched modules (ranked by score)
    matched_modules: List[ModuleMatch] = field(default_factory=list)
    
    # Best module details
    best_module: Optional[ModuleMatch] = None
    
    # Candidate genes within best module
    candidate_genes: List[GeneCandidate] = field(default_factory=list)
    
    # Alternative top candidates from other modules (safeguard)
    alternative_genes: List[GeneCandidate] = field(default_factory=list)
    
    # Predicted missing phenotypes
    predicted_phenotypes: List[PhenotypePrediction] = field(default_factory=list)
    
    # Discriminative questions (Top 2 module differentiation)
    discriminative_questions: List[PhenotypePrediction] = field(default_factory=list)
    
    # Score explanation
    explanation: List[ExplainabilityItem] = field(default_factory=list)
    
    # Input summary
    observed_phenotypes: List[str] = field(default_factory=list)
    excluded_phenotypes: List[str] = field(default_factory=list)
    unmatched_inputs: List[str] = field(default_factory=list)
    
    def summary(self) -> str:
        """Return a human-readable summary of results."""
        lines = []
        if self.best_module:
            lines.append(f"Best Module: {self.best_module.module_id} "
                        f"(score: {self.best_module.score:.3f}, "
                        f"confidence: {self.best_module.confidence:.2%})")
            lines.append(f"  Genes in module: {self.best_module.gene_count}")
        
        if self.candidate_genes:
            lines.append(f"\nTop Candidate Genes ({len(self.candidate_genes)} total):")
            for gene in self.candidate_genes[:5]:
                lines.append(f"  - {gene.gene} ({gene.classification}, "
                           f"support: {gene.support_score:.3f})")
        
        if self.predicted_phenotypes:
            lines.append(f"\nPredicted Missing Phenotypes:")
            for pheno in self.predicted_phenotypes[:5]:
                lines.append(f"  - {pheno.name} (prevalence: {pheno.prevalence:.1f}%)")
        
        if self.unmatched_inputs:
            lines.append(f"\nUnmatched inputs: {', '.join(self.unmatched_inputs)}")
        
        return "\n".join(lines)


@dataclass
class GeneQueryResult:
    """Result from a single-gene query."""
    gene: str
    module_id: int
    stability_score: float
    classification: str
    
    # Other genes in the same module
    module_genes: List[GeneCandidate] = field(default_factory=list)
    
    # Characteristic phenotypes for this module
    characteristic_phenotypes: List[PhenotypeInfo] = field(default_factory=list)
    
    def summary(self) -> str:
        """Return a human-readable summary."""
        lines = [
            f"Gene: {self.gene}",
            f"Module: {self.module_id}",
            f"Classification: {self.classification} (stability: {self.stability_score:.3f})",
        ]
        
        if self.module_genes:
            core_genes = [g.gene for g in self.module_genes if g.classification == 'core'][:5]
            if core_genes:
                lines.append(f"\nCore genes in module: {', '.join(core_genes)}")
        
        if self.characteristic_phenotypes:
            lines.append(f"\nCharacteristic phenotypes:")
            for pheno in self.characteristic_phenotypes[:5]:
                lines.append(f"  - {pheno.name} ({pheno.prevalence:.1f}%)")
        
        return "\n".join(lines)
