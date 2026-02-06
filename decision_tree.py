"""
Interactive Decision Tree for Clinical Decision Support.

Implements a soft decision-tree mode where users can answer Yes/No/Unknown
questions about phenotypes, with scores updated incrementally.
"""

from typing import List, Set, Optional, Dict, Tuple
from dataclasses import dataclass, field
from enum import Enum

from data_loader import DataLoader
from scoring_engine import ScoringEngine, ScoringConfig
from prediction_engine import PredictionEngine, PredictionConfig
from output_models import (
    ModuleMatch, GeneCandidate, PhenotypePrediction, QueryResult
)


class Response(Enum):
    """User response to a phenotype question."""
    YES = "yes"
    NO = "no"
    UNKNOWN = "unknown"


@dataclass
class SessionState:
    """Current state of an interactive session."""
    observed: Set[str] = field(default_factory=set)  # HPO IDs marked as present
    excluded: Set[str] = field(default_factory=set)  # HPO IDs marked as absent
    unknown: Set[str] = field(default_factory=set)  # HPO IDs marked as unknown
    history: List[Tuple[str, str, Response]] = field(default_factory=list)  # (name, hpo_id, response)
    
    @property
    def asked(self) -> Set[str]:
        """All phenotypes that have been asked about."""
        return self.observed | self.excluded | self.unknown


class InteractiveSession:
    """
    Manages an interactive clinical decision support session.
    
    Usage:
        session = InteractiveSession()
        
        # Get first question
        question = session.get_next_question()
        print(f"Is {question.name} present? (y/n/u)")
        
        # Record answer and get next question
        session.answer(question.hpo_id, Response.YES)
        next_question = session.get_next_question()
        
        # Get current best result at any time
        result = session.get_current_result()
    """
    
    def __init__(self,
                 loader: Optional[DataLoader] = None,
                 scoring_config: Optional[ScoringConfig] = None,
                 prediction_config: Optional[PredictionConfig] = None):
        """Initialize interactive session."""
        if loader is None:
            loader = DataLoader()
            loader.load()
        
        self.loader = loader
        self.scorer = ScoringEngine(loader, scoring_config)
        self.predictor = PredictionEngine(loader, prediction_config)
        self.state = SessionState()
        
        # Cache current rankings
        self._cached_modules: List[ModuleMatch] = []
        self._cache_valid = False
    
    def reset(self) -> None:
        """Reset session to initial state."""
        self.state = SessionState()
        self._cache_valid = False
    
    def _update_cache(self) -> None:
        """Update cached module rankings."""
        if not self._cache_valid:
            self._cached_modules = self.scorer.rank_modules(
                self.state.observed, 
                self.state.excluded
            )
            self._cache_valid = True
    
    def answer(self, 
               identifier: str, 
               response: Response,
               phenotype_name: Optional[str] = None) -> None:
        """
        Record an answer to a phenotype question.
        
        Args:
            identifier: HPO ID or phenotype name
            response: YES, NO, or UNKNOWN
            phenotype_name: Optional name for history tracking
        """
        # Resolve to HPO ID if needed
        hpo_id = self.loader.resolve_phenotype(identifier)
        if hpo_id is None:
            hpo_id = identifier  # Keep as-is if not resolved
        
        # Remove from any previous category
        self.state.observed.discard(hpo_id)
        self.state.excluded.discard(hpo_id)
        self.state.unknown.discard(hpo_id)
        
        # Add to appropriate category
        if response == Response.YES:
            self.state.observed.add(hpo_id)
        elif response == Response.NO:
            self.state.excluded.add(hpo_id)
        else:
            self.state.unknown.add(hpo_id)
        
        # Record in history
        name = phenotype_name or identifier
        self.state.history.append((name, hpo_id, response))
        
        # Invalidate cache
        self._cache_valid = False
    
    def answer_yes(self, identifier: str) -> None:
        """Shorthand for answering YES."""
        self.answer(identifier, Response.YES)
    
    def answer_no(self, identifier: str) -> None:
        """Shorthand for answering NO."""
        self.answer(identifier, Response.NO)
    
    def answer_unknown(self, identifier: str) -> None:
        """Shorthand for answering UNKNOWN."""
        self.answer(identifier, Response.UNKNOWN)
    
    def get_ranked_modules(self) -> List[ModuleMatch]:
        """Get current module rankings."""
        self._update_cache()
        return self._cached_modules
    
    def get_best_module(self) -> Optional[ModuleMatch]:
        """Get the current best-matching module."""
        self._update_cache()
        return self._cached_modules[0] if self._cached_modules else None
    
    def get_candidate_genes(self, top_n: int = 20) -> List[GeneCandidate]:
        """Get candidate genes from the best-matching module."""
        best = self.get_best_module()
        if best is None or best.score <= 0:
            return []
        
        return self.scorer.rank_genes(best.module_id, self.state.observed)[:top_n]
    
    def get_next_question(self) -> Optional[PhenotypePrediction]:
        """
        Get the next most informative phenotype to ask about.
        
        Uses the predictor's discrimination logic to find phenotypes
        that best separate top candidate modules.
        """
        self._update_cache()
        
        if not self._cached_modules:
            return None
        
        return self.predictor.suggest_next_question(
            self._cached_modules,
            self.state.observed,
            self.state.excluded | self.state.unknown  # Treat unknown as "don't ask again"
        )
    
    def get_predicted_phenotypes(self, top_n: int = 10) -> List[PhenotypePrediction]:
        """Get predicted missing phenotypes based on current best module."""
        best = self.get_best_module()
        if best is None or best.score <= 0:
            return []
        
        return self.predictor.predict_missing_phenotypes(
            best.module_id,
            self.state.observed,
            self.state.excluded
        )[:top_n]
    
    def get_current_result(self) -> QueryResult:
        """Get full result based on current session state."""
        self._update_cache()
        
        best = self.get_best_module()
        candidates = self.get_candidate_genes() if best else []
        predictions = self.get_predicted_phenotypes() if best else []
        
        explanation = []
        if best and best.score > 0:
            explanation = self.predictor.explain_scoring(
                best.module_id,
                self.state.observed,
                self.state.excluded
            )
        
        return QueryResult(
            matched_modules=self._cached_modules,
            best_module=best,
            candidate_genes=candidates,
            predicted_phenotypes=predictions,
            explanation=explanation,
            observed_phenotypes=list(self.state.observed),
            excluded_phenotypes=list(self.state.excluded),
            unmatched_inputs=[]
        )
    
    def get_session_summary(self) -> str:
        """Get a summary of the current session."""
        lines = [
            f"Questions answered: {len(self.state.history)}",
            f"  Observed (YES): {len(self.state.observed)}",
            f"  Excluded (NO): {len(self.state.excluded)}",
            f"  Unknown: {len(self.state.unknown)}",
        ]
        
        best = self.get_best_module()
        if best:
            lines.append(f"\nCurrent best match: Module {best.module_id}")
            lines.append(f"  Score: {best.score:.3f}")
            lines.append(f"  Confidence: {best.confidence:.1%}")
        
        return "\n".join(lines)


def run_interactive_demo():
    """Run an interactive demo session."""
    print("=" * 70)
    print("INTERACTIVE CLINICAL DECISION SUPPORT - Demo")
    print("=" * 70)
    print("\nThis demo simulates an interactive Q&A session.")
    print("Answers: Y = Yes (observed), N = No (excluded), U = Unknown\n")
    
    session = InteractiveSession()
    
    # Simulate a session with BBS-like phenotypes
    demo_answers = [
        ("Rod-cone dystrophy", Response.YES),
        ("Obesity", Response.YES),
        ("Hearing impairment", Response.NO),
        ("Polydactyly", Response.YES),
    ]
    
    for pheno_name, response in demo_answers:
        hpo_id = session.loader.resolve_phenotype(pheno_name)
        if hpo_id:
            print(f"Q: Is '{pheno_name}' present?")
            print(f"A: {response.value.upper()}")
            session.answer(hpo_id, response, pheno_name)
            
            best = session.get_best_module()
            if best:
                print(f"   -> Current best: Module {best.module_id} "
                      f"(score: {best.score:.3f}, conf: {best.confidence:.1%})")
            print()
    
    print("-" * 70)
    print("FINAL RESULT")
    print("-" * 70)
    print(session.get_session_summary())
    
    print("\nTop candidate genes:")
    for gene in session.get_candidate_genes()[:5]:
        print(f"  {gene.gene} ({gene.classification})")
    
    next_q = session.get_next_question()
    if next_q:
        print(f"\nSuggested next question: {next_q.name}")


if __name__ == "__main__":
    run_interactive_demo()
