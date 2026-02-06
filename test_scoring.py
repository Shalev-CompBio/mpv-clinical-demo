"""
Automated tests for the Clinical Decision Support Framework.

Tests cover:
- Data loading and indexing
- Module scoring
- Gene ranking
- Phenotype prediction
- Edge cases
"""

import pytest
from pathlib import Path

from data_loader import DataLoader
from scoring_engine import ScoringEngine, ScoringConfig
from prediction_engine import PredictionEngine
from clinical_support import ClinicalSupportEngine
from decision_tree import InteractiveSession, Response


# Fixtures
@pytest.fixture(scope="module")
def loader():
    """Load data once for all tests."""
    loader = DataLoader()
    loader.load()
    return loader


@pytest.fixture(scope="module") 
def engine(loader):
    """Create engine with loaded data."""
    return ClinicalSupportEngine()


@pytest.fixture(scope="module")
def scorer(loader):
    """Create scorer with loaded data."""
    return ScoringEngine(loader)


class TestDataLoader:
    """Tests for data_loader.py"""
    
    def test_load_genes(self, loader):
        """Test gene loading."""
        assert len(loader.gene_mapping) > 400
        assert "RPGR" in loader.gene_mapping
        assert "BBS1" in loader.gene_mapping
    
    def test_load_modules(self, loader):
        """Test module loading."""
        assert len(loader.module_profiles) == 13
        for i in range(13):
            assert i in loader.module_profiles
    
    def test_phenotype_indexing(self, loader):
        """Test phenotype index."""
        assert len(loader.phenotype_to_modules) > 1000
        assert "HP:0000510" in loader.phenotype_to_modules  # Rod-cone dystrophy
    
    def test_phenotype_resolution_by_name(self, loader):
        """Test phenotype resolution by name."""
        # Exact name match
        hpo_id = loader.resolve_phenotype("Visual impairment")
        assert hpo_id is not None
        assert hpo_id.startswith("HP:")
    
    def test_phenotype_resolution_by_id(self, loader):
        """Test phenotype resolution by HPO ID."""
        hpo_id = loader.resolve_phenotype("HP:0000510")
        assert hpo_id == "HP:0000510"
    
    def test_gene_info(self, loader):
        """Test gene info retrieval."""
        info = loader.get_gene_info("RPGR")
        assert info is not None
        assert info.module_id == 12
        assert info.classification == "core"


class TestScoringEngine:
    """Tests for scoring_engine.py"""
    
    def test_bbs_module_match(self, scorer, loader):
        """Test that BBS phenotypes match Module 7."""
        observed = set()
        for name in ["Obesity", "Polydactyly", "Rod-cone dystrophy"]:
            hpo_id = loader.resolve_phenotype(name)
            if hpo_id:
                observed.add(hpo_id)
        
        matches = scorer.rank_modules(observed, set())
        
        assert len(matches) > 0
        assert matches[0].module_id == 7  # BBS module
        assert matches[0].score > 0
    
    def test_usher_module_match(self, scorer, loader):
        """Test that Usher phenotypes match Module 11."""
        observed = set()
        for name in ["Sensorineural hearing impairment", "Rod-cone dystrophy"]:
            hpo_id = loader.resolve_phenotype(name)
            if hpo_id:
                observed.add(hpo_id)
        
        matches = scorer.rank_modules(observed, set())
        
        assert len(matches) > 0
        assert matches[0].module_id == 11  # Usher module
    
    def test_exclusion_penalty(self, scorer, loader):
        """Test that excluded phenotypes reduce score."""
        observed = {"HP:0000510"}  # Rod-cone dystrophy
        
        # Score without exclusion
        matches_no_exc = scorer.rank_modules(observed, set())
        
        # Score with exclusion (add a common phenotype as excluded)
        excluded = {"HP:0001513"}  # Obesity
        matches_with_exc = scorer.rank_modules(observed, excluded)
        
        # Module 7 should score lower with obesity excluded
        mod7_no_exc = next((m for m in matches_no_exc if m.module_id == 7), None)
        mod7_with_exc = next((m for m in matches_with_exc if m.module_id == 7), None)
        
        assert mod7_no_exc is not None
        assert mod7_with_exc is not None
        assert mod7_with_exc.score <= mod7_no_exc.score
    
    def test_gene_ranking_bbs(self, scorer, loader):
        """Test gene ranking within BBS module."""
        observed = set()
        for name in ["Obesity", "Polydactyly"]:
            hpo_id = loader.resolve_phenotype(name)
            if hpo_id:
                observed.add(hpo_id)
        
        genes = scorer.rank_genes(7, observed)
        
        assert len(genes) > 0
        # BBS genes should be top-ranked
        top_gene_names = [g.gene for g in genes[:10]]
        bbs_genes = ["BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BBS10", "BBS12"]
        
        # At least some BBS genes should be in top 10
        overlap = set(top_gene_names) & set(bbs_genes)
        assert len(overlap) >= 3
    
    def test_core_genes_ranked_higher(self, scorer):
        """Test that core genes are ranked higher than unstable."""
        genes = scorer.rank_genes(0, set())  # Module 0, no phenotype filter
        
        core_genes = [g for g in genes if g.classification == "core"]
        unstable_genes = [g for g in genes if g.classification == "unstable"]
        
        if core_genes and unstable_genes:
            # Core genes should generally score higher
            avg_core = sum(g.support_score for g in core_genes) / len(core_genes)
            avg_unstable = sum(g.support_score for g in unstable_genes) / len(unstable_genes)
            assert avg_core >= avg_unstable


class TestClinicalSupport:
    """Tests for clinical_support.py main API."""
    
    def test_phenotype_query(self, engine):
        """Test basic phenotype query."""
        result = engine.query(
            observed=["Obesity", "Rod-cone dystrophy"],
            excluded=[]
        )
        
        assert result is not None
        assert result.best_module is not None
        assert len(result.matched_modules) == 13
        assert len(result.candidate_genes) > 0
    
    def test_gene_query(self, engine):
        """Test gene query."""
        result = engine.query_gene("RPGR")
        
        assert result is not None
        assert result.gene == "RPGR"
        assert result.module_id == 12
        assert len(result.characteristic_phenotypes) > 0
    
    def test_gene_query_not_found(self, engine):
        """Test gene query with invalid gene."""
        result = engine.query_gene("NOTAREALGENE")
        assert result is None
    
    def test_empty_query(self, engine):
        """Test query with no phenotypes."""
        result = engine.query(observed=[], excluded=[])
        
        assert result is not None
        assert result.best_module is not None
        # With no input, scores should all be 0
        assert result.best_module.score == 0
    
    def test_unmatched_phenotypes(self, engine):
        """Test handling of unrecognized phenotype names."""
        result = engine.query(
            observed=["Rod-cone dystrophy", "NotARealPhenotype"],
            excluded=[]
        )
        
        assert "NotARealPhenotype" in result.unmatched_inputs
        assert len(result.observed_phenotypes) == 1


class TestInteractiveSession:
    """Tests for decision_tree.py interactive mode."""
    
    def test_session_creation(self, loader):
        """Test session creation."""
        session = InteractiveSession(loader)
        assert session is not None
        assert len(session.state.observed) == 0
    
    def test_answer_recording(self, loader):
        """Test answer recording."""
        session = InteractiveSession(loader)
        
        session.answer_yes("HP:0000510")
        assert "HP:0000510" in session.state.observed
        
        session.answer_no("HP:0001513")
        assert "HP:0001513" in session.state.excluded
        
        session.answer_unknown("HP:0000505")
        assert "HP:0000505" in session.state.unknown
    
    def test_incremental_scoring(self, loader):
        """Test that scores update incrementally."""
        session = InteractiveSession(loader)
        
        # Initial state
        initial = session.get_best_module()
        assert initial is None or initial.score == 0
        
        # Add a phenotype
        session.answer_yes("HP:0000510")  # Rod-cone dystrophy
        after_one = session.get_best_module()
        assert after_one.score > 0
        
        # Add another
        session.answer_yes("HP:0001513")  # Obesity
        after_two = session.get_best_module()
        assert after_two.score >= after_one.score
    
    def test_session_reset(self, loader):
        """Test session reset."""
        session = InteractiveSession(loader)
        session.answer_yes("HP:0000510")
        
        session.reset()
        
        assert len(session.state.observed) == 0
        assert len(session.state.history) == 0


class TestEdgeCases:
    """Edge case tests."""
    
    def test_single_phenotype(self, engine):
        """Test with single phenotype input."""
        result = engine.query(observed=["Rod-cone dystrophy"])
        
        assert result.best_module is not None
        assert len(result.candidate_genes) > 0
    
    def test_all_excluded(self, engine):
        """Test with only excluded phenotypes."""
        result = engine.query(
            observed=[],
            excluded=["Obesity", "Polydactyly"]
        )
        
        assert result is not None
        # Scores should be 0 or negative (but we floor at 0)
        assert result.best_module.score == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
