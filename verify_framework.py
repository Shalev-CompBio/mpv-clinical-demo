"""
Verification script for the Clinical Decision Support Framework.

Runs comprehensive tests without requiring pytest.
Run with: python verify_framework.py
"""

import sys


def run_test(name, test_func):
    """Run a single test and report result."""
    try:
        test_func()
        print(f"  [PASS] {name}")
        return True
    except AssertionError as e:
        print(f"  [FAIL] {name}: {e}")
        return False
    except Exception as e:
        print(f"  [ERROR] {name}: {type(e).__name__}: {e}")
        return False


def main():
    """Run all verification tests."""
    print("=" * 70)
    print("CLINICAL DECISION SUPPORT FRAMEWORK - Verification")
    print("=" * 70)
    
    # Import modules
    print("\n[1] Loading modules...")
    try:
        from data_loader import DataLoader
        from scoring_engine import ScoringEngine
        from prediction_engine import PredictionEngine
        from clinical_support import ClinicalSupportEngine
        from decision_tree import InteractiveSession, Response
        print("  All modules imported successfully")
    except Exception as e:
        print(f"  FAILED to import modules: {e}")
        return False
    
    # Load data
    print("\n[2] Loading data...")
    loader = DataLoader()
    loader.load()
    
    engine = ClinicalSupportEngine()
    scorer = ScoringEngine(loader)
    predictor = PredictionEngine(loader)
    
    passed = 0
    failed = 0
    
    # Test Suite 1: Data Loading
    print("\n[3] Data Loading Tests")
    print("-" * 40)
    
    def test_gene_count():
        assert len(loader.gene_mapping) >= 400, f"Expected 400+ genes, got {len(loader.gene_mapping)}"
    
    def test_module_count():
        assert len(loader.module_profiles) == 13, f"Expected 13 modules, got {len(loader.module_profiles)}"
    
    def test_phenotype_count():
        assert len(loader.phenotype_to_modules) >= 1000, f"Expected 1000+ phenotypes"
    
    def test_rpgr_in_module12():
        info = loader.get_gene_info("RPGR")
        assert info is not None, "RPGR not found"
        assert info.module_id == 12, f"RPGR should be in module 12, got {info.module_id}"
    
    for name, func in [
        ("Load 400+ genes", test_gene_count),
        ("Load 13 modules", test_module_count),
        ("Index 1000+ phenotypes", test_phenotype_count),
        ("RPGR in module 12", test_rpgr_in_module12),
    ]:
        if run_test(name, func):
            passed += 1
        else:
            failed += 1
    
    # Test Suite 2: Module Scoring
    print("\n[4] Module Scoring Tests")
    print("-" * 40)
    
    def test_bbs_module_match():
        """BBS phenotypes should match Module 7"""
        observed = set()
        for name in ["Obesity", "Polydactyly", "Rod-cone dystrophy"]:
            hpo_id = loader.resolve_phenotype(name)
            if hpo_id:
                observed.add(hpo_id)
        
        matches = scorer.rank_modules(observed, set())
        assert matches[0].module_id == 7, f"Expected module 7, got {matches[0].module_id}"
    
    def test_usher_module_match():
        """Usher phenotypes should match Module 11"""
        observed = set()
        for name in ["Sensorineural hearing impairment", "Rod-cone dystrophy"]:
            hpo_id = loader.resolve_phenotype(name)
            if hpo_id:
                observed.add(hpo_id)
        
        matches = scorer.rank_modules(observed, set())
        assert matches[0].module_id == 11, f"Expected module 11, got {matches[0].module_id}"
    
    def test_bbs_genes_ranked():
        """BBS genes should be top-ranked in Module 7"""
        observed = set()
        for name in ["Obesity", "Polydactyly"]:
            hpo_id = loader.resolve_phenotype(name)
            if hpo_id:
                observed.add(hpo_id)
        
        genes = scorer.rank_genes(7, observed)
        top_names = [g.gene for g in genes[:10]]
        bbs_genes = ["BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BBS10", "BBS12"]
        overlap = set(top_names) & set(bbs_genes)
        assert len(overlap) >= 3, f"Expected at least 3 BBS genes in top 10, got {overlap}"
    
    for name, func in [
        ("BBS phenotypes -> Module 7", test_bbs_module_match),
        ("Usher phenotypes -> Module 11", test_usher_module_match),
        ("BBS genes top-ranked", test_bbs_genes_ranked),
    ]:
        if run_test(name, func):
            passed += 1
        else:
            failed += 1
    
    # Test Suite 3: Clinical API
    print("\n[5] Clinical API Tests")
    print("-" * 40)
    
    def test_query_basic():
        result = engine.query(observed=["Rod-cone dystrophy"])
        assert result.best_module is not None
        assert len(result.candidate_genes) > 0
    
    def test_gene_query():
        result = engine.query_gene("RPGR")
        assert result is not None
        assert result.module_id == 12
        assert len(result.characteristic_phenotypes) > 0
    
    def test_empty_query():
        result = engine.query(observed=[], excluded=[])
        assert result is not None
        assert result.best_module.score == 0
    
    def test_unmatched_handling():
        result = engine.query(observed=["FakeNotRealPhenotype"])
        assert "FakeNotRealPhenotype" in result.unmatched_inputs
    
    for name, func in [
        ("Basic phenotype query", test_query_basic),
        ("Gene query (RPGR)", test_gene_query),
        ("Empty query handling", test_empty_query),
        ("Unmatched phenotype handling", test_unmatched_handling),
    ]:
        if run_test(name, func):
            passed += 1
        else:
            failed += 1
    
    # Test Suite 4: Interactive Session
    print("\n[6] Interactive Session Tests")
    print("-" * 40)
    
    def test_session_creation():
        session = InteractiveSession(loader)
        assert len(session.state.observed) == 0
    
    def test_session_answers():
        session = InteractiveSession(loader)
        session.answer_yes("HP:0000510")
        session.answer_no("HP:0001513")
        assert "HP:0000510" in session.state.observed
        assert "HP:0001513" in session.state.excluded
    
    def test_incremental_scoring():
        session = InteractiveSession(loader)
        initial = session.get_best_module()
        assert initial is None or initial.score == 0
        
        session.answer_yes("HP:0000510")
        after = session.get_best_module()
        assert after.score > 0
    
    for name, func in [
        ("Session creation", test_session_creation),
        ("Answer recording", test_session_answers),
        ("Incremental scoring", test_incremental_scoring),
    ]:
        if run_test(name, func):
            passed += 1
        else:
            failed += 1
    
    # Summary
    print("\n" + "=" * 70)
    print(f"VERIFICATION COMPLETE: {passed} passed, {failed} failed")
    print("=" * 70)
    
    if failed == 0:
        print("\nAll tests passed! Framework is ready for use.")
        return True
    else:
        print(f"\n{failed} test(s) failed. Please review.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
