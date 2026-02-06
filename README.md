# Modular Phenotype-Driven Variant Prioritization (MPV)

## A Clinical Decision Support Framework for Inherited Retinal Diseases

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Table of Contents

1. [Overview](#overview)
2. [Conceptual Framework](#conceptual-framework)
3. [Data Foundation](#data-foundation)
4. [System Architecture](#system-architecture)
5. [Module Descriptions](#module-descriptions)
6. [Scoring Algorithm](#scoring-algorithm)
7. [Installation](#installation)
8. [Usage Guide](#usage-guide)
9. [API Reference](#api-reference)
10. [Validation Results](#validation-results)
11. [Limitations](#limitations)
12. [Citation](#citation)

---

## Overview

The **Modular Phenotype-driven Variant prioritization (MPV)** framework is a transparent, data-driven clinical decision support system designed to assist in the molecular diagnosis of **Inherited Retinal Diseases (IRDs)**.

Given a partial clinical phenotype description, the system:

1. **Identifies the most compatible disease module** from a set of 13 phenotypically-defined gene clusters
2. **Prioritizes candidate genes** within the matched module based on phenotype support
3. **Predicts missing or expected phenotypes** to guide further clinical evaluation
4. **Provides full explainability** for all scoring decisions

### Key Design Principles

- **Transparency**: All scoring logic is explicit and interpretable
- **Data-driven**: Built exclusively on pre-computed, study-specific phenotype profiles
- **No external dependencies**: Does not require real-time ontology traversal or external APIs
- **Robust to incomplete data**: Designed to function under incomplete or noisy clinical input

---

## Conceptual Framework

### Modular Disease Architecture

IRDs exhibit significant phenotypic and genetic heterogeneity. This framework operationalizes the hypothesis that IRD genes cluster into **phenotypically coherent modules** that can be leveraged for clinical inference.

Each module represents a group of genes that share:
- Common phenotypic presentations
- Related biological pathways
- Similar clinical trajectories

### From Phenotype to Gene

The clinical workflow is:

```
Clinical Phenotypes → Module Matching → Gene Prioritization → Missing Phenotype Prediction
```

Unlike machine learning classifiers, MPV is a **transparent scoring engine** where:
- Each phenotype contributes an explicit, traceable score
- Exclusions (absent phenotypes) apply explicit penalties
- Confidence is derived from score separation between competing modules

---

## Data Foundation

The framework is built on two pre-computed data files generated from upstream analysis:

### 1. Module–Phenotype Profiles

**File**: `module_top500_HPO_background_comparison_20260128_1907.xlsx`

An Excel workbook containing **13 sheets** (one per disease module), each with **500 phenotype entries**.

| Column | Description |
|--------|-------------|
| `phenotype_name` | Human-readable phenotype name |
| `hpo_id` | HPO identifier (e.g., HP:0000510) |
| `target_module_genes_with_phenotype` | Genes in module annotated with this phenotype |
| `target_module_gene_count_with_phenotype` | Count of module genes with phenotype |
| `target_module_genes_without_phenotype` | Genes in module NOT annotated with phenotype |
| `target_module_phenotype_prevalence_percent` | **Prevalence**: % of module genes with phenotype |
| `non_target_ird_genes_with_phenotype` | Non-module IRD genes with phenotype |
| `non_target_ird_gene_count_with_phenotype` | Count of non-module IRD genes |
| `non_ird_genes_with_phenotype` | Non-IRD genes with phenotype (~20,000 background) |
| `non_ird_gene_count_with_phenotype` | Count of non-IRD genes |
| `target_module_share_of_phenotype_percent` | **Specificity**: % of all genes with phenotype that are in module |

### 2. Gene–Module Mapping

**File**: `gene_classification_20251230_1900.csv`

A CSV file containing **442 IRD genes** with module assignments.

| Column | Description |
|--------|-------------|
| `gene` | HGNC gene symbol |
| `module_id` | Assigned module (0–12) |
| `stability_score` | Clustering stability score (0–1) |
| `classification` | Gene stability class: `core`, `peripheral`, or `unstable` |

### Module Summary

| Module | Primary Association | Example Genes | Gene Count |
|--------|---------------------|---------------|------------|
| 0 | Isolated RP, splicing | PRPF31, PRPF8, RP2, CHM, ABCA4 | 73 |
| 1 | Neurodegeneration | CLN3, CLN5, SPG7, MFN2 | 70 |
| 2 | Metabolic/mitochondrial | NDUFA12, MECR | 51 |
| 3 | Syndromic, peroxisomal | PEX1, PEX6, COL2A1 | 49 |
| 4 | Ciliopathies, Joubert | CEP290, CC2D2A, TMEM67 | 33 |
| 5 | Vitreoretinopathies | FZD4, NDP, COL2A1 | 35 |
| 6 | LCA, photoreceptor | RPE65, GUCY2D, CRB1, NMNAT1 | 26 |
| 7 | **Bardet-Biedl Syndrome** | BBS1-12, ALMS1, MKKS | 26 |
| 8 | Albinism, pigmentation | TYR, OCA2, HPS1-6 | 25 |
| 9 | **Mitochondrial** | MT-ND1-6, MT-ATP6, OPA1 | 30 |
| 10 | **CSNB**, phototransduction | NYX, GRM6, TRPM1, RHO | 20 |
| 11 | **Usher Syndrome** | MYO7A, USH1C, USH2A, CDH23 | 16 |
| 12 | **Achromatopsia**, cone | CNGA3, CNGB3, PDE6C, RPGR | 8 |

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        INPUT LAYER                               │
├─────────────────────────────────────────────────────────────────┤
│  Observed Phenotypes    Excluded Phenotypes    Single Gene      │
│  (present features)     (absent features)      (for lookup)     │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│                        DATA LAYER                                │
├─────────────────────────────────────────────────────────────────┤
│  data_loader.py                                                  │
│  ├── ModuleProfile: phenotype profiles per module               │
│  ├── GeneInfo: gene-to-module mappings                          │
│  └── PhenotypeIndex: name/HPO ID → module lookup                │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│                      SCORING LAYER                               │
├─────────────────────────────────────────────────────────────────┤
│  scoring_engine.py                                               │
│  ├── Module scoring: prevalence × specificity                   │
│  ├── Exclusion penalties                                         │
│  ├── Confidence calculation (score separation)                  │
│  └── Gene ranking within modules                                 │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│                     PREDICTION LAYER                             │
├─────────────────────────────────────────────────────────────────┤
│  prediction_engine.py                                            │
│  ├── Missing phenotype prediction                                │
│  ├── Next-question suggestion (discriminative phenotypes)       │
│  └── Score explainability                                        │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│                       OUTPUT LAYER                               │
├─────────────────────────────────────────────────────────────────┤
│  output_models.py                                                │
│  ├── QueryResult: full query response                            │
│  ├── ModuleMatch: module with score and confidence              │
│  ├── GeneCandidate: ranked gene with support details            │
│  └── PhenotypePrediction: predicted missing phenotypes          │
└─────────────────────────────────────────────────────────────────┘
```

---

## Module Descriptions

### File Structure

```
MPV/
├── data_loader.py          # Data loading and indexing
├── scoring_engine.py       # Module and gene scoring
├── prediction_engine.py    # Missing phenotype prediction
├── clinical_support.py     # Main API entry point
├── decision_tree.py        # Interactive Q&A mode
├── output_models.py        # Result dataclasses
├── verify_framework.py     # Standalone test suite
├── test_scoring.py         # Pytest test suite
├── classic_cases.py        # Classic IRD case demonstrations
├── challenge_tests.py      # Edge case testing
└── README.md               # This documentation
```

### Core Modules

#### `data_loader.py`

Responsible for loading and indexing all input data.

**Classes:**
- `DataLoader`: Main loader class
- `ModuleProfile`: Contains all phenotypes for a module
- `PhenotypeProfile`: Single phenotype with prevalence/specificity
- `GeneInfo`: Gene with module assignment and stability

**Key Methods:**
```python
loader = DataLoader()
loader.load()  # Load Excel and CSV files

# Resolve phenotype by name or HPO ID
hpo_id = loader.resolve_phenotype("Retinitis pigmentosa")
hpo_id = loader.resolve_phenotype("HP:0000510")

# Get gene information
gene_info = loader.get_gene_info("RPGR")
# Returns: GeneInfo(gene='RPGR', module_id=12, stability_score=0.857, classification='core')

# Get module profile
profile = loader.get_module_profile(7)  # BBS module
```

#### `scoring_engine.py`

Implements the core scoring algorithm.

**Classes:**
- `ScoringEngine`: Main scoring logic
- `ScoringConfig`: Configurable weights and penalties

**Key Methods:**
```python
scorer = ScoringEngine(loader)

# Score all modules
matches = scorer.rank_modules(observed_hpo_ids, excluded_hpo_ids)

# Rank genes within a module
genes = scorer.rank_genes(module_id=7, observed_hpo_ids)
```

#### `prediction_engine.py`

Handles phenotype prediction and explainability.

**Classes:**
- `PredictionEngine`: Prediction logic
- `PredictionConfig`: Thresholds for prediction

**Key Methods:**
```python
predictor = PredictionEngine(loader)

# Predict missing phenotypes
missing = predictor.predict_missing_phenotypes(module_id, observed, excluded)

# Get expected phenotypes for a module
expected = predictor.get_expected_phenotypes(module_id, top_n=20)

# Suggest next discriminating question
next_q = predictor.suggest_next_question(ranked_modules, observed, excluded)
```

#### `clinical_support.py`

Main API entry point combining all functionality.

**Classes:**
- `ClinicalSupportEngine`: Unified interface

**Key Methods:**
```python
engine = ClinicalSupportEngine()

# Query with phenotypes
result = engine.query(
    observed=["Obesity", "Polydactyly", "Rod-cone dystrophy"],
    excluded=["Hearing impairment"]
)

# Query with gene
gene_result = engine.query_gene("RPGR")

# Get next question suggestion
next_q = engine.suggest_next_phenotype(observed=["Rod-cone dystrophy"])
```

#### `decision_tree.py`

Interactive session management for step-by-step Q&A.

**Classes:**
- `InteractiveSession`: Manages session state
- `Response`: Enum for YES/NO/UNKNOWN

**Key Methods:**
```python
session = InteractiveSession()

session.answer_yes("Rod-cone dystrophy")
session.answer_no("Hearing impairment")
session.answer_unknown("Obesity")

# Get current best result
result = session.get_current_result()

# Get next suggested question
next_q = session.get_next_question()
```

---

## Scoring Algorithm

### Module Scoring Formula

For each module, the score is computed as:

```
Module_Score = Σ (prevalence_i + specificity_i) / 200     for observed phenotypes
             - penalty × Σ (prevalence_j + specificity_j) / 400   for excluded phenotypes
```

Where:
- **Prevalence** = % of genes in module annotated with phenotype
- **Specificity** = % of all genes with phenotype that belong to module
- **Penalty** = 0.5 (configurable)

### Gene Ranking Formula

Within a selected module, genes are ranked by:

```
Gene_Score = Σ phenotype_contribution     for observed phenotypes the gene has
           + stability_bonus              if gene is 'core'
           - stability_penalty            if gene is 'unstable'
```

Where:
- `stability_bonus` = 0.1
- `stability_penalty` = 0.05

### Confidence Calculation

Confidence is based on score separation between the top two modules:

```
Confidence = (Score_1 - Score_2) / Score_1
```

High confidence indicates a clear winning module; low confidence suggests multiple plausible diagnoses.

---

## Installation

### Requirements

- Python 3.11+
- pandas
- openpyxl (for Excel file reading)

### Setup

```bash
# Clone or download the project
cd MPV

# Install dependencies
pip install pandas openpyxl

# Verify installation
python verify_framework.py
```

Expected output:
```
VERIFICATION COMPLETE: 14 passed, 0 failed
All tests passed! Framework is ready for use.
```

---

## Usage Guide

### Basic Phenotype Query

```python
from clinical_support import ClinicalSupportEngine

engine = ClinicalSupportEngine()

# Query with observed phenotypes
result = engine.query(
    observed=["Obesity", "Polydactyly", "Rod-cone dystrophy"]
)

print(result.summary())
```

Output:
```
Best Module: 7 (score: 2.940, confidence: 45.91%)
  Genes in module: 26

Top Candidate Genes (20 total):
  - BBS1 (core, support: 3.040)
  - BBS2 (core, support: 3.040)
  - BBS10 (core, support: 3.040)
  ...

Predicted Missing Phenotypes:
  - Impaired fasting glucose (prevalence: 80.8%)
  - Fifth finger distal phalanx clinodactyly (prevalence: 80.8%)
  ...
```

### Query with Exclusions

```python
result = engine.query(
    observed=["Rod-cone dystrophy"],
    excluded=["Hearing impairment", "Obesity", "Polydactyly"]
)
```

### Gene Query

```python
gene_result = engine.query_gene("RPGR")
print(gene_result.summary())
```

Output:
```
Gene: RPGR
Module: 12
Classification: core (stability: 0.857)

Core genes in module: ATF6, PDE6C, PDE6H, CNGA3, CNGB3

Characteristic phenotypes:
  - Monochromacy (100.0%)
  - Undetectable light-adapted electroretinogram (87.5%)
  - Color vision test abnormality (87.5%)
```

### Interactive Decision Tree

```python
from decision_tree import InteractiveSession

session = InteractiveSession()

# Answer questions incrementally
session.answer_yes("Rod-cone dystrophy")
print(session.get_session_summary())

session.answer_yes("Obesity")
print(session.get_session_summary())

# Get next recommended question
next_q = session.get_next_question()
print(f"Next question: {next_q.name}")

# Get current candidate genes
for gene in session.get_candidate_genes()[:5]:
    print(f"  {gene.gene} ({gene.classification})")
```

### Accessing Detailed Results

```python
result = engine.query(observed=["Obesity", "Rod-cone dystrophy"])

# All ranked modules
for match in result.matched_modules[:5]:
    print(f"Module {match.module_id}: score={match.score:.3f}")

# Score explanation
for exp in result.explanation:
    print(f"{exp.phenotype_name}: {exp.contribution:+.3f}")
```

---

## API Reference

### ClinicalSupportEngine

```python
class ClinicalSupportEngine:
    def __init__(self, data_dir=None, scoring_config=None, prediction_config=None)
    
    def query(self, observed=None, excluded=None, top_genes=20, top_predictions=10) -> QueryResult
    
    def query_gene(self, gene: str, top_phenotypes=20) -> GeneQueryResult
    
    def suggest_next_phenotype(self, observed=None, excluded=None) -> PhenotypePrediction
    
    def get_module_summary(self, module_id: int) -> dict
```

### QueryResult

```python
@dataclass
class QueryResult:
    matched_modules: List[ModuleMatch]      # All modules ranked by score
    best_module: ModuleMatch                # Top-scoring module
    candidate_genes: List[GeneCandidate]    # Genes ranked within best module
    predicted_phenotypes: List[PhenotypePrediction]  # Missing phenotype predictions
    explanation: List[ExplainabilityItem]   # Score breakdown
    observed_phenotypes: List[str]          # Resolved HPO IDs
    excluded_phenotypes: List[str]          # Resolved HPO IDs
    unmatched_inputs: List[str]             # Inputs that couldn't be resolved
```

### ModuleMatch

```python
@dataclass
class ModuleMatch:
    module_id: int
    score: float
    confidence: float                       # 0-1 based on score separation
    gene_count: int
    contributing_phenotypes: List[PhenotypeInfo]
    penalized_phenotypes: List[PhenotypeInfo]
```

### GeneCandidate

```python
@dataclass
class GeneCandidate:
    gene: str
    module_id: int
    support_score: float
    stability_score: float
    classification: str                     # 'core', 'peripheral', 'unstable'
    supporting_phenotypes: List[str]
```

---

## Validation Results

### Automated Test Suite

```
======================================================================
VERIFICATION COMPLETE: 14 passed, 0 failed
======================================================================
```

**Tests:**
- Data Loading (4): gene count, module count, phenotype indexing, gene lookup
- Module Scoring (3): BBS match, Usher match, gene ranking
- Clinical API (4): basic query, gene query, empty query, unmatched handling
- Interactive Session (3): creation, answers, incremental scoring

### Classic Case Validation

| Syndrome | Input Phenotypes | Expected | Result | Top Genes |
|----------|------------------|----------|--------|-----------|
| **Bardet-Biedl** | Obesity + Polydactyly + RP | Module 7 | ✅ Module 7 (49.1%) | BBS1-12 |
| **Usher Type 1** | Hearing loss + RP + Vestibular | Module 11 | ✅ Module 11 | MYO7A, USH1C, USH2A |
| **Achromatopsia** | Color blindness + Photophobia + Nystagmus | Module 12 | ✅ Module 12 (32.0%) | CNGA3, CNGB3, PDE6C |
| **CSNB** | Night blindness + Abnormal ERG | Module 10 | ✅ Module 10 | NYX, GRM6, TRPM1 |
| **Isolated RP** | Night blindness + Bone spicules | Module 0 | ✅ Module 0 (23.8%) | PRPF31, RP2, CHM |
| **Choroideremia** | Night blindness + Choroidal atrophy | Module 0 | ✅ Module 0 | CHM (first) |
| **LCA** | Severe vision loss + Absent ERG | Module 4 | ✅ Module 4 | CEP290 |
| **Alström** | Cone-rod + Obesity + Hearing loss + Cardiomyopathy | Module 7 | ✅ Module 7 (37.4%) | ALMS1 (first) |
| **Mitochondrial** | Optic atrophy + Ataxia + Neuropathy | Module 9 | ✅ Module 9 (42.8%) | MT-ND1-6, OPA1 |

---

## Limitations

### Methodological Limitations

1. **Phenotype hierarchy not explicitly modeled**: While hierarchical relationships are implicitly encoded through prevalence weights, explicit parent-child HPO relationships are not traversed at runtime.

2. **Binary phenotype model**: Phenotypes are treated as present/absent; severity gradations are not captured.

3. **Module-level granularity**: The system operates at the module level; within-module genetic heterogeneity is handled only through gene stability scores.

### Data Limitations

1. **Fixed phenotype vocabulary**: Only phenotypes present in the pre-computed profiles can be matched.

2. **Gene set fixed at 442**: Newly discovered IRD genes are not included.

3. **No variant-level information**: The system prioritizes genes, not specific variants.

### Intended Use

This framework is designed as:
- A **demonstrative clinical decision support tool**
- A **proof-of-concept** for modular disease architecture
- A **research tool** for phenotype-driven gene prioritization

**It is NOT intended to replace clinical diagnosis** but to support structured clinical reasoning under incomplete information.

---

## Citation

If you use this framework in your research, please cite:

```
[Citation to be added upon publication]
```

---

## License

This project is licensed under the MIT License.

---

## Contact

For questions or issues, please contact the repository maintainers.

---

## Acknowledgments

This work was developed as part of research into modular disease architecture in Inherited Retinal Diseases at [Institution].

Special thanks to:
- The Human Phenotype Ontology (HPO) project
- The Retinal Information Network (RetNet)
- [Additional acknowledgments]
