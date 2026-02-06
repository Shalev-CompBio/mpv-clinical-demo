"""
Challenging test cases for the Clinical Decision Support Framework.
"""
from clinical_support import ClinicalSupportEngine

engine = ClinicalSupportEngine()

print("="*70)
print("CHALLENGE 1: Single vague phenotype")
print("Input: Only 'Visual impairment'")
print("="*70)
result = engine.query(observed=["Visual impairment"])
print(f"Best Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("Top 3 modules:")
for m in result.matched_modules[:3]:
    print(f"  Module {m.module_id}: score={m.score:.3f}")
print()

print("="*70)
print("CHALLENGE 2: Conflicting phenotypes (BBS + Usher features)")
print("Input: Obesity + Hearing impairment + Rod-cone dystrophy")
print("="*70)
result = engine.query(observed=["Obesity", "Sensorineural hearing impairment", "Rod-cone dystrophy"])
print(f"Best Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("Top 3 modules:")
for m in result.matched_modules[:3]:
    print(f"  Module {m.module_id}: score={m.score:.3f}")
print("Top genes:")
for g in result.candidate_genes[:5]:
    print(f"  {g.gene} ({g.classification})")
print()

print("="*70)
print("CHALLENGE 3: Rare phenotype combination")
print("Input: Optic atrophy + Ataxia + Peripheral neuropathy")
print("="*70)
result = engine.query(observed=["Optic atrophy", "Ataxia", "Peripheral neuropathy"])
print(f"Best Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("Top genes:")
for g in result.candidate_genes[:5]:
    print(f"  {g.gene} ({g.classification})")
print()

print("="*70)
print("CHALLENGE 4: Exclusion-heavy query")
print("Input: Rod-cone dystrophy, EXCLUDE: Hearing, Obesity, Polydactyly, Ataxia")
print("="*70)
result = engine.query(
    observed=["Rod-cone dystrophy"],
    excluded=["Hearing impairment", "Obesity", "Polydactyly", "Ataxia"]
)
print(f"Best Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print(f"Score after exclusions: {result.best_module.score:.3f}")
print("Top 3 modules:")
for m in result.matched_modules[:3]:
    print(f"  Module {m.module_id}: score={m.score:.3f}")
print()

print("="*70)
print("CHALLENGE 5: Mitochondrial-like features")  
print("Input: Optic atrophy + Hearing impairment + Muscle weakness")
print("="*70)
result = engine.query(observed=["Optic atrophy", "Sensorineural hearing impairment", "Muscle weakness"])
print(f"Best Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("Top genes:")
for g in result.candidate_genes[:8]:
    print(f"  {g.gene} ({g.classification})")
print()

print("="*70)
print("CHALLENGE 6: Very specific - Achromatopsia features")
print("Input: Photophobia + Color blindness + Nystagmus")
print("="*70)
result = engine.query(observed=["Photophobia", "Color blindness", "Nystagmus"])
print(f"Best Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("Top genes:")
for g in result.candidate_genes[:5]:
    print(f"  {g.gene} ({g.classification})")
print()

print("="*70)
print("CHALLENGE 7: Ciliopathy overlap - Joubert vs BBS")
print("Input: Intellectual disability + Polydactyly + Renal anomaly")
print("="*70)
result = engine.query(observed=["Intellectual disability", "Polydactyly", "Renal abnormality"])
print(f"Best Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("Top 5 modules (ciliopathy battle):")
for m in result.matched_modules[:5]:
    print(f"  Module {m.module_id}: score={m.score:.3f}")
