"""
Classic IRD Cases - Textbook Presentations
Cases that every ophthalmologist/geneticist would recognize immediately.
"""
from clinical_support import ClinicalSupportEngine

engine = ClinicalSupportEngine()

print("="*70)
print("      CLASSIC IRD CASES - TEXTBOOK PRESENTATIONS")
print("="*70)

# Case 1: Classic Usher Type 1
print("\n" + "="*70)
print("CASE 1: USHER SYNDROME TYPE 1 (Classic)")
print("Congenital profound deafness + RP onset in first decade + Vestibular areflexia")
print("="*70)
result = engine.query(observed=[
    "Profound sensorineural hearing impairment",
    "Rod-cone dystrophy",
    "Vestibular areflexia",
    "Delayed motor development"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes:")
for g in result.candidate_genes[:8]:
    print(f"  {g.gene} ({g.classification})")

# Case 2: Classic Leber Congenital Amaurosis
print("\n" + "="*70)
print("CASE 2: LEBER CONGENITAL AMAUROSIS (LCA)")
print("Severe vision loss from birth + Nystagmus + Absent ERG + Eye poking")
print("="*70)
result = engine.query(observed=[
    "Severe visual impairment",
    "Nystagmus",
    "Absent electroretinogram",
    "Oculodigital reflex"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes (LCA genes expected):")
for g in result.candidate_genes[:10]:
    print(f"  {g.gene} ({g.classification})")

# Case 3: Classic Stargardt Disease
print("\n" + "="*70)
print("CASE 3: STARGARDT DISEASE")
print("Central vision loss in teens + Yellow flecks + Dark choroid on FA")
print("="*70)
result = engine.query(observed=[
    "Macular degeneration",
    "Central scotoma",
    "Photophobia"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes:")
for g in result.candidate_genes[:5]:
    print(f"  {g.gene} ({g.classification})")
# Check for ABCA4
gene_result = engine.query_gene("ABCA4")
if gene_result:
    print(f"\n>>> ABCA4 is in Module {gene_result.module_id} ({gene_result.classification})")

# Case 4: Classic Bardet-Biedl Syndrome
print("\n" + "="*70)
print("CASE 4: BARDET-BIEDL SYNDROME (BBS)")  
print("RP + Obesity + Polydactyly + Hypogonadism + Renal anomaly + Cognitive issues")
print("="*70)
result = engine.query(observed=[
    "Rod-cone dystrophy",
    "Obesity", 
    "Postaxial polydactyly",
    "Hypogonadism",
    "Renal abnormality",
    "Intellectual disability"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes (BBS genes expected):")
for g in result.candidate_genes[:12]:
    print(f"  {g.gene} ({g.classification})")

# Case 5: Classic Retinitis Pigmentosa (isolated)
print("\n" + "="*70)
print("CASE 5: ISOLATED RETINITIS PIGMENTOSA")
print("Night blindness + Peripheral vision loss + Bone spicules + Attenuated vessels")
print("="*70)
result = engine.query(observed=[
    "Rod-cone dystrophy",
    "Nyctalopia",
    "Peripheral visual field loss",
    "Bone spicule pigmentation of the retina"
], excluded=[
    "Hearing impairment",
    "Obesity",
    "Polydactyly"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes:")
for g in result.candidate_genes[:10]:
    print(f"  {g.gene} ({g.classification})")

# Case 6: Classic Achromatopsia
print("\n" + "="*70)
print("CASE 6: ACHROMATOPSIA (Rod Monochromacy)")
print("Complete color blindness + Photophobia + Nystagmus + Reduced VA")
print("="*70)
result = engine.query(observed=[
    "Achromatopsia",
    "Photophobia",
    "Nystagmus",
    "Reduced visual acuity"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes (CNGA3, CNGB3, PDE6C expected):")
for g in result.candidate_genes[:6]:
    print(f"  {g.gene} ({g.classification})")

# Case 7: Classic Choroideremia
print("\n" + "="*70)
print("CASE 7: CHOROIDEREMIA")
print("X-linked + Night blindness + Progressive peripheral loss + Choroidal atrophy")
print("="*70)
result = engine.query(observed=[
    "Nyctalopia",
    "Chorioretinal atrophy",
    "Peripheral visual field loss"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
gene_result = engine.query_gene("CHM")
if gene_result:
    print(f"\n>>> CHM gene is in Module {gene_result.module_id} ({gene_result.classification})")
print("\nTop candidate genes:")
for g in result.candidate_genes[:8]:
    print(f"  {g.gene} ({g.classification})")

# Case 8: LHON (Leber Hereditary Optic Neuropathy)
print("\n" + "="*70)
print("CASE 8: LHON (Leber Hereditary Optic Neuropathy)")
print("Acute bilateral vision loss + Central scotoma + Mitochondrial inheritance")
print("="*70)
result = engine.query(observed=[
    "Optic atrophy",
    "Central scotoma",
    "Reduced visual acuity"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes (MT-ND genes expected for LHON):")
for g in result.candidate_genes[:10]:
    print(f"  {g.gene} ({g.classification})")

# Case 9: Alstrom Syndrome
print("\n" + "="*70)
print("CASE 9: ALSTROM SYNDROME")
print("Cone-rod dystrophy + Obesity + Hearing loss + Cardiomyopathy + Diabetes")
print("="*70)
result = engine.query(observed=[
    "Cone/cone-rod dystrophy",
    "Obesity",
    "Sensorineural hearing impairment",
    "Cardiomyopathy",
    "Diabetes mellitus"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
gene_result = engine.query_gene("ALMS1")
if gene_result:
    print(f"\n>>> ALMS1 gene is in Module {gene_result.module_id} ({gene_result.classification})")
print("\nTop candidate genes:")
for g in result.candidate_genes[:5]:
    print(f"  {g.gene} ({g.classification})")

# Case 10: Congenital Stationary Night Blindness
print("\n" + "="*70)
print("CASE 10: CONGENITAL STATIONARY NIGHT BLINDNESS (CSNB)")
print("Non-progressive night blindness + Abnormal ERG + Normal fundus")
print("="*70)
result = engine.query(observed=[
    "Nyctalopia",
    "Abnormal electroretinogram",
    "Myopia"
])
print(f"\nBest Module: {result.best_module.module_id} (confidence: {result.best_module.confidence:.1%})")
print("\nTop candidate genes (CSNB genes: NYX, GRM6, TRPM1, etc.):")
for g in result.candidate_genes[:10]:
    print(f"  {g.gene} ({g.classification})")

print("\n" + "="*70)
print("END OF CLASSIC CASES")
print("="*70)
