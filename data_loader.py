"""
Data Loader Module for the Clinical Decision Support Framework.

Loads and indexes:
- Module phenotype profiles from Excel file
- Gene-to-module mappings from CSV file
- Builds reverse indexes for efficient lookup
"""

import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Set, Optional, Tuple
import re


@dataclass
class PhenotypeProfile:
    """Profile of a single phenotype within a module."""
    name: str
    hpo_id: str
    genes_with: List[str]  # Genes in module with this phenotype
    genes_without: List[str]  # Genes in module without this phenotype
    prevalence: float  # target_module_phenotype_prevalence_percent
    specificity: float  # target_module_share_of_phenotype_percent
    non_target_ird_genes: List[str]  # IRD genes outside module with this phenotype
    non_ird_genes: List[str]  # Non-IRD genes with this phenotype


@dataclass
class ModuleProfile:
    """Complete profile for a disease module."""
    module_id: int
    phenotypes: Dict[str, PhenotypeProfile]  # keyed by HPO ID
    all_genes: Set[str] = field(default_factory=set)
    
    def get_phenotype_by_name(self, name: str) -> Optional[PhenotypeProfile]:
        """Find phenotype by name (case-insensitive)."""
        name_lower = name.lower()
        for pheno in self.phenotypes.values():
            if pheno.name.lower() == name_lower:
                return pheno
        return None
    
    def get_phenotype(self, identifier: str) -> Optional[PhenotypeProfile]:
        """Find phenotype by HPO ID or name."""
        # Try HPO ID first
        if identifier in self.phenotypes:
            return self.phenotypes[identifier]
        # Try name lookup
        return self.get_phenotype_by_name(identifier)


@dataclass
class GeneInfo:
    """Information about a single gene."""
    gene: str
    module_id: int
    stability_score: float
    classification: str  # core, peripheral, unstable


class DataLoader:
    """Loads and manages clinical decision support data."""
    
    def __init__(self, data_dir: Optional[Path] = None):
        """Initialize with optional data directory path."""
        if data_dir is None:
            data_dir = Path(__file__).parent
        self.data_dir = Path(data_dir)
        
        # Data storage
        self.module_profiles: Dict[int, ModuleProfile] = {}
        self.gene_mapping: Dict[str, GeneInfo] = {}
        
        # Indexes
        self.phenotype_to_modules: Dict[str, Set[int]] = {}  # HPO ID -> modules
        self.phenotype_name_to_hpo: Dict[str, str] = {}  # name (lower) -> HPO ID
        self.module_genes: Dict[int, Set[str]] = {}  # module_id -> genes
        
        self._loaded = False
    
    def load(self, 
             excel_file: str = "module_top500_HPO_background_comparison_20260128_1907.xlsx",
             csv_file: str = "gene_classification_20251230_1900.csv") -> None:
        """Load all data files and build indexes."""
        self._load_gene_mapping(csv_file)
        self._load_module_profiles(excel_file)
        self._build_indexes()
        self._loaded = True
    
    def _parse_gene_list(self, value) -> List[str]:
        """Parse a gene list from Excel cell (comma or space separated)."""
        if pd.isna(value) or value == '':
            return []
        if isinstance(value, str):
            # Split by comma or whitespace
            genes = re.split(r'[,\s]+', value.strip())
            return [g.strip() for g in genes if g.strip()]
        return []
    
    def _load_gene_mapping(self, csv_file: str) -> None:
        """Load gene-to-module mapping from CSV."""
        filepath = self.data_dir / csv_file
        df = pd.read_csv(filepath)
        
        for _, row in df.iterrows():
            gene = row['gene']
            if pd.isna(gene) or gene == '':
                continue
            
            self.gene_mapping[gene] = GeneInfo(
                gene=gene,
                module_id=int(row['module_id']),
                stability_score=float(row['stability_score']),
                classification=row['classification']
            )
            
            # Track module genes
            module_id = int(row['module_id'])
            if module_id not in self.module_genes:
                self.module_genes[module_id] = set()
            self.module_genes[module_id].add(gene)
        
        print(f"Loaded {len(self.gene_mapping)} genes across {len(self.module_genes)} modules")
    
    def _load_module_profiles(self, excel_file: str) -> None:
        """Load module phenotype profiles from Excel file."""
        filepath = self.data_dir / excel_file
        xl = pd.ExcelFile(filepath)
        
        for sheet_name in xl.sheet_names:
            # Extract module ID from sheet name (e.g., "module_0" -> 0)
            match = re.match(r'module_(\d+)', sheet_name)
            if not match:
                continue
            
            module_id = int(match.group(1))
            df = pd.read_excel(xl, sheet_name=sheet_name)
            
            phenotypes = {}
            all_genes = set()
            
            for _, row in df.iterrows():
                hpo_id = row.get('hpo_id', '')
                if pd.isna(hpo_id) or hpo_id == '':
                    continue
                
                genes_with = self._parse_gene_list(row.get('target_module_genes_with_phenotype', ''))
                genes_without = self._parse_gene_list(row.get('target_module_genes_without_phenotype', ''))
                
                all_genes.update(genes_with)
                all_genes.update(genes_without)
                
                phenotypes[hpo_id] = PhenotypeProfile(
                    name=str(row.get('phenotype_name', '')),
                    hpo_id=hpo_id,
                    genes_with=genes_with,
                    genes_without=genes_without,
                    prevalence=float(row.get('target_module_phenotype_prevalence_percent', 0)),
                    specificity=float(row.get('target_module_share_of_phenotype_percent', 0)),
                    non_target_ird_genes=self._parse_gene_list(
                        row.get('non_target_ird_genes_with_phenotype', '')),
                    non_ird_genes=self._parse_gene_list(
                        row.get('non_ird_genes_with_phenotype', ''))
                )
            
            self.module_profiles[module_id] = ModuleProfile(
                module_id=module_id,
                phenotypes=phenotypes,
                all_genes=all_genes
            )
        
        print(f"Loaded {len(self.module_profiles)} module profiles")
    
    def _build_indexes(self) -> None:
        """Build reverse indexes for efficient lookup."""
        for module_id, profile in self.module_profiles.items():
            for hpo_id, pheno in profile.phenotypes.items():
                # HPO ID -> modules index
                if hpo_id not in self.phenotype_to_modules:
                    self.phenotype_to_modules[hpo_id] = set()
                self.phenotype_to_modules[hpo_id].add(module_id)
                
                # Name -> HPO ID index (case-insensitive)
                name_lower = pheno.name.lower()
                if name_lower not in self.phenotype_name_to_hpo:
                    self.phenotype_name_to_hpo[name_lower] = hpo_id
        
        print(f"Indexed {len(self.phenotype_to_modules)} unique phenotypes")
    
    def resolve_phenotype(self, identifier: str) -> Optional[str]:
        """
        Resolve a phenotype identifier to HPO ID.
        Accepts HPO ID (HP:XXXXXXX) or phenotype name.
        Returns None if not found.
        """
        # Check if it's already an HPO ID
        if identifier.startswith('HP:'):
            if identifier in self.phenotype_to_modules:
                return identifier
            return None
        
        # Try name lookup
        name_lower = identifier.lower()
        return self.phenotype_name_to_hpo.get(name_lower)
    
    def get_gene_info(self, gene: str) -> Optional[GeneInfo]:
        """Get information about a gene."""
        return self.gene_mapping.get(gene)
    
    def get_module_profile(self, module_id: int) -> Optional[ModuleProfile]:
        """Get the phenotype profile for a module."""
        return self.module_profiles.get(module_id)
    
    def get_module_genes(self, module_id: int) -> Set[str]:
        """Get all genes in a module."""
        return self.module_genes.get(module_id, set())
    
    def get_all_module_ids(self) -> List[int]:
        """Get list of all module IDs."""
        return sorted(self.module_profiles.keys())
    
    def get_phenotype_info(self, hpo_id: str, module_id: int) -> Optional[PhenotypeProfile]:
        """Get phenotype information within a specific module."""
        profile = self.module_profiles.get(module_id)
        if profile:
            return profile.phenotypes.get(hpo_id)
        return None
    
    def get_all_phenotypes(self) -> Dict[str, str]:
        """
        Get all known phenotypes for autocomplete.
        Returns dict: {"Name (HPO:ID)": "HPO:ID"} sorted by name.
        """
        phenotypes = {}
        for profile in self.module_profiles.values():
            for pheno in profile.phenotypes.values():
                display_name = f"{pheno.name} ({pheno.hpo_id})"
                phenotypes[display_name] = pheno.hpo_id
        
        # Sort by display name
        return dict(sorted(phenotypes.items()))


# Singleton instance for convenience
_default_loader: Optional[DataLoader] = None


def get_loader(data_dir: Optional[Path] = None) -> DataLoader:
    """Get or create the default data loader."""
    global _default_loader
    if _default_loader is None:
        _default_loader = DataLoader(data_dir)
        _default_loader.load()
    return _default_loader


if __name__ == "__main__":
    # Test loading
    loader = DataLoader()
    loader.load()
    
    print(f"\nModule IDs: {loader.get_all_module_ids()}")
    print(f"\nSample gene lookup (RPGR): {loader.get_gene_info('RPGR')}")
    
    # Test phenotype resolution
    test_phenotypes = ["Retinitis pigmentosa", "HP:0000510", "Visual impairment"]
    for pheno in test_phenotypes:
        resolved = loader.resolve_phenotype(pheno)
        print(f"'{pheno}' -> {resolved}")
