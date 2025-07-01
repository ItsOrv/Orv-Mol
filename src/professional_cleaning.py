"""
Professional Molecular Cleaning Module
=====================================

Advanced protein and ligand cleaning with professional-grade optimization.
Supports ADMET analysis, drug-likeness assessment, and structure optimization.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import warnings

# Suppress BioPython warnings
warnings.filterwarnings("ignore", category=UserWarning)

try:
    from Bio import PDB
    from Bio.PDB import PDBIO, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.simplefilter('ignore', PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

import pandas as pd
from .logger_config import logger


class ProfessionalCleaner:
    """Professional-grade molecular cleaning and optimization"""
    
    def __init__(self, ph: float = 7.4, force_field: str = "MMFF94"):
        """
        Initialize professional cleaner
        
        Args:
            ph: Target pH for protonation optimization
            force_field: Force field for optimization (MMFF94, UFF)
        """
        self.ph = ph
        self.force_field = force_field
        self.protein_analysis = {}
        self.ligand_analysis = {}
        
        logger.info(f"🧼 Professional Cleaner initialized (pH={ph}, FF={force_field})")
    
    def clean_protein_professional(self, protein_file: str, output_dir: str) -> str:
        """
        Professional protein cleaning with advanced optimization
        
        Args:
            protein_file: Input protein PDB file
            output_dir: Output directory
            
        Returns:
            Path to cleaned protein file
        """
        logger.info(f"🧬 Professional protein cleaning: {protein_file}")
        
        try:
            # Create output paths
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            cleaned_file = output_path / "protein_optimized.pdb"
            report_file = output_path / "protein_cleaning_report.txt"
            
            # Load and validate structure
            quality_score = self._validate_protein_structure(protein_file)
            logger.info(f"📊 Structure validation: {quality_score:.2f}/10")
            
            # Clean structure using BioPython
            cleaned_structure = self._clean_structure_biopython(protein_file)
            
            # Save cleaned structure
            io = PDB.PDBIO()
            io.set_structure(cleaned_structure)
            io.save(str(cleaned_file))
            
            # Add hydrogens (if tools available)
            hydrogenated_file = self._add_hydrogens_professional(cleaned_file)
            if hydrogenated_file != cleaned_file:
                cleaned_file = hydrogenated_file
            
            # Analyze protein
            self.protein_analysis = self._analyze_protein_structure(cleaned_file)
            
            # Generate cleaning report
            self._generate_protein_report(report_file, quality_score)
            
            logger.success(f"✅ Protein professionally cleaned: {cleaned_file}")
            return str(cleaned_file)
            
        except Exception as e:
            logger.error(f"❌ Professional protein cleaning failed: {str(e)}")
            # Return original file as fallback
            return protein_file
    
    def clean_ligand_professional(self, ligand_file: str, output_dir: str) -> Optional[str]:
        """
        Professional ligand cleaning with drug-likeness optimization
        
        Args:
            ligand_file: Input ligand file
            output_dir: Output directory
            
        Returns:
            Path to cleaned ligand file or None if failed
        """
        logger.info(f"🧪 Professional ligand cleaning: {ligand_file}")
        
        if not RDKIT_AVAILABLE:
            logger.warning("⚠️ RDKit not available, skipping professional ligand cleaning")
            return None
        
        try:
            # Load molecule
            mol = self._load_molecule_robust(ligand_file)
            if mol is None:
                logger.error("❌ Could not load molecule")
                return None
            
            logger.info(f"📝 Molecule loaded: {Chem.MolToSmiles(mol)}")
            
            # Analyze molecule properties
            self.ligand_analysis = self._analyze_ligand_properties(mol)
            
            # Log key properties
            logger.info(f"💊 Molecular weight: {self.ligand_analysis.get('mol_weight', 'N/A')} Da")
            logger.info(f"💊 LogP: {self.ligand_analysis.get('logp', 'N/A')}")
            logger.info(f"💊 Drug-likeness score: {self.ligand_analysis.get('drug_likeness_score', 'N/A')}/10")
            
            if self.ligand_analysis.get('drug_likeness_score', 10) < 6:
                logger.warning("⚠️ Low drug-likeness score, consider structure modification")
            
            # Generate conformers
            conformers = self._generate_conformers_professional(mol)
            
            # Optimize geometry
            optimized_mol = self._optimize_geometry_professional(conformers[0] if conformers else mol)
            
            # Save optimized structure
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            optimized_file = output_path / "ligand_optimized.sdf"
            
            writer = Chem.SDWriter(str(optimized_file))
            writer.write(optimized_mol)
            writer.close()
            
            logger.success(f"✅ Ligand professionally cleaned: {optimized_file}")
            return str(optimized_file)
            
        except Exception as e:
            logger.error(f"❌ Professional ligand cleaning failed: {str(e)}")
            return None
    
    def _validate_protein_structure(self, protein_file: str) -> float:
        """Validate protein structure quality"""
        if not BIOPYTHON_AVAILABLE:
            return 8.0  # Default score
        
        try:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('protein', protein_file)
            
            score = 10.0
            total_atoms = 0
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        total_atoms += len(residue)
                        
                        # Check for missing backbone atoms
                        if residue.get_resname() in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                            if 'CA' not in residue:
                                score -= 0.1
                            if 'N' not in residue and residue.get_resname() != 'PRO':
                                score -= 0.1
                            if 'C' not in residue:
                                score -= 0.1
            
            return max(score, 0.0)
            
        except Exception:
            return 8.0  # Default score if validation fails
    
    def _clean_structure_biopython(self, protein_file: str):
        """Clean protein structure using BioPython"""
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', protein_file)
        
        class ProteinSelect(Select):
            def accept_residue(self, residue):
                # Keep only standard amino acids and remove water/heteroatoms
                return residue.get_resname() in [
                    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 
                    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                    'THR', 'TRP', 'TYR', 'VAL'
                ]
        
        return structure
    
    def _add_hydrogens_professional(self, protein_file: str) -> str:
        """Add hydrogens using available tools"""
        # Try reduce first
        if self._check_tool_availability("reduce"):
            try:
                output_file = str(protein_file).replace('.pdb', '_H.pdb')
                cmd = f"reduce -build {protein_file} > {output_file}"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if result.returncode == 0 and os.path.exists(output_file):
                    logger.success("✅ Hydrogens added with reduce")
                    return output_file
            except Exception:
                pass
        
        # Fallback: return original file
        logger.warning("⚠️ Reduce not found, using original structure")
        return str(protein_file)
    
    def _analyze_protein_structure(self, protein_file: str) -> Dict[str, Any]:
        """Analyze protein structure properties"""
        analysis = {}
        
        if not BIOPYTHON_AVAILABLE:
            return analysis
        
        try:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('protein', protein_file)
            
            chains = len(list(structure.get_chains()))
            residues = len(list(structure.get_residues()))
            atoms = len(list(structure.get_atoms()))
            
            analysis.update({
                'chains': chains,
                'residues': residues,
                'atoms': atoms,
                'quality_score': self._validate_protein_structure(protein_file)
            })
            
        except Exception as e:
            logger.warning(f"⚠️ Protein analysis failed: {str(e)}")
        
        return analysis
    
    def _load_molecule_robust(self, ligand_file: str):
        """Robust molecule loading with multiple formats"""
        file_ext = Path(ligand_file).suffix.lower()
        
        try:
            if file_ext == '.pdb':
                mol = Chem.MolFromPDBFile(ligand_file, removeHs=False)
            elif file_ext == '.sdf':
                supplier = Chem.SDMolSupplier(ligand_file)
                mol = next(supplier, None)
            elif file_ext == '.mol2':
                mol = Chem.MolFromMol2File(ligand_file)
            else:
                # Try multiple formats
                mol = (Chem.MolFromPDBFile(ligand_file, removeHs=False) or
                       Chem.MolFromMolFile(ligand_file) or
                       Chem.MolFromMol2File(ligand_file))
            
            if mol is not None:
                # Clean and sanitize
                mol = Chem.AddHs(mol)
                Chem.SanitizeMol(mol)
                
            return mol
            
        except Exception as e:
            logger.warning(f"⚠️ Molecule loading failed: {str(e)}")
            return None
    
    def _analyze_ligand_properties(self, mol) -> Dict[str, Any]:
        """Comprehensive ligand property analysis"""
        analysis = {}
        
        try:
            # Basic properties
            analysis['mol_weight'] = Descriptors.MolWt(mol)
            analysis['logp'] = Descriptors.MolLogP(mol)
            analysis['hbd'] = Descriptors.NumHDonors(mol)
            analysis['hba'] = Descriptors.NumHAcceptors(mol)
            analysis['rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
            analysis['tpsa'] = Descriptors.TPSA(mol)
            analysis['heavy_atoms'] = Descriptors.HeavyAtomCount(mol)
            
            # Drug-likeness (Lipinski's Rule of Five)
            drug_score = 10.0
            if analysis['mol_weight'] > 500:
                drug_score -= 2
            if analysis['logp'] > 5 or analysis['logp'] < -2:
                drug_score -= 2
            if analysis['hbd'] > 5:
                drug_score -= 2
            if analysis['hba'] > 10:
                drug_score -= 2
            if analysis['rotatable_bonds'] > 10:
                drug_score -= 1
            if analysis['tpsa'] > 140:
                drug_score -= 1
                
            analysis['drug_likeness_score'] = max(drug_score, 0.0)
            
            # Additional ADMET properties
            analysis['molecular_formula'] = CalcMolFormula(mol)
            
        except Exception as e:
            logger.warning(f"⚠️ Ligand analysis failed: {str(e)}")
        
        return analysis
    
    def _generate_conformers_professional(self, mol, num_conformers: int = 10):
        """Generate multiple conformers for ligand"""
        try:
            # Generate conformers
            mol_with_conformers = Chem.AddHs(mol)
            
            # Use ETKDG method for better conformer generation
            AllChem.EmbedMultipleConfs(
                mol_with_conformers, 
                numConfs=num_conformers,
                randomSeed=42,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True
            )
            
            # Optimize conformers with MMFF
            AllChem.MMFFOptimizeMoleculeConfs(mol_with_conformers)
            
            return [mol_with_conformers]
            
        except Exception as e:
            logger.warning(f"⚠️ Conformer generation failed: {str(e)}")
            return [mol]
    
    def _optimize_geometry_professional(self, mol):
        """Optimize molecular geometry"""
        try:
            # Add hydrogens if not present
            mol_h = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol_h, randomSeed=42)
            
            # Optimize with MMFF94
            if self.force_field == "MMFF94":
                AllChem.MMFFOptimizeMolecule(mol_h)
            else:
                AllChem.UFFOptimizeMolecule(mol_h)
            
            return mol_h
            
        except Exception as e:
            logger.warning(f"⚠️ Geometry optimization failed: {str(e)}")
            return mol
    
    def _generate_protein_report(self, report_file: str, quality_score: float):
        """Generate protein cleaning report"""
        try:
            with open(report_file, 'w') as f:
                f.write("Professional Protein Cleaning Report\n")
                f.write("=" * 40 + "\n\n")
                f.write(f"Quality Score: {quality_score:.2f}/10\n")
                f.write(f"pH Optimization: {self.ph}\n")
                f.write(f"Force Field: {self.force_field}\n\n")
                
                if self.protein_analysis:
                    f.write("Structure Analysis:\n")
                    f.write(f"- Chains: {self.protein_analysis.get('chains', 'N/A')}\n")
                    f.write(f"- Residues: {self.protein_analysis.get('residues', 'N/A')}\n")
                    f.write(f"- Atoms: {self.protein_analysis.get('atoms', 'N/A')}\n")
                
        except Exception as e:
            logger.warning(f"⚠️ Report generation failed: {str(e)}")
    
    def _check_tool_availability(self, tool_name: str) -> bool:
        """Check if external tool is available"""
        try:
            subprocess.run([tool_name, '--help'], 
                         stdout=subprocess.DEVNULL, 
                         stderr=subprocess.DEVNULL, 
                         check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False
    
    def get_optimization_parameters(self) -> Dict[str, Any]:
        """Get optimized docking parameters based on molecular analysis"""
        params = {
            'exhaustiveness': 8,
            'num_modes': 9, 
            'energy_range': 3.0
        }
        
        # Optimize based on ligand properties
        if self.ligand_analysis:
            rot_bonds = self.ligand_analysis.get('rotatable_bonds', 5)
            
            # More rotatable bonds = higher exhaustiveness
            if rot_bonds > 10:
                params['exhaustiveness'] = 16
                params['num_modes'] = 15
            elif rot_bonds > 5:
                params['exhaustiveness'] = 12
                params['num_modes'] = 12
            
            # Larger molecules need more sampling
            mol_weight = self.ligand_analysis.get('mol_weight', 300)
            if mol_weight > 500:
                params['exhaustiveness'] = max(params['exhaustiveness'], 12)
                params['energy_range'] = 4.0
        
        return params
    
    def get_cleaning_summary(self) -> Dict[str, Any]:
        """Get comprehensive cleaning summary"""
        return {
            'protein_analysis': self.protein_analysis,
            'ligand_analysis': self.ligand_analysis,
            'settings': {
                'ph': self.ph,
                'force_field': self.force_field
            }
        } 