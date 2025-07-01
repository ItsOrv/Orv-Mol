"""
Enhanced Preprocessing module for molecular docking pipeline.
Handles ligand and protein preparation, format conversion, binding box calculation,
and professional-grade cleaning and optimization.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple, Optional

import numpy as np
from .logger_config import logger
from Bio.PDB import PDBParser, PDBIO
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Import professional cleaning capabilities
from .professional_cleaning import ProfessionalCleaner


class PreprocessingEngine:
    """
    Enhanced preprocessing engine with professional-grade cleaning capabilities.
    
    Features:
    - Standard preprocessing (backward compatibility)
    - Professional-grade protein and ligand cleaning
    - Advanced molecular property analysis
    - Optimized docking parameter suggestions
    """
    
    def __init__(self, professional_mode: bool = True, ph: float = 7.4):
        """
        Initialize preprocessing engine.
        
        Args:
            professional_mode: Enable professional cleaning features
            ph: Target pH for protonation optimization
        """
        self.parser = PDBParser(QUIET=True)
        self.professional_mode = professional_mode
        
        if professional_mode:
            self.professional_cleaner = ProfessionalCleaner(ph=ph)
            logger.info("🧼 Professional cleaning mode enabled")
        else:
            logger.info("⚙️ Standard preprocessing mode")
        
    def prepare_ligand(
        self, 
        ligand_file: str, 
        output_dir: str,
        professional_cleaning: bool = None
    ) -> str:
        """
        Enhanced ligand preparation with optional professional cleaning.
        
        Args:
            ligand_file: Path to input ligand file (.mol2, .sdf, .pdb)
            output_dir: Directory to save processed files
            professional_cleaning: Override professional mode for this ligand
            
        Returns:
            Path to prepared PDBQT file
        """
        ligand_path = Path(ligand_file)
        output_path = Path(output_dir)
        
        # Generate output filename
        ligand_pdbqt = output_path / f"{ligand_path.stem}.pdbqt"
        
        logger.info(f"🧪 Preparing ligand: {ligand_path.name}")
        
        # Determine if we should use professional cleaning
        use_professional = (
            professional_cleaning 
            if professional_cleaning is not None 
            else self.professional_mode
        )
        
        try:
            if use_professional and hasattr(self, 'professional_cleaner'):
                # Professional ligand cleaning
                logger.info("🧼 Using professional ligand cleaning...")
                
                cleaned_ligand_file = self.professional_cleaner.clean_ligand_professional(
                    ligand_file=ligand_file,
                    output_dir=output_dir
                )
                
                if cleaned_ligand_file:
                    # Convert cleaned file to PDBQT
                    if cleaned_ligand_file.endswith('.sdf'):
                        self._convert_to_pdbqt_ligand(cleaned_ligand_file, ligand_pdbqt)
                    else:
                        # Convert via MOL2 first
                        mol2_file = self._convert_to_mol2(cleaned_ligand_file, output_dir)
                        self._convert_to_pdbqt_ligand(mol2_file, ligand_pdbqt)
                    
                    # Store cleaning results for later use
                    self._last_ligand_properties = self.professional_cleaner.ligand_analysis
                    
                    drug_score = self._last_ligand_properties.get('drug_likeness_score', 'N/A')
                    mol_weight = self._last_ligand_properties.get('mol_weight', 'N/A')
                    logger.info(f"💊 Drug-likeness: {drug_score}/10")
                    logger.info(f"💊 Molecular weight: {mol_weight} Da")
                else:
                    # Fallback to standard method
                    logger.warning("⚠️ Professional cleaning failed, using standard method...")
                    mol2_file = self._convert_to_mol2(ligand_file, output_dir)
                    minimized_mol2 = self._minimize_ligand(mol2_file, output_dir)
                    self._convert_to_pdbqt_ligand(minimized_mol2, ligand_pdbqt)
                
            else:
                # Standard preprocessing (backward compatibility)
                logger.info("⚙️ Using standard ligand preparation...")
                
                # Step 1: Convert to MOL2 if needed and minimize
                mol2_file = self._convert_to_mol2(ligand_file, output_dir)
                
                # Step 2: Minimize ligand structure
                minimized_mol2 = self._minimize_ligand(mol2_file, output_dir)
                
                # Step 3: Convert to PDBQT
                self._convert_to_pdbqt_ligand(minimized_mol2, ligand_pdbqt)
            
            logger.success(f"✅ Ligand prepared: {ligand_pdbqt}")
            return str(ligand_pdbqt)
            
        except Exception as e:
            logger.error(f"❌ Failed to prepare ligand: {e}")
            raise
    
    def prepare_receptor(
        self, 
        protein_file: str, 
        output_dir: str,
        professional_cleaning: bool = None
    ) -> str:
        """
        Enhanced protein receptor preparation with optional professional cleaning.
        
        Args:
            protein_file: Path to input protein PDB file
            output_dir: Directory to save processed files
            professional_cleaning: Override professional mode for this protein
            
        Returns:
            Path to prepared PDBQT file
        """
        protein_path = Path(protein_file)
        output_path = Path(output_dir)
        
        # Generate output filenames
        clean_pdb = output_path / f"{protein_path.stem}_clean.pdb"
        protein_pdbqt = output_path / f"{protein_path.stem}.pdbqt"
        
        logger.info(f"🧬 Preparing protein receptor: {protein_path.name}")
        
        # Determine if we should use professional cleaning
        use_professional = (
            professional_cleaning 
            if professional_cleaning is not None 
            else self.professional_mode
        )
        
        try:
            if use_professional and hasattr(self, 'professional_cleaner'):
                # Professional protein cleaning
                logger.info("🧼 Using professional protein cleaning...")
                
                cleaned_file = self.professional_cleaner.clean_protein_professional(
                    protein_file=protein_file,
                    output_dir=output_dir
                )
                
                # Convert cleaned PDB to PDBQT
                self._convert_to_pdbqt_receptor(cleaned_file, protein_pdbqt)
                
                # Store cleaning results for later use
                self._last_protein_analysis = self.professional_cleaner.protein_analysis
                
                quality_score = self._last_protein_analysis.get('quality_score', 0)
                chains = self._last_protein_analysis.get('chains', 'N/A')
                residues = self._last_protein_analysis.get('residues', 'N/A')
                
                logger.info(f"📊 Quality score: {quality_score}/10")
                logger.info(f"🧬 Chains: {chains}")
                logger.info(f"🧬 Residues: {residues}")
                
            else:
                # Standard preprocessing (backward compatibility)
                logger.info("⚙️ Using standard protein preparation...")
                
                # Step 1: Clean protein (remove waters, add hydrogens)
                self._clean_protein(protein_file, clean_pdb)
                
                # Step 2: Convert to PDBQT
                self._convert_to_pdbqt_receptor(clean_pdb, protein_pdbqt)
            
            logger.success(f"✅ Receptor prepared: {protein_pdbqt}")
            return str(protein_pdbqt)
            
        except Exception as e:
            logger.error(f"❌ Failed to prepare receptor: {e}")
            raise
    
    def calculate_blind_docking_box(self, protein_pdbqt: str) -> Dict[str, float]:
        """
        Calculate binding box parameters for blind docking.
        
        Args:
            protein_pdbqt: Path to prepared protein PDBQT file
            
        Returns:
            Dictionary with box center and size parameters
        """
        logger.info("📏 Calculating blind docking box...")
        
        try:
            # Read coordinates from PDBQT file
            coords = self._extract_coordinates_from_pdbqt(protein_pdbqt)
            
            if len(coords) == 0:
                raise ValueError("No coordinates found in protein file")
            
            coords = np.array(coords)
            
            # Calculate bounding box
            min_coords = np.min(coords, axis=0)
            max_coords = np.max(coords, axis=0)
            
            # Calculate center
            center = (min_coords + max_coords) / 2
            
            # Calculate size with some padding
            size = max_coords - min_coords
            padding = 5.0  # Angstroms padding
            size += padding
            
            # Ensure minimum size for blind docking
            min_size = 20.0
            size = np.maximum(size, min_size)
            
            box_params = {
                'center_x': float(center[0]),
                'center_y': float(center[1]),
                'center_z': float(center[2]),
                'size_x': float(size[0]),
                'size_y': float(size[1]),
                'size_z': float(size[2])
            }
            
            logger.info(f"📦 Box calculated - Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
            logger.info(f"📦 Box size: ({size[0]:.2f}, {size[1]:.2f}, {size[2]:.2f})")
            
            return box_params
            
        except Exception as e:
            logger.error(f"❌ Failed to calculate binding box: {e}")
            raise
    
    def get_optimized_docking_parameters(self) -> Dict[str, float]:
        """
        Get optimized docking parameters based on professional analysis.
        
        Returns:
            Dictionary with optimized docking parameters
        """
        if not self.professional_mode or not hasattr(self, 'professional_cleaner'):
            logger.warning("Professional mode not enabled, returning default parameters")
            return {
                'exhaustiveness': 8,
                'num_modes': 9,
                'energy_range': 3.0
            }
        
        # Get stored analysis results
        protein_analysis = getattr(self, '_last_protein_analysis', {})
        ligand_properties = getattr(self, '_last_ligand_properties', {})
        
        if not protein_analysis and not ligand_properties:
            logger.warning("No analysis data available, returning default parameters")
            return {
                'exhaustiveness': 8,
                'num_modes': 9,
                'energy_range': 3.0
            }
        
        # Use professional cleaner to optimize parameters
        optimized_params = self.professional_cleaner.get_optimization_parameters()
        
        return optimized_params
    
    def get_cleaning_summary(self) -> Dict[str, any]:
        """
        Get summary of professional cleaning results.
        
        Returns:
            Dictionary with cleaning summary
        """
        summary = {
            'professional_mode': self.professional_mode,
            'protein_analysis': getattr(self, '_last_protein_analysis', None),
            'ligand_properties': getattr(self, '_last_ligand_properties', None)
        }
        
        return summary
    
    def _convert_to_mol2(self, input_file: str, output_dir: str) -> str:
        """Convert ligand to MOL2 format using Open Babel."""
        input_path = Path(input_file)
        output_path = Path(output_dir) / f"{input_path.stem}.mol2"
        
        if input_path.suffix.lower() == '.mol2':
            # Already MOL2, just copy
            import shutil
            shutil.copy2(input_file, output_path)
            return str(output_path)
        
        # Use RDKit for conversion
        try:
            if input_path.suffix.lower() == '.sdf':
                mol = Chem.SDMolSupplier(str(input_path))[0]
            elif input_path.suffix.lower() == '.pdb':
                mol = Chem.MolFromPDBFile(str(input_path))
            else:
                raise ValueError(f"Unsupported input format: {input_path.suffix}")
            
            if mol is None:
                raise ValueError("Failed to read molecule")
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Write to MOL2 using Open Babel if available
            cmd = [
                'obabel',
                str(input_file),
                '-O', str(output_path),
                '--gen3d'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.warning("Open Babel failed, using RDKit fallback")
                # Fallback: use RDKit to generate SDF then convert
                sdf_file = output_path.with_suffix('.sdf')
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)
                writer = Chem.SDWriter(str(sdf_file))
                writer.write(mol)
                writer.close()
                
                # Convert SDF to MOL2
                cmd = ['obabel', str(sdf_file), '-O', str(output_path)]
                subprocess.run(cmd, check=True)
            
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Conversion to MOL2 failed: {e}")
            raise
    
    def _minimize_ligand(self, mol2_file: str, output_dir: str) -> str:
        """Minimize ligand structure using RDKit."""
        mol2_path = Path(mol2_file)
        output_path = Path(output_dir) / f"{mol2_path.stem}_minimized.mol2"
        
        try:
            # Load molecule
            mol = Chem.MolFromMol2File(str(mol2_path))
            if mol is None:
                # Fallback: try to load with Open Babel and convert
                logger.warning("RDKit failed to load MOL2, using original file")
                return mol2_file
            
            # Add hydrogens if not present
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # Minimize with MMFF
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
            except:
                logger.warning("MMFF minimization failed, using UFF")
                AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
            
            # Write minimized structure
            # Note: RDKit doesn't directly write MOL2, so we'll use SDF and convert
            sdf_file = output_path.with_suffix('.sdf')
            writer = Chem.SDWriter(str(sdf_file))
            writer.write(mol)
            writer.close()
            
            # Convert to MOL2
            cmd = ['obabel', str(sdf_file), '-O', str(output_path)]
            result = subprocess.run(cmd, capture_output=True)
            
            if result.returncode == 0:
                return str(output_path)
            else:
                logger.warning("MOL2 conversion failed, using original")
                return mol2_file
                
        except Exception as e:
            logger.warning(f"Ligand minimization failed: {e}, using original")
            return mol2_file
    
    def _convert_to_pdbqt_ligand(self, mol2_file: str, output_pdbqt: str):
        """Convert ligand to PDBQT format using meeko or prepare_ligand4.py."""
        try:
            # Try using meeko first (more modern)
            from meeko import MoleculePreparation
            from rdkit import Chem
            
            mol = Chem.MolFromMol2File(mol2_file)
            if mol is not None:
                preparator = MoleculePreparation()
                preparator.prepare(mol)
                preparator.write_pdbqt_file(output_pdbqt)
                logger.info("Used meeko for PDBQT conversion")
                return
                
        except ImportError:
            logger.info("Meeko not available, trying prepare_ligand4.py")
        except Exception as e:
            logger.warning(f"Meeko failed: {e}, trying prepare_ligand4.py")
        
        # Fallback to prepare_ligand4.py from AutoDockTools
        try:
            cmd = [
                'prepare_ligand4.py',
                '-l', mol2_file,
                '-o', output_pdbqt,
                '-A', 'hydrogens',
                '-U', 'nphs_lps'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
                
            logger.info("Used prepare_ligand4.py for PDBQT conversion")
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Final fallback: use Open Babel
            logger.info("prepare_ligand4.py not available, using Open Babel")
            cmd = ['obabel', mol2_file, '-O', output_pdbqt]
            subprocess.run(cmd, check=True)
    
    def _clean_protein(self, input_pdb: str, output_pdb: str):
        """Clean protein structure by removing waters and adding hydrogens."""
        try:
            # Use BioPython to clean the structure
            structure = self.parser.get_structure('protein', input_pdb)
            
            # Remove water molecules
            for model in structure:
                chains_to_remove = []
                for chain in model:
                    residues_to_remove = []
                    for residue in chain:
                        if residue.get_resname() in ['HOH', 'WAT', 'H2O']:
                            residues_to_remove.append(residue.get_id())
                    
                    for res_id in residues_to_remove:
                        chain.detach_child(res_id)
                    
                    if len(chain) == 0:
                        chains_to_remove.append(chain.get_id())
                
                for chain_id in chains_to_remove:
                    model.detach_child(chain_id)
            
            # Save cleaned structure
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(output_pdb))
            
            logger.info(f"🧹 Cleaned protein structure: {output_pdb}")
            
        except Exception as e:
            logger.warning(f"BioPython cleaning failed: {e}, copying original")
            import shutil
            shutil.copy2(input_pdb, output_pdb)
    
    def _convert_to_pdbqt_receptor(self, input_pdb: str, output_pdbqt: str):
        """Convert protein to PDBQT format using prepare_receptor4.py."""
        try:
            cmd = [
                'prepare_receptor4.py',
                '-r', input_pdb,
                '-o', output_pdbqt,
                '-A', 'hydrogens',
                '-U', 'nphs_lps_waters_nonstdres'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
                
            logger.info("Used prepare_receptor4.py for PDBQT conversion")
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Fallback: use Open Babel with basic conversion
            logger.warning("prepare_receptor4.py not available, using Open Babel fallback")
            cmd = ['obabel', input_pdb, '-O', output_pdbqt, '-p', '7.4']
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError:
                # Final fallback: manual conversion
                logger.warning("Open Babel failed, using manual conversion")
                self._manual_pdb_to_pdbqt(input_pdb, output_pdbqt)
    
    def _manual_pdb_to_pdbqt(self, input_pdb: str, output_pdbqt: str):
        """Manual PDB to PDBQT conversion with proper atom types."""
        # Map PDB element to AutoDock atom types
        atom_type_map = {
            'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'P': 'P',
            'H': 'H', 'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I',
            'FE': 'Fe', 'ZN': 'Zn', 'MG': 'Mg', 'CA': 'Ca', 'MN': 'Mn'
        }
        
        with open(input_pdb, 'r') as infile, open(output_pdbqt, 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    # Extract element from PDB line (columns 77-78)
                    element = line[76:78].strip() if len(line) > 77 else ''
                    if not element:
                        # Try to guess from atom name (columns 13-16)
                        atom_name = line[12:16].strip()
                        element = atom_name[0] if atom_name else 'C'
                    
                    # Map to AutoDock type
                    ad_type = atom_type_map.get(element.upper(), 'C')
                    
                    # Format PDBQT line: PDB line + charge + AutoDock type  
                    # Ensure proper spacing for PDBQT format
                    base_line = line[:66].ljust(66)
                    pdbqt_line = f"{base_line}  0.000  0.000    {ad_type}\n"
                    outfile.write(pdbqt_line)
                elif line.startswith(('REMARK', 'HEADER', 'TITLE')):
                    outfile.write(line)
        
        logger.warning("Used basic manual conversion - charges may not be accurate")
    
    def _extract_coordinates_from_pdbqt(self, pdbqt_file: str) -> list:
        """Extract atomic coordinates from PDBQT file."""
        coords = []
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
        
        return coords 
    def _convert_ligand_simple(self, ligand_file: str, output_dir: str) -> str:
        """تبدیل ساده لیگاند به PDBQT بدون تمیزکاری حرفه‌ای"""
        try:
            output_file = Path(output_dir) / f"{Path(ligand_file).stem}.pdbqt"
            
            # استفاده از OpenBabel برای تبدیل
            cmd = f"obabel {ligand_file} -O {output_file} -h"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and output_file.exists():
                logger.success(f"✅ لیگاند تبدیل شد: {output_file}")
                return str(output_file)
            else:
                logger.error("❌ تبدیل لیگاند ناموفق")
                return None
                
        except Exception as e:
            logger.error(f"❌ خطا در تبدیل لیگاند: {str(e)}")
            return None
