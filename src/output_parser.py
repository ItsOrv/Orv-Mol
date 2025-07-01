"""
Output parser module for processing AutoDock Vina results.
Extracts binding poses, affinity scores, and converts to various formats.
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
from .logger_config import logger

try:
    from scipy.spatial.distance import cdist
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    logger.warning("SciPy not available, some analysis features will be limited")

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available, some analysis features will be limited")


class OutputParser:
    """Handles parsing and processing of AutoDock Vina output files."""
    
    def __init__(self):
        """Initialize output parser with analysis capabilities."""
        self.interaction_cutoffs = {
            'hydrogen_bond': 3.5,  # Å
            'hydrophobic': 4.0,    # Å
            'electrostatic': 5.0,  # Å
            'van_der_waals': 4.5   # Å
        }
        self.amino_acid_properties = self._load_amino_acid_properties()
    
    def parse_vina_output(
        self,
        log_file: str,
        output_pdbqt: str,
        output_dir: str
    ) -> Dict:
        """
        Parse AutoDock Vina output and extract results.
        """
        logger.info("📊 Parsing AutoDock Vina results...")
        
        try:
            # Parse binding affinities from PDBQT file (modern Vina versions)
            affinities = self._parse_binding_affinities(output_pdbqt)
            
            # Extract individual poses from PDBQT
            poses = self._extract_poses_from_pdbqt(output_pdbqt, output_dir)
            
            # Add affinity data to poses if available
            if affinities and len(affinities) == len(poses):
                for i, pose in enumerate(poses):
                    if i < len(affinities):
                        pose['affinity'] = affinities[i]['affinity']
                        pose['rmsd_lb'] = affinities[i]['rmsd_lb']
                        pose['rmsd_ub'] = affinities[i]['rmsd_ub']
            
            # Convert best pose to different formats
            best_pose_files = self._convert_best_pose(poses[0] if poses else None, output_dir)
            
            # Generate summary statistics
            summary = self._generate_summary(affinities, poses)
            
            results = {
                'affinities': affinities,
                'poses': poses,
                'best_affinity': affinities[0]['affinity'] if affinities else None,
                'best_pose_pdbqt': poses[0]['file'] if poses else None,
                'best_pose_pdb': best_pose_files.get('pdb'),
                'best_pose_sdf': best_pose_files.get('sdf'),
                'summary': summary,
                'num_poses': len(poses)
            }
            
            # Save results to CSV
            self._save_results_csv(results, output_dir)
            
            # Save detailed JSON results
            self._save_results_json(results, output_dir)
            
            # Success message with safe formatting
            if results['best_affinity'] is not None:
            logger.success(f"✅ Parsed {len(poses)} poses with best affinity: {results['best_affinity']:.2f} kcal/mol")
            else:
                logger.success(f"✅ Parsed {len(poses)} poses (no affinity data found)")
            
            return results
            
        except Exception as e:
            logger.error(f"❌ Failed to parse results: {e}")
            raise
    
    def analyze_binding_interactions(self, detailed_results: Dict, protein_file: str = None, ligand_file: str = None) -> Dict:
        """
        Analyze protein-ligand binding interactions in detail.
        
        Args:
            detailed_results: Results from parse_vina_output
            protein_file: Path to protein PDB file (optional, for enhanced analysis)
            ligand_file: Path to ligand file (optional, for enhanced analysis)
            
        Returns:
            Dictionary with detailed interaction analysis
        """
        logger.info("🔬 Analyzing binding interactions...")
        
        analysis = {
            'interaction_summary': {},
            'binding_site_residues': [],
            'interaction_types': {},
            'geometric_properties': {},
            'pharmacological_analysis': {}
        }
        
        try:
            # Basic analysis from poses
            if detailed_results.get('best_pose_pdb'):
                pose_analysis = self._analyze_pose_interactions(detailed_results['best_pose_pdb'])
                analysis['interaction_summary'] = pose_analysis
            
            # Enhanced analysis if protein and ligand files provided
            if protein_file and ligand_file:
                enhanced_analysis = self._enhanced_interaction_analysis(
                    protein_file, detailed_results.get('best_pose_pdb', ligand_file)
                )
                analysis.update(enhanced_analysis)
            
            # Calculate binding site properties
            if detailed_results.get('affinities'):
                analysis['binding_affinity_analysis'] = self._analyze_binding_affinity_distribution(
                    detailed_results['affinities']
                )
            
            # Pharmacological predictions
            if RDKIT_AVAILABLE and ligand_file:
                analysis['pharmacological_analysis'] = self._analyze_drug_properties(ligand_file)
            
            logger.success("✅ Interaction analysis completed")
            return analysis
            
        except Exception as e:
            logger.error(f"❌ Failed to analyze interactions: {e}")
            return analysis

    def _load_amino_acid_properties(self) -> Dict:
        """Load amino acid properties for interaction analysis."""
        return {
            'hydrophobic': ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'],
            'polar': ['SER', 'THR', 'TYR', 'ASN', 'GLN', 'CYS'],
            'charged_positive': ['LYS', 'ARG', 'HIS'],
            'charged_negative': ['ASP', 'GLU'],
            'aromatic': ['PHE', 'TYR', 'TRP', 'HIS'],
            'hydrogen_bond_donors': ['SER', 'THR', 'TYR', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'TRP'],
            'hydrogen_bond_acceptors': ['ASP', 'GLU', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'HIS']
        }

    def _analyze_pose_interactions(self, pose_file: str) -> Dict:
        """Analyze interactions from a single pose file."""
        analysis = {
            'total_atoms': 0,
            'heavy_atoms': 0,
            'coordinates_center': [0, 0, 0],
            'bounding_box': {}
        }
        
        try:
            coords = []
            heavy_atom_coords = []
            
            with open(pose_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        analysis['total_atoms'] += 1
                        
                        # Extract coordinates
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                            
                            # Check if heavy atom (not hydrogen)
                            atom_name = line[12:16].strip()
                            if not atom_name.startswith('H'):
                                analysis['heavy_atoms'] += 1
                                heavy_atom_coords.append([x, y, z])
                                
                        except (ValueError, IndexError):
                            continue
            
            if coords:
                coords_array = np.array(coords)
                analysis['coordinates_center'] = list(np.mean(coords_array, axis=0))
                
                # Calculate bounding box
                min_coords = np.min(coords_array, axis=0)
                max_coords = np.max(coords_array, axis=0)
                analysis['bounding_box'] = {
                    'min': list(min_coords),
                    'max': list(max_coords),
                    'size': list(max_coords - min_coords)
                }
                
                # Calculate molecular volume approximation
                if len(coords) > 0:
                    box_volume = np.prod(max_coords - min_coords)
                    analysis['approximate_volume'] = float(box_volume)
                
        except Exception as e:
            logger.warning(f"⚠️ Pose analysis failed: {e}")
            
        return analysis

    def _enhanced_interaction_analysis(self, protein_file: str, ligand_file: str) -> Dict:
        """Enhanced interaction analysis with protein-ligand contact analysis."""
        analysis = {
            'binding_site_residues': [],
            'interaction_types': {
                'hydrogen_bonds': [],
                'hydrophobic_contacts': [],
                'electrostatic_interactions': [],
                'van_der_waals': []
            },
            'geometric_properties': {}
        }
        
        try:
            # Extract coordinates and residue information
            protein_coords, protein_residues = self._extract_protein_residues(protein_file)
            ligand_coords = self._extract_coordinates_simple(ligand_file)
            
            if not protein_coords or not ligand_coords:
                logger.warning("⚠️ Could not extract coordinates for interaction analysis")
                # Try alternative extraction method
                protein_coords = self._extract_coordinates_simple(protein_file)
                if not protein_coords or not ligand_coords:
                    return analysis
                # Create dummy residue information
                protein_residues = [{'name': 'UNK', 'number': i+1, 'chain': 'A'} for i in range(len(protein_coords))]
            
            if not SCIPY_AVAILABLE:
                logger.warning("⚠️ SciPy not available, using simplified distance calculation")
                # Simple distance calculation without scipy
                binding_site_residues = self._find_binding_site_residues_simple(
                    protein_coords, ligand_coords, protein_residues
                )
                analysis['binding_site_residues'] = binding_site_residues
            else:
                # Calculate distance matrix using scipy
                distances = cdist(ligand_coords, protein_coords)
                
                # Find close contacts (within 5 Å)
                close_contacts = np.where(distances < 5.0)
                
                # Analyze binding site residues
                unique_residue_indices = np.unique(close_contacts[1])
                binding_site_residues = []
                
                for res_idx in unique_residue_indices:
                    if res_idx < len(protein_residues):
                        residue = protein_residues[res_idx]
                        min_distance = np.min(distances[:, res_idx])
                        
                        binding_site_residues.append({
                            'residue_name': residue.get('name', 'UNK'),
                            'residue_number': residue.get('number', res_idx + 1),
                            'chain': residue.get('chain', 'A'),
                            'min_distance': float(min_distance),
                            'interaction_type': self._classify_interaction(residue.get('name', 'UNK'), min_distance)
                        })
                
                analysis['binding_site_residues'] = binding_site_residues
                
                # Calculate geometric properties
                if len(ligand_coords) > 0:
                    ligand_center = np.mean(ligand_coords, axis=0)
                    analysis['geometric_properties'] = {
                        'ligand_center': ligand_center.tolist(),
                        'binding_site_volume': self._calculate_binding_site_volume(protein_coords, ligand_coords),
                        'buried_surface_area': self._estimate_buried_surface_area(protein_coords, ligand_coords),
                        'num_binding_site_residues': len(binding_site_residues)
                    }
            
            logger.info(f"🔍 Found {len(analysis['binding_site_residues'])} binding site residues")
            return analysis
            
        except Exception as e:
            logger.error(f"❌ Failed to analyze interactions: {e}")
            return analysis

    def _find_binding_site_residues_simple(self, protein_coords, ligand_coords, protein_residues):
        """Simple binding site residue finding without scipy."""
        binding_site_residues = []
        
        try:
            cutoff = 5.0  # 5 Angstrom cutoff
            
            for i, protein_coord in enumerate(protein_coords):
                min_distance = float('inf')
                
                # Calculate minimum distance to any ligand atom
                for ligand_coord in ligand_coords:
                    distance = ((protein_coord[0] - ligand_coord[0])**2 + 
                              (protein_coord[1] - ligand_coord[1])**2 + 
                              (protein_coord[2] - ligand_coord[2])**2)**0.5
                    min_distance = min(min_distance, distance)
                
                # If within cutoff, add to binding site
                if min_distance <= cutoff:
                    residue = protein_residues[i] if i < len(protein_residues) else {'name': 'UNK', 'number': i+1, 'chain': 'A'}
                    binding_site_residues.append({
                        'residue_name': residue.get('name', 'UNK'),
                        'residue_number': residue.get('number', i + 1),
                        'chain': residue.get('chain', 'A'),
                        'min_distance': min_distance,
                        'interaction_type': self._classify_interaction(residue.get('name', 'UNK'), min_distance)
                    })
            
            return binding_site_residues
            
        except Exception as e:
            logger.warning(f"⚠️ Simple binding site analysis failed: {e}")
            return []

    def _extract_protein_residues(self, protein_file: str) -> Tuple[List, List]:
        """Extract protein coordinates and residue information from PDB file."""
        coords = []
        residues = []
        
        try:
            with open(protein_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        try:
                            # Extract coordinates
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                            
                            # Extract residue information
                            res_name = line[17:20].strip()
                            res_number = int(line[22:26].strip())
                            chain = line[21:22].strip()
                            
                            residues.append({
                                'name': res_name,
                                'number': res_number,
                                'chain': chain
                            })
                            
                        except (ValueError, IndexError):
                            continue
            
            logger.info(f"📍 Extracted {len(coords)} protein atoms from {len(set(r['number'] for r in residues))} residues")
            return coords, residues
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to extract protein residues: {e}")
            return [], []

    def _extract_coordinates_simple(self, file_path: str) -> List:
        """Extract coordinates from PDB/PDBQT file."""
        coords = []
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                        except (ValueError, IndexError):
                            continue
            
            logger.info(f"📍 Extracted {len(coords)} coordinates from {file_path}")
            return coords
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to extract coordinates: {e}")
            return []

    def _classify_interaction(self, residue_name: str, distance: float) -> str:
        """Classify interaction type based on residue and distance."""
        if distance <= self.interaction_cutoffs['hydrogen_bond']:
            if residue_name in self.amino_acid_properties['hydrogen_bond_donors'] or \
               residue_name in self.amino_acid_properties['hydrogen_bond_acceptors']:
                return 'hydrogen_bond'
        
        if distance <= self.interaction_cutoffs['electrostatic']:
            if residue_name in self.amino_acid_properties['charged_positive'] or \
               residue_name in self.amino_acid_properties['charged_negative']:
                return 'electrostatic'
        
        if distance <= self.interaction_cutoffs['hydrophobic']:
            if residue_name in self.amino_acid_properties['hydrophobic']:
                return 'hydrophobic'
        
        if distance <= self.interaction_cutoffs['van_der_waals']:
            return 'van_der_waals'
        
        return 'weak_contact'

    def _calculate_binding_site_volume(self, protein_coords: List, ligand_coords: List) -> float:
        """Estimate binding site volume."""
        if not protein_coords or not ligand_coords:
            return 0.0
        
        try:
            # Simple approximation using ligand bounding box
            ligand_array = np.array(ligand_coords)
            min_coords = np.min(ligand_array, axis=0)
            max_coords = np.max(ligand_array, axis=0)
            
            # Add padding for binding site
            padding = 2.0  # Angstroms
            volume = np.prod(max_coords - min_coords + 2 * padding)
            
            return float(volume)
            
        except Exception:
            return 0.0

    def _estimate_buried_surface_area(self, protein_coords: List, ligand_coords: List) -> float:
        """Estimate buried surface area (simplified calculation)."""
        if not protein_coords or not ligand_coords or not SCIPY_AVAILABLE:
            return 0.0
        
        try:
            # Calculate contacts within 4 Å
            distances = cdist(ligand_coords, protein_coords)
            close_contacts = np.sum(distances < 4.0)
            
            # Rough approximation: 15 Ų per close contact
            estimated_bsa = close_contacts * 15.0
            
            return float(estimated_bsa)
            
        except Exception:
            return 0.0

    def _analyze_binding_affinity_distribution(self, affinities: List[Dict]) -> Dict:
        """Analyze binding affinity distribution across poses."""
        if not affinities:
            return {}
        
        energies = [a['affinity'] for a in affinities]
        
        analysis = {
            'best_affinity': min(energies),
            'worst_affinity': max(energies),
            'mean_affinity': np.mean(energies),
            'std_affinity': np.std(energies),
            'median_affinity': np.median(energies),
            'poses_within_1_kcal': len([e for e in energies if e <= min(energies) + 1.0]),
            'poses_within_2_kcal': len([e for e in energies if e <= min(energies) + 2.0]),
            'poses_within_3_kcal': len([e for e in energies if e <= min(energies) + 3.0]),
            'energy_range': max(energies) - min(energies)
        }
        
        # Classify binding strength
        best_energy = min(energies)
        if best_energy <= -10.0:
            analysis['binding_strength'] = 'very_strong'
        elif best_energy <= -8.0:
            analysis['binding_strength'] = 'strong'
        elif best_energy <= -6.0:
            analysis['binding_strength'] = 'moderate'
        elif best_energy <= -4.0:
            analysis['binding_strength'] = 'weak'
        else:
            analysis['binding_strength'] = 'very_weak'
        
        return analysis

    def _analyze_drug_properties(self, ligand_file: str) -> Dict:
        """Analyze drug-like properties of the ligand."""
        if not RDKIT_AVAILABLE:
            return {}
        
        analysis = {}
        
        try:
            # Load molecule
            if ligand_file.endswith('.sdf'):
                mol = Chem.SDMolSupplier(ligand_file)[0]
            elif ligand_file.endswith('.pdb'):
                mol = Chem.MolFromPDBFile(ligand_file)
            else:
                mol = Chem.MolFromMolFile(ligand_file)
            
            if mol is None:
                return {}
            
            # Calculate molecular descriptors
            analysis['molecular_weight'] = Descriptors.MolWt(mol)
            analysis['logp'] = Descriptors.MolLogP(mol)
            analysis['hbd'] = Descriptors.NumHDonors(mol)
            analysis['hba'] = Descriptors.NumHAcceptors(mol)
            analysis['rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
            analysis['tpsa'] = Descriptors.TPSA(mol)
            analysis['heavy_atoms'] = Descriptors.HeavyAtomCount(mol)
            
            # Lipinski's Rule of Five
            lipinski_violations = 0
            if analysis['molecular_weight'] > 500:
                lipinski_violations += 1
            if analysis['logp'] > 5:
                lipinski_violations += 1
            if analysis['hbd'] > 5:
                lipinski_violations += 1
            if analysis['hba'] > 10:
                lipinski_violations += 1
            
            analysis['lipinski_violations'] = lipinski_violations
            analysis['lipinski_compliant'] = lipinski_violations <= 1
            
            # Drug-likeness score
            drug_score = 10.0
            drug_score -= lipinski_violations * 2
            if analysis['rotatable_bonds'] > 10:
                drug_score -= 1
            if analysis['tpsa'] > 140:
                drug_score -= 1
            
            analysis['drug_likeness_score'] = max(drug_score, 0.0)
            
        except Exception as e:
            logger.warning(f"⚠️ Drug property analysis failed: {e}")
            
        return analysis

    def _parse_binding_affinities(self, pdbqt_file: str) -> List[Dict]:
        """Parse binding affinities from PDBQT file using REMARK VINA RESULT format."""
        affinities = []
        
        try:
            with open(pdbqt_file, 'r') as f:
                content = f.read()
            
            # Find all REMARK VINA RESULT lines
            vina_results = re.findall(r'REMARK VINA RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', content)
            
            for i, result in enumerate(vina_results, 1):
                try:
                    affinity = float(result[0])
                    rmsd_lb = float(result[1])
                    rmsd_ub = float(result[2])
                    
                    affinities.append({
                        'mode': i,
                        'affinity': affinity,
                        'rmsd_lb': rmsd_lb,
                        'rmsd_ub': rmsd_ub
                    })
                except ValueError:
                            continue
            
            # Sort by affinity (best first)
            if affinities:
            affinities.sort(key=lambda x: x['affinity'])
            
            logger.info(f"📈 Parsed {len(affinities)} binding affinities from PDBQT")
            
            return affinities
            
        except Exception as e:
            logger.error(f"Failed to parse affinities from PDBQT: {e}")
            return []
    
    def _extract_poses_from_pdbqt(self, output_pdbqt: str, output_dir: str) -> List[Dict]:
        """Extract individual poses from multi-model PDBQT file."""
        poses = []
        output_path = Path(output_dir)
        
        try:
            with open(output_pdbqt, 'r') as f:
                content = f.read()
            
            # Split by MODEL/ENDMDL keywords
            models = re.split(r'MODEL\s+\d+', content)[1:]  # Skip first empty part
            
            for i, model_content in enumerate(models, 1):
                # Remove ENDMDL if present
                model_content = re.sub(r'ENDMDL.*', '', model_content, flags=re.DOTALL)
                
                # Create filename for this pose
                pose_file = output_path / f"pose_{i:02d}.pdbqt"
                
                # Write pose to file
                with open(pose_file, 'w') as f:
                    f.write(f"MODEL {i}\n")
                    f.write(model_content.strip())
                    f.write("\nENDMDL\n")
                
                poses.append({
                    'mode': i,
                    'file': str(pose_file)
                })
            
            logger.info(f"📁 Extracted {len(poses)} poses")
            
            return poses
            
        except Exception as e:
            logger.error(f"Failed to extract poses: {e}")
            return []
    
    def _convert_best_pose(self, best_pose: Optional[Dict], output_dir: str) -> Dict[str, str]:
        """Convert best pose to PDB and SDF formats."""
        if not best_pose:
            return {}
        
        output_path = Path(output_dir)
        converted_files = {}
        
        try:
            # Convert to PDB using Open Babel
            pdb_file = output_path / "best_pose.pdb"
            self._convert_pdbqt_to_pdb(best_pose['file'], pdb_file)
            converted_files['pdb'] = str(pdb_file)
            
            # Convert to SDF using Open Babel
            sdf_file = output_path / "best_pose.sdf"
            self._convert_pdbqt_to_sdf(best_pose['file'], sdf_file)
            converted_files['sdf'] = str(sdf_file)
            
            logger.info(f"💾 Converted best pose to PDB and SDF formats")
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to convert best pose: {e}")
        
        return converted_files
    
    def _convert_pdbqt_to_pdb(self, pdbqt_file: str, pdb_file: Path):
        """Convert PDBQT to PDB format."""
        try:
            cmd = ['obabel', pdbqt_file, '-O', str(pdb_file)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
                
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("Open Babel not available, using manual conversion")
            self._manual_pdbqt_to_pdb(pdbqt_file, pdb_file)
    
    def _convert_pdbqt_to_sdf(self, pdbqt_file: str, sdf_file: Path):
        """Convert PDBQT to SDF format."""
        try:
            cmd = ['obabel', pdbqt_file, '-O', str(sdf_file)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
                
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("Open Babel not available for SDF conversion, skipping")
    
    def _manual_pdbqt_to_pdb(self, pdbqt_file: str, pdb_file: Path):
        """Manual PDBQT to PDB conversion."""
        with open(pdbqt_file, 'r') as infile, open(pdb_file, 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    pdb_line = line[:66] + '\n'
                    outfile.write(pdb_line)
                elif line.startswith(('MODEL', 'ENDMDL', 'REMARK', 'HEADER')):
                    outfile.write(line)
    
    def _generate_summary(self, affinities: List[Dict], poses: List[Dict]) -> Dict:
        """Generate summary statistics from results."""
        if not affinities:
            return {}
        
        affinity_values = [a['affinity'] for a in affinities]
        
        summary = {
            'best_affinity': min(affinity_values),
            'worst_affinity': max(affinity_values),
            'mean_affinity': sum(affinity_values) / len(affinity_values),
            'affinity_range': max(affinity_values) - min(affinity_values),
            'num_poses': len(poses),
            'poses_within_1_kcal': len([a for a in affinity_values if a <= min(affinity_values) + 1.0]),
            'poses_within_2_kcal': len([a for a in affinity_values if a <= min(affinity_values) + 2.0])
        }
        
        return summary
    
    def _save_results_csv(self, results: Dict, output_dir: str):
        """Save results to CSV file."""
        csv_file = Path(output_dir) / "docking_results.csv"
        
        try:
            if results['affinities']:
                df = pd.DataFrame(results['affinities'])
                df.to_csv(csv_file, index=False)
                logger.info(f"💾 Saved results to: {csv_file}")
        except Exception as e:
            logger.warning(f"⚠️ Failed to save CSV: {e}")

    def _save_results_json(self, results: Dict, output_dir: str):
        """Save detailed results to JSON file."""
        import json
        
        json_file = Path(output_dir) / "docking_results.json"
        
        try:
            # Convert numpy types to native Python types for JSON serialization
            def convert_numpy_types(obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                return obj
            
            # Deep copy and convert results
            json_results = {}
            for key, value in results.items():
                if isinstance(value, dict):
                    json_results[key] = {k: convert_numpy_types(v) for k, v in value.items()}
                elif isinstance(value, list):
                    json_results[key] = [convert_numpy_types(item) if not isinstance(item, dict) 
                                       else {k: convert_numpy_types(v) for k, v in item.items()} 
                                       for item in value]
                else:
                    json_results[key] = convert_numpy_types(value)
            
            with open(json_file, 'w') as f:
                json.dump(json_results, f, indent=2, ensure_ascii=False)
                
            logger.info(f"💾 Saved detailed results to: {json_file}")
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to save JSON: {e}")
