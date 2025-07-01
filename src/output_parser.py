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
from loguru import logger

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
            # Parse binding affinities from log file
            affinities = self._parse_binding_affinities(log_file)
            
            # Extract individual poses from PDBQT
            poses = self._extract_poses_from_pdbqt(output_pdbqt, output_dir)
            
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
            
            logger.success(f"✅ Parsed {len(poses)} poses with best affinity: {results['best_affinity']:.2f} kcal/mol")
            
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
        
        if not SCIPY_AVAILABLE:
            logger.warning("⚠️ SciPy not available, limited interaction analysis")
            return analysis
            
        try:
            # Extract coordinates
            protein_coords, protein_residues = self._extract_protein_residues(protein_file)
            ligand_coords = self._extract_coordinates_simple(ligand_file)
            
            if not protein_coords or not ligand_coords:
                logger.warning("⚠️ Could not extract coordinates for interaction analysis")
                return analysis
            
            # Calculate distance matrix
            distances = cdist(ligand_coords, protein_coords)
            
            # Find close contacts
            close_contacts = np.where(distances < self.interaction_cutoffs['hydrophobic'])
            
            # Analyze binding site residues
            unique_residue_indices = np.unique(close_contacts[1])
            binding_site_residues = []
            
            for idx in unique_residue_indices:
                if idx < len(protein_residues):
                    residue_info = protein_residues[idx]
                    min_distance = np.min(distances[:, idx])
                    
                    binding_site_residues.append({
                        'residue': residue_info['name'],
                        'chain': residue_info['chain'],
                        'number': residue_info['number'],
                        'min_distance': float(min_distance),
                        'interaction_type': self._classify_interaction(
                            residue_info['name'], min_distance
                        )
                    })
            
            analysis['binding_site_residues'] = binding_site_residues
            
            # Calculate geometric properties
            analysis['geometric_properties'] = {
                'binding_site_volume': self._calculate_binding_site_volume(
                    protein_coords, ligand_coords
                ),
                'buried_surface_area': self._estimate_buried_surface_area(
                    protein_coords, ligand_coords
                ),
                'ligand_protein_contacts': len(close_contacts[0])
            }
            
        except Exception as e:
            logger.warning(f"⚠️ Enhanced interaction analysis failed: {e}")
            
        return analysis

    def _extract_protein_residues(self, protein_file: str) -> Tuple[List, List]:
        """Extract protein coordinates and residue information."""
        coords = []
        residues = []
        
        try:
            with open(protein_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                            
                            residue_info = {
                                'name': line[17:20].strip(),
                                'chain': line[21:22].strip(),
                                'number': int(line[22:26].strip())
                            }
                            residues.append(residue_info)
                            
                        except (ValueError, IndexError):
                            continue
                            
        except Exception as e:
            logger.warning(f"⚠️ Failed to extract protein residues: {e}")
            
        return coords, residues

    def _extract_coordinates_simple(self, file_path: str) -> List:
        """Simple coordinate extraction from PDB-like files."""
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
                            
        except Exception as e:
            logger.warning(f"⚠️ Failed to extract coordinates: {e}")
            
        return coords

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

    def _parse_binding_affinities(self, log_file: str) -> List[Dict]:
        """Parse binding affinities from Vina log file."""
        affinities = []
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            # Look for the results table
            table_pattern = r'mode.*?affinity.*?kcal/mol.*?\n.*?-+\n(.*?)(?=\n\n|\Z)'
            table_match = re.search(table_pattern, content, re.DOTALL)
            
            if table_match:
                table_content = table_match.group(1)
                
                # Parse each result line
                for line in table_content.strip().split('\n'):
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = [p.strip() for p in line.split('|')]
                    
                    if len(parts) >= 2:
                        try:
                            mode = int(parts[0])
                            affinity = float(parts[1])
                            
                            affinities.append({
                                'mode': mode,
                                'affinity': affinity,
                                'rmsd_lb': float(parts[2]) if len(parts) > 2 else None,
                                'rmsd_ub': float(parts[3]) if len(parts) > 3 else None
                            })
                        except ValueError:
                            continue
            
            # Sort by affinity (best first)
            affinities.sort(key=lambda x: x['affinity'])
            
            logger.info(f"📈 Parsed {len(affinities)} binding affinities")
            
            return affinities
            
        except Exception as e:
            logger.error(f"Failed to parse binding affinities: {e}")
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
