"""
Output parser module for processing AutoDock Vina results.
Extracts binding poses, affinity scores, and converts to various formats.
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from loguru import logger


class OutputParser:
    """Handles parsing and processing of AutoDock Vina output files."""
    
    def __init__(self):
        pass
    
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
            
            logger.success(f"✅ Parsed {len(poses)} poses with best affinity: {results['best_affinity']:.2f} kcal/mol")
            
            return results
            
        except Exception as e:
            logger.error(f"❌ Failed to parse results: {e}")
            raise
    
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
