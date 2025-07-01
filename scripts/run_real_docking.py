#!/usr/bin/env python3
"""
Real Docking Runner with Professional Cleaning System
"""

import os
import sys
import warnings
from pathlib import Path

# Suppress warnings that cause interruptions
warnings.filterwarnings("ignore")
os.environ['PYTHONWARNINGS'] = 'ignore'

# Add src to path
sys.path.append('src')

from loguru import logger
from src.preprocessing import PreprocessingEngine
from src.docking import DockingEngine
from src.utils import create_output_structure

def run_professional_docking():
    """Run professional docking with error management"""
    
    logger.info("Starting professional molecular docking")
    logger.info("=" * 50)
    
    # Input files
    protein_file = "data/input/protein.pdb"
    ligand_file = "data/input/ligand_small.pdb"  # Use small ligand
    output_dir = "results/real_docking_results"
    
    # Check file existence
    if not os.path.exists(protein_file):
        logger.error(f"❌ Protein file not found: {protein_file}")
        return False
        
    if not os.path.exists(ligand_file):
        logger.error(f"❌ Ligand file not found: {ligand_file}")
        return False
    
    logger.info(f"Protein: {protein_file}")
    logger.info(f"Ligand: {ligand_file}")
    logger.info(f"Output: {output_dir}")
    
    try:
        # Create output structure
        logger.info("Creating directory structure...")
        output_paths = create_output_structure(output_dir)
        logger.success("✅ Directory structure created")
        
        # Initialize professional cleaning system
        logger.info("Initializing professional cleaning system...")
        preprocessor = PreprocessingEngine(
            professional_mode=True,
            ph=7.4
        )
        logger.success("✅ Cleaning system ready")
        
        # Prepare protein
        logger.info("Preparing protein...")
        protein_pdbqt = preprocessor.prepare_receptor(
            protein_file, output_paths['temp']
        )
        logger.success(f"✅ Protein prepared: {protein_pdbqt}")
        
        # Prepare ligand
        logger.info("Preparing ligand...")
        ligand_pdbqt = preprocessor.prepare_ligand(
            ligand_file, output_paths['temp']
        )
        
        if ligand_pdbqt:
            logger.success(f"✅ Ligand prepared: {ligand_pdbqt}")
        else:
            logger.warning("⚠️ Professional ligand cleaning failed, using standard method...")
            # Simple ligand conversion
            ligand_pdbqt = preprocessor._convert_ligand_simple(ligand_file, output_paths['temp'])
            logger.info(f"Ligand prepared with standard method: {ligand_pdbqt}")
        
        # Calculate docking box
        logger.info("Calculating docking box...")
        box_params = preprocessor.calculate_blind_docking_box(protein_pdbqt)
        logger.info(f"Box center: ({box_params['center_x']:.1f}, {box_params['center_y']:.1f}, {box_params['center_z']:.1f})")
        logger.info(f"Box size: ({box_params['size_x']:.1f}, {box_params['size_y']:.1f}, {box_params['size_z']:.1f})")
        logger.success("✅ Docking box calculated")
        
        # Optimize parameters
        logger.info("Optimizing docking parameters...")
        optimized_params = preprocessor.get_optimized_docking_parameters()
        logger.info(f"Optimized parameters: {optimized_params}")
        logger.success("✅ Parameters optimized")
        
        # Run docking
        logger.info("Starting docking with AutoDock Vina...")
        docker = DockingEngine()
        
        # Run docking with optimized parameters
        docking_results = docker.run_vina_docking(
            receptor=protein_pdbqt,
            ligand=ligand_pdbqt,
            box_params=box_params,
            output_dir=output_paths['poses'],
            exhaustiveness=optimized_params.get('exhaustiveness', 8),
            num_modes=optimized_params.get('num_modes', 9),
            energy_range=optimized_params.get('energy_range', 3.0)
        )
        
        # Check results
        if docking_results.get('success', False):
            logger.success("✅ Docking completed successfully!")
            
            # Display results
            best_affinity = docking_results.get('best_affinity', 'N/A')
            poses = docking_results.get('poses', [])
            
            logger.info(f"Best binding energy: {best_affinity} kcal/mol")
            logger.info(f"Total poses: {len(poses)}")
            
            if poses:
                logger.info("Top 3 poses:")
                for i, pose in enumerate(poses[:3], 1):
                    affinity = pose.get('affinity', 'N/A')
                    logger.info(f"  Pose {i}: {affinity} kcal/mol")
            
            # Cleaning report
            logger.info("Professional cleaning report:")
            summary = preprocessor.get_cleaning_summary()
            
            if summary.get('protein_analysis'):
                protein_info = summary['protein_analysis']
                logger.info(f"  Protein quality: {protein_info.get('quality_score', 'N/A')}/10")
                logger.info(f"  Chain count: {protein_info.get('chains', 'N/A')}")
                logger.info(f"  Residue count: {protein_info.get('residues', 'N/A')}")
            
            if summary.get('ligand_analysis'):
                ligand_info = summary['ligand_analysis']
                logger.info(f"  Molecular weight: {ligand_info.get('mol_weight', 'N/A')} Da")
                logger.info(f"  Drug-likeness score: {ligand_info.get('drug_likeness_score', 'N/A')}/10")
                logger.info(f"  LogP: {ligand_info.get('logp', 'N/A')}")
            
            logger.success("✅ Professional docking completed successfully!")
            return True
            
        else:
            error_msg = docking_results.get('error', 'Unknown error')
            logger.error(f"❌ Docking failed: {error_msg}")
            return False
            
    except Exception as e:
        logger.error(f"❌ Error in docking process: {str(e)}")
        import traceback
        logger.debug(f"Error details: {traceback.format_exc()}")
        return False

def add_simple_ligand_conversion_method():
    """Add simple ligand conversion method to preprocessing"""
    
    # Check if method exists in preprocessing
    preprocessing_file = "src/preprocessing.py"
    
    try:
        with open(preprocessing_file, 'r') as f:
            content = f.read()
        
        # If method doesn't exist, add it
        if '_convert_ligand_simple' not in content:
            simple_method = '''
    def _convert_ligand_simple(self, ligand_file: str, output_dir: str) -> str:
        """Simple ligand conversion to PDBQT without professional cleaning"""
        try:
            output_file = Path(output_dir) / f"{Path(ligand_file).stem}.pdbqt"
            
            # Use OpenBabel for conversion
            cmd = f"obabel {ligand_file} -O {output_file} -h"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and output_file.exists():
                logger.success(f"✅ Ligand converted: {output_file}")
                return str(output_file)
            else:
                logger.error("❌ Ligand conversion failed")
                return None
                
        except Exception as e:
            logger.error(f"❌ Error in ligand conversion: {str(e)}")
            return None
'''
            
            # Add method to class
            class_end = content.rfind('class PreprocessingEngine')
            if class_end != -1:
                next_class = content.find('\nclass ', class_end + 1)
                if next_class == -1:
                    next_class = len(content)
                
                new_content = content[:next_class] + simple_method + content[next_class:]
                
                with open(preprocessing_file, 'w') as f:
                    f.write(new_content)
                    
                logger.info("✅ Simple ligand conversion method added")
            
    except Exception as e:
        logger.warning(f"⚠️ Could not add simple conversion method: {str(e)}")

if __name__ == "__main__":
    # Add simple conversion method if needed
    add_simple_ligand_conversion_method()
    
    # Run docking
    success = run_professional_docking()
    
    if success:
        logger.info("\n✅ Success: Professional docking system worked successfully!")
        logger.info("Results saved in results/real_docking_results directory")
    else:
        logger.error("\n❌ Failed: There was a problem with the docking system") 