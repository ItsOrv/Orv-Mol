#!/usr/bin/env python3
"""
Test Professional Docking with Small Ligand
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent / 'src'))

from loguru import logger
from preprocessing import PreprocessingEngine
from docking import DockingEngine
from output_parser import OutputParser
from utils import create_output_structure

def test_small_ligand_docking():
    """Test professional docking system with small ligand"""
    
    logger.info("Testing Professional Docking with Small Ligand")
    logger.info("=" * 50)
    
    # Test files
    project_root = Path(__file__).parent.parent
    protein_file = project_root / "data" / "input" / "protein.pdb"
    ligand_file = project_root / "data" / "input" / "ligand_small.pdb"
    output_dir = project_root / "results" / "test_results_small"
    
    # Verify input files exist
    if not protein_file.exists():
        logger.error(f"❌ Protein file not found: {protein_file}")
        return False
        
    if not ligand_file.exists():
        logger.error(f"❌ Ligand file not found: {ligand_file}")
        return False
    
    logger.info(f"Protein: {protein_file}")
    logger.info(f"Ligand: {ligand_file} (small molecule)")
    logger.info(f"Output: {output_dir}")
    
    try:
        # Create output structure
        logger.info("Creating output directory...")
        output_paths = create_output_structure(str(output_dir), overwrite=True)
        
        # Initialize professional preprocessing
        logger.info("Initializing Professional Preprocessing...")
        preprocessor = PreprocessingEngine(
            professional_mode=True,
            ph=7.4
        )
        logger.success("✅ Professional Preprocessing initialized")
        
        # Prepare protein with professional cleaning
        logger.info("Preparing protein with professional cleaning...")
        protein_pdbqt = preprocessor.prepare_receptor(
            str(protein_file), output_paths['temp']
        )
        logger.success(f"✅ Protein prepared: {protein_pdbqt}")
        
        # Prepare ligand with professional cleaning
        logger.info("Preparing ligand with professional cleaning...")
        ligand_pdbqt = preprocessor.prepare_ligand(
            str(ligand_file), output_paths['temp']
        )
        logger.success(f"✅ Ligand prepared: {ligand_pdbqt}")
        
        # Calculate binding box
        logger.info("Calculating binding box...")
        box_params = preprocessor.calculate_blind_docking_box(protein_pdbqt)
        
        logger.success("✅ Binding box calculated")
        logger.info(f"Box center: ({box_params['center_x']:.2f}, {box_params['center_y']:.2f}, {box_params['center_z']:.2f})")
        logger.info(f"Box size: ({box_params['size_x']:.2f}, {box_params['size_y']:.2f}, {box_params['size_z']:.2f})")
        
        # Get optimized parameters
        optimized_params = preprocessor.get_optimized_docking_parameters()
        logger.success("✅ Parameters optimized")
        
        # Run AutoDock Vina docking
        logger.info("Running AutoDock Vina docking...")
        docker = DockingEngine()
        
        docking_results = docker.run_vina_docking(
            receptor=protein_pdbqt,
            ligand=ligand_pdbqt,
            box_params=box_params,
            output_dir=output_paths['poses'],
            exhaustiveness=optimized_params.get('exhaustiveness', 8),
            num_modes=optimized_params.get('num_modes', 9),
            energy_range=optimized_params.get('energy_range', 3.0)
        )
        
        if docking_results and 'log_file' in docking_results:
            logger.success("✅ Docking completed successfully!")
            
            # Parse results
            parser = OutputParser()
            parsed_results = parser.parse_vina_output(
                docking_results['log_file'],
                docking_results['output_pdbqt'],
                output_paths['poses']
            )
            
            logger.info("Top 3 binding poses:")
            for i, pose in enumerate(parsed_results['poses'][:3], 1):
                logger.info(f"  Pose {i}: {pose['affinity']:.2f} kcal/mol")
        else:
            logger.error("❌ Docking failed")
            return False
        
        # Get cleaning summary
        logger.info("Generating cleaning summary...")
        cleaning_summary = preprocessor.get_cleaning_summary()
        
        if cleaning_summary.get('protein_analysis'):
            protein_info = cleaning_summary['protein_analysis']
            logger.info(f"Protein quality score: {protein_info.get('quality_score', 'N/A')}/10")
            logger.info(f"Protein chains: {protein_info.get('chain_count', 'N/A')}")
        
        if cleaning_summary.get('ligand_properties'):
            ligand_info = cleaning_summary['ligand_properties']
            logger.info(f"Ligand MW: {ligand_info.get('MW', 'N/A')} Da")
            logger.info(f"Drug-likeness: {ligand_info.get('druglike_score', 'N/A')}/10")
        
        logger.success("✅ Small ligand test completed successfully!")
        return True
        
    except Exception as e:
        logger.error(f"❌ Test failed: {str(e)}")
        import traceback
        logger.debug(traceback.format_exc())
        return False

if __name__ == "__main__":
    success = test_small_ligand_docking()
    
    if success:
        print("\n✅ Test PASSED: Small ligand docking works correctly")
    else:
        print("\n❌ Test FAILED: Check logs for details")
        sys.exit(1) 