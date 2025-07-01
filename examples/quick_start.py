#!/usr/bin/env python3
"""
Quick Start Example for Orv-Mol
===============================

This script demonstrates how to use the Orv-Mol docking system with minimal setup.
"""

import os
import sys
from pathlib import Path

# Add src to path so we can import modules
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from docking import MolecularDocking
from preprocessing import PreprocessingEngine
from professional_cleaning import ProfessionalCleaner
from output_parser import OutputParser
from visualization import VisualizationEngine


def quick_docking_example():
    """Run a quick docking example with provided test data."""
    print("Orv-Mol Quick Start Example")
    print("=" * 40)
    
    # Define file paths
    project_root = Path(__file__).parent.parent
    protein_file = project_root / "data" / "input" / "protein.pdb"
    ligand_file = project_root / "data" / "input" / "ligand_small.pdb"
    output_dir = project_root / "example_output"
    
    # Check if input files exist
    if not protein_file.exists():
        print(f"❌ Protein file not found: {protein_file}")
        print("Please ensure test data is available in data/input/")
        return
    
    if not ligand_file.exists():
        print(f"❌ Ligand file not found: {ligand_file}")
        print("Please ensure test data is available in data/input/")
        return
    
    print(f"Protein: {protein_file.name}")
    print(f"Ligand: {ligand_file.name}")
    print(f"Output: {output_dir}")
    print()
    
    try:
        # Create output directory
        output_dir.mkdir(exist_ok=True)
        
        # Initialize engines
        print("Initializing engines...")
        cleaner = ProfessionalCleaner()
        preprocessor = PreprocessingEngine(professional_mode=True)
        docker = MolecularDocking()
        parser = OutputParser()
        
        # Step 1: Professional cleaning
        print("Step 1: Professional cleaning...")
        
        # Clean protein
        cleaned_protein = cleaner.clean_protein(str(protein_file))
        print(f"  ✅ Protein cleaned (quality: {cleaned_protein.get('quality_score', 'N/A'):.1f}/10)")
        
        # Clean ligand
        cleaned_ligand = cleaner.clean_ligand(str(ligand_file))
        print(f"  ✅ Ligand cleaned (drug-likeness: {cleaned_ligand.get('druglike_score', 'N/A'):.1f}/10)")
        
        # Step 2: Preprocessing
        print("Step 2: Preprocessing...")
        temp_dir = output_dir / "temp"
        temp_dir.mkdir(exist_ok=True)
        
        # Prepare files for docking
        ligand_pdbqt = preprocessor.prepare_ligand(str(ligand_file), str(temp_dir))
        protein_pdbqt = preprocessor.prepare_receptor(str(protein_file), str(temp_dir))
        
        # Calculate binding box (automatic)
        box_params = preprocessor.calculate_blind_docking_box(protein_pdbqt)
        print(f"  ✅ Binding box calculated: center=({box_params['center_x']:.1f}, {box_params['center_y']:.1f}, {box_params['center_z']:.1f})")
        
        # Step 3: Docking
        print("Step 3: Running molecular docking...")
        poses_dir = output_dir / "poses"
        poses_dir.mkdir(exist_ok=True)
        
        docking_results = docker.run_vina_docking(
            receptor=protein_pdbqt,
            ligand=ligand_pdbqt,
            box_params=box_params,
            output_dir=str(poses_dir),
            exhaustiveness=16,  # Good balance of speed vs accuracy
            num_modes=10
        )
        
        # Step 4: Parse results
        print("Step 4: Parsing results...")
        parsed_results = parser.parse_vina_output(
            docking_results['log_file'],
            docking_results['output_pdbqt'],
            str(poses_dir)
        )
        
        # Display results
        print("\nDOCKING RESULTS")
        print("=" * 30)
        print(f"Best binding affinity: {parsed_results['best_affinity']:.2f} kcal/mol")
        print(f"Number of poses: {len(parsed_results['poses'])}")
        
        print("\nTop 3 poses:")
        for i, pose in enumerate(parsed_results['poses'][:3], 1):
            print(f"  {i}. {pose['affinity']:.2f} kcal/mol")
        
        # Step 5: Optional visualization
        create_viz = input("\nCreate visualizations? (y/N): ").lower().strip()
        if create_viz == 'y':
            print("Step 5: Creating visualizations...")
            visualizer = VisualizationEngine()
            images_dir = output_dir / "images"
            images_dir.mkdir(exist_ok=True)
            
            visualizer.create_docking_images(
                protein_file=str(protein_file),
                ligand_pose=parsed_results['best_pose_pdb'],
                output_dir=str(images_dir),
                binding_site_center=(box_params['center_x'], box_params['center_y'], box_params['center_z'])
            )
            print("  ✅ Images created")
        
        print(f"\n✅ Example completed successfully!")
        print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        print("Please check the error message and try again.")


def professional_cleaning_example():
    """Example demonstrating professional cleaning features."""
    print("\nProfessional Cleaning Example")
    print("=" * 40)
    
    project_root = Path(__file__).parent.parent
    protein_file = project_root / "data" / "input" / "protein.pdb"
    ligand_file = project_root / "data" / "input" / "ligand_small.pdb"
    
    if not protein_file.exists() or not ligand_file.exists():
        print("❌ Test data not found")
        return
    
    try:
        # Initialize professional cleaner
        cleaner = ProfessionalCleaner()
        
        # Analyze protein
        print("Analyzing protein...")
        protein_analysis = cleaner.analyze_protein_structure(str(protein_file))
        
        print(f"  • Quality score: {protein_analysis.get('quality_score', 'N/A'):.1f}/10")
        print(f"  • Number of chains: {protein_analysis.get('chain_count', 'N/A')}")
        print(f"  • Number of residues: {protein_analysis.get('residue_count', 'N/A')}")
        print(f"  • Number of atoms: {protein_analysis.get('atom_count', 'N/A')}")
        
        # Analyze ligand
        print("\nAnalyzing ligand...")
        ligand_properties = cleaner.calculate_ligand_properties(str(ligand_file))
        
        print(f"  • Molecular weight: {ligand_properties.get('MW', 'N/A'):.1f} Da")
        print(f"  • LogP: {ligand_properties.get('LogP', 'N/A'):.2f}")
        print(f"  • Drug-likeness: {ligand_properties.get('druglike_score', 'N/A'):.1f}/10")
        print(f"  • Rotatable bonds: {ligand_properties.get('RotatableBonds', 'N/A')}")
        
        # Get optimized parameters
        print("\nRecommended docking parameters:")
        optimized_params = cleaner.optimize_docking_parameters(ligand_properties)
        
        print(f"  • Exhaustiveness: {optimized_params.get('exhaustiveness', 8)}")
        print(f"  • Number of modes: {optimized_params.get('num_modes', 9)}")
        print(f"  • Energy range: {optimized_params.get('energy_range', 3.0)} kcal/mol")
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")


if __name__ == "__main__":
    print("Welcome to Orv-Mol!")
    print("\nAvailable examples:")
    print("1. Quick docking example")
    print("2. Professional cleaning example")
    print("3. Both examples")
    
    choice = input("\nSelect option (1-3) or press Enter for quick example: ").strip()
    
    if choice == "2":
        professional_cleaning_example()
    elif choice == "3":
        professional_cleaning_example()
        print("\n" + "="*50)
        quick_docking_example()
    else:
        quick_docking_example() 