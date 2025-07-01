#!/usr/bin/env python3
"""
Example usage script for the Automated Molecular Docking Pipeline.
This script demonstrates various ways to use the pipeline programmatically.
"""

import os
import sys
from pathlib import Path

# Add src to path to import modules
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from preprocessing import PreprocessingEngine
from docking import DockingEngine
from output_parser import OutputParser
from visualization import VisualizationEngine
from utils import setup_logging, create_output_structure

def example_basic_docking():
    """Example: Basic docking workflow"""
    print("🧪 Example 1: Basic Docking Workflow")
    print("=" * 50)
    
    # Setup logging
    setup_logging("INFO")
    
    # Input files (replace with your actual files)
    protein_file = "examples/protein.pdb"
    ligand_file = "examples/ligand.sdf"
    output_dir = "example_results_basic"
    
    if not Path(protein_file).exists() or not Path(ligand_file).exists():
        print("❌ Example files not found. Please provide protein.pdb and ligand.sdf")
        return
    
    try:
        # Create output structure
        output_paths = create_output_structure(output_dir)
        
        # Initialize engines
        preprocessor = PreprocessingEngine()
        docker = DockingEngine()
        parser = OutputParser()
        
        # Step 1: Preprocess files
        print("📦 Preprocessing files...")
        ligand_pdbqt = preprocessor.prepare_ligand(ligand_file, output_paths['temp'])
        protein_pdbqt = preprocessor.prepare_receptor(protein_file, output_paths['temp'])
        
        # Step 2: Calculate binding box
        print("📏 Calculating binding box...")
        box_params = preprocessor.calculate_blind_docking_box(protein_pdbqt)
        
        # Step 3: Run docking
        print("⚙️ Running docking...")
        docking_results = docker.run_vina_docking(
            receptor=protein_pdbqt,
            ligand=ligand_pdbqt,
            box_params=box_params,
            output_dir=output_paths['poses'],
            exhaustiveness=8,
            num_modes=9
        )
        
        # Step 4: Parse results
        print("📊 Parsing results...")
        parsed_results = parser.parse_vina_output(
            docking_results['log_file'],
            docking_results['output_pdbqt'],
            output_paths['poses']
        )
        
        # Print summary
        print("\n🎯 Docking Results:")
        print(f"Best binding affinity: {parsed_results['best_affinity']:.2f} kcal/mol")
        print(f"Number of poses: {parsed_results['num_poses']}")
        
    except Exception as e:
        print(f"❌ Error: {e}")

def example_targeted_docking():
    """Example: Targeted docking with known binding site"""
    print("\n🧪 Example 2: Targeted Docking")
    print("=" * 50)
    
    # Known binding site coordinates (example)
    binding_site = {
        'center_x': 25.0,
        'center_y': 30.0,
        'center_z': 15.0,
        'size_x': 20.0,
        'size_y': 20.0,
        'size_z': 20.0
    }
    
    print(f"Using binding site: center=({binding_site['center_x']}, {binding_site['center_y']}, {binding_site['center_z']})")
    print("This approach is faster and more focused than blind docking.")

def example_high_performance_docking():
    """Example: High-performance docking settings"""
    print("\n🧪 Example 3: High-Performance Docking")
    print("=" * 50)
    
    # High-performance parameters
    hp_params = {
        'exhaustiveness': 32,  # Very thorough search
        'num_modes': 50,       # Many binding modes
        'energy_range': 5.0    # Larger energy range
    }
    
    print("High-performance parameters:")
    for key, value in hp_params.items():
        print(f"  {key}: {value}")
    
    print("⚠️  Warning: This will take significantly longer to run!")

def example_batch_processing():
    """Example: Batch processing multiple ligands"""
    print("\n🧪 Example 4: Batch Processing")
    print("=" * 50)
    
    # Simulate multiple ligand files
    ligand_files = [
        "ligand_1.sdf",
        "ligand_2.mol2", 
        "ligand_3.pdb"
    ]
    
    protein_file = "protein.pdb"
    
    print("Batch processing workflow:")
    print(f"Protein: {protein_file}")
    print(f"Ligands: {len(ligand_files)} files")
    
    for i, ligand in enumerate(ligand_files, 1):
        output_dir = f"batch_results_{i:03d}"
        print(f"  {i}. {ligand} -> {output_dir}/")
    
    print("\n# Example batch script:")
    print("for ligand in ligands/*.sdf; do")
    print("    name=$(basename $ligand .sdf)")
    print("    python dock.py protein.pdb $ligand -o results_$name")
    print("done")

def example_visualization_only():
    """Example: Visualization of existing results"""
    print("\n🧪 Example 5: Visualization Only")
    print("=" * 50)
    
    # Assume we have existing docking results
    protein_file = "protein.pdb"
    best_pose_file = "results/poses/best_pose.pdb"
    output_dir = "visualization_only"
    
    if Path(protein_file).exists() and Path(best_pose_file).exists():
        try:
            visualizer = VisualizationEngine()
            
            print("🖼️ Creating images...")
            visualizer.create_docking_images(
                protein_file=protein_file,
                ligand_pose=best_pose_file,
                output_dir=output_dir,
                binding_site_center=(25.0, 30.0, 15.0)
            )
            
            print("🎬 Creating animation...")
            visualizer.create_3d_animation(
                protein_file=protein_file,
                ligand_pose=best_pose_file,
                output_dir=output_dir,
                binding_site_center=(25.0, 30.0, 15.0)
            )
            
        except Exception as e:
            print(f"❌ Visualization error: {e}")
            print("Note: Requires PyMOL or ChimeraX")
    else:
        print("⚠️  Files not found. Run docking first or provide existing results.")

def example_custom_analysis():
    """Example: Custom analysis of docking results"""
    print("\n🧪 Example 6: Custom Analysis")
    print("=" * 50)
    
    # Example of custom post-processing
    results_file = "results/poses/docking_results.csv"
    
    if Path(results_file).exists():
        try:
            import pandas as pd
            
            # Load results
            df = pd.read_csv(results_file)
            
            print("📊 Custom Analysis Results:")
            print(f"Total poses: {len(df)}")
            print(f"Best affinity: {df['affinity'].min():.2f} kcal/mol")
            print(f"Worst affinity: {df['affinity'].max():.2f} kcal/mol")
            print(f"Mean affinity: {df['affinity'].mean():.2f} kcal/mol")
            print(f"Standard deviation: {df['affinity'].std():.2f} kcal/mol")
            
            # Count poses by affinity ranges
            strong = len(df[df['affinity'] <= -10.0])
            good = len(df[(df['affinity'] > -10.0) & (df['affinity'] <= -7.0)])
            weak = len(df[df['affinity'] > -7.0])
            
            print("\nBinding strength distribution:")
            print(f"  Strong (≤ -10.0 kcal/mol): {strong}")
            print(f"  Good (-10.0 to -7.0 kcal/mol): {good}")
            print(f"  Weak (> -7.0 kcal/mol): {weak}")
            
        except ImportError:
            print("❌ pandas not available for analysis")
        except Exception as e:
            print(f"❌ Analysis error: {e}")
    else:
        print("⚠️  Results file not found. Run docking first.")

def main():
    """Run all examples"""
    print("🧪 Automated Molecular Docking Pipeline - Examples")
    print("=" * 60)
    print()
    
    # Check if this is being run directly
    if len(sys.argv) > 1:
        example_name = sys.argv[1]
        
        examples = {
            'basic': example_basic_docking,
            'targeted': example_targeted_docking,
            'performance': example_high_performance_docking,
            'batch': example_batch_processing,
            'visualization': example_visualization_only,
            'analysis': example_custom_analysis
        }
        
        if example_name in examples:
            examples[example_name]()
        else:
            print(f"❌ Unknown example: {example_name}")
            print(f"Available examples: {', '.join(examples.keys())}")
    else:
        # Run all examples as demonstrations
        example_basic_docking()
        example_targeted_docking()
        example_high_performance_docking()
        example_batch_processing()
        example_visualization_only()
        example_custom_analysis()
        
        print("\n🎯 To run a specific example:")
        print("  python example_usage.py basic")
        print("  python example_usage.py targeted")
        print("  python example_usage.py performance")
        print("  python example_usage.py batch")
        print("  python example_usage.py visualization")
        print("  python example_usage.py analysis")

if __name__ == "__main__":
    main()
