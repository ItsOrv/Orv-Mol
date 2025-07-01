#!/usr/bin/env python3
"""
Orv-Mol: Professional Molecular Docking System
==============================================

A comprehensive automated pipeline for molecular docking using AutoDock Vina
with professional-grade cleaning and optimization features.

Usage:
    python dock.py protein.pdb ligand.sdf [options]

Author: Orv Development Team
License: MIT
"""

import os
import sys
import argparse
import traceback
from pathlib import Path
from datetime import datetime

from loguru import logger
from tqdm import tqdm

# Import pipeline modules
try:
    from src.preprocessing import PreprocessingEngine
    from src.docking import DockingEngine
    from src.output_parser import OutputParser
    from src.visualization import VisualizationEngine
    from src.utils import setup_logging, validate_inputs, create_output_structure
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Please ensure you're running from the project root and all dependencies are installed.")
    sys.exit(1)


def create_parser():
    """Create comprehensive command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Orv-Mol: Professional Molecular Docking System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic docking
  python dock.py protein.pdb ligand.sdf
  
  # Professional mode with high accuracy
  python dock.py protein.pdb ligand.sdf --professional-cleaning --exhaustiveness 32
  
  # Targeted docking with known binding site
  python dock.py protein.pdb ligand.sdf --center-x 168.67 --center-y 29.82 --center-z 182.29
  
  # Quick mode for screening
  python dock.py protein.pdb ligand.sdf --quick --skip-visualization
  
  # Research-grade with all features
  python dock.py protein.pdb ligand.sdf --research-grade
        """
    )
    
    # Required arguments
    parser.add_argument("protein", help="Protein structure file (.pdb)")
    parser.add_argument("ligand", help="Ligand structure file (.mol2, .sdf, .pdb)")
    
    # Output Options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument("--output-dir", "-o", default="results", 
                            help="Output directory (default: results)")
    output_group.add_argument("--prefix", default="", 
                            help="File prefix for all outputs")
    output_group.add_argument("--timestamp", action="store_true",
                            help="Add timestamp to output directory name")
    output_group.add_argument("--overwrite", action="store_true",
                            help="Overwrite existing output directory")
    
    # Docking Parameters
    docking_group = parser.add_argument_group('Docking Parameters')
    docking_group.add_argument("--exhaustiveness", "-e", type=int, default=8,
                             help="AutoDock Vina exhaustiveness (default: 8, max: 32)")
    docking_group.add_argument("--num-modes", "-n", type=int, default=9,
                             help="Number of binding modes (default: 9)")
    docking_group.add_argument("--energy-range", type=float, default=3.0,
                             help="Energy range for binding modes (default: 3.0 kcal/mol)")
    docking_group.add_argument("--cpu", type=int, default=0,
                             help="Number of CPUs to use (0 = auto-detect)")
    docking_group.add_argument("--seed", type=int,
                             help="Random seed for reproducible results")
    
    # Binding Site Definition
    site_group = parser.add_argument_group('Binding Site Definition')
    site_group.add_argument("--center-x", type=float, help="Grid center X coordinate")
    site_group.add_argument("--center-y", type=float, help="Grid center Y coordinate")
    site_group.add_argument("--center-z", type=float, help="Grid center Z coordinate")
    site_group.add_argument("--size-x", type=float, default=22.5, 
                          help="Grid size X (default: 22.5 Å)")
    site_group.add_argument("--size-y", type=float, default=22.5, 
                          help="Grid size Y (default: 22.5 Å)")
    site_group.add_argument("--size-z", type=float, default=22.5, 
                          help="Grid size Z (default: 22.5 Å)")
    site_group.add_argument("--blind-docking", action="store_true", default=True,
                          help="Use blind docking (search entire protein surface)")
    site_group.add_argument("--binding-site-residues", 
                          help="Comma-separated list of binding site residues (e.g., 'A:123,A:124')")
    
    # Preprocessing Options
    prep_group = parser.add_argument_group('Preprocessing Options')
    prep_group.add_argument("--professional-cleaning", action="store_true", default=True,
                          help="Enable professional-grade cleaning (default: True)")
    prep_group.add_argument("--standard-mode", action="store_true",
                          help="Use standard preprocessing only (disable professional features)")
    prep_group.add_argument("--ph", type=float, default=7.4,
                          help="Target pH for protonation optimization (default: 7.4)")
    prep_group.add_argument("--auto-optimize-params", action="store_true", default=True,
                          help="Automatically optimize docking parameters")
    prep_group.add_argument("--add-hydrogens", action="store_true", default=True,
                          help="Add hydrogens to structures")
    prep_group.add_argument("--remove-waters", action="store_true", default=True,
                          help="Remove water molecules from protein")
    prep_group.add_argument("--remove-heteroatoms", action="store_true",
                          help="Remove heteroatoms from protein")
    prep_group.add_argument("--optimize-ligand", action="store_true", default=True,
                          help="Optimize ligand geometry")
    prep_group.add_argument("--generate-conformers", type=int, default=1,
                          help="Number of conformers to generate (default: 1)")
    
    # Analysis Options
    analysis_group = parser.add_argument_group('Analysis Options')
    analysis_group.add_argument("--drug-likeness", action="store_true", default=True,
                               help="Calculate drug-likeness properties")
    analysis_group.add_argument("--admet-properties", action="store_true",
                               help="Calculate ADMET properties")
    analysis_group.add_argument("--protein-quality", action="store_true", default=True,
                               help="Assess protein structure quality")
    analysis_group.add_argument("--interaction-analysis", action="store_true",
                               help="Analyze protein-ligand interactions")
    analysis_group.add_argument("--binding-site-analysis", action="store_true",
                               help="Perform binding site cavity analysis")
    
    # Visualization Options
    viz_group = parser.add_argument_group('Visualization Options')
    viz_group.add_argument("--visualize", action="store_true", default=True,
                         help="Generate visualization images")
    viz_group.add_argument("--animation", action="store_true",
                         help="Create 3D rotation animation")
    viz_group.add_argument("--interaction-plots", action="store_true",
                         help="Generate interaction diagrams")
    viz_group.add_argument("--surface-plots", action="store_true",
                         help="Generate surface representation plots")
    viz_group.add_argument("--image-format", choices=['png', 'jpg', 'svg', 'pdf'], 
                         default='png', help="Output image format")
    viz_group.add_argument("--image-dpi", type=int, default=300,
                         help="Image resolution in DPI (default: 300)")
    viz_group.add_argument("--image-size", default="1200x900",
                         help="Image size in pixels (default: 1200x900)")
    
    # Performance Presets
    perf_group = parser.add_argument_group('Performance Presets')
    perf_group.add_argument("--quick", action="store_true",
                          help="Quick mode: low exhaustiveness, no visualization")
    perf_group.add_argument("--fast", action="store_true",
                          help="Fast mode: moderate exhaustiveness")
    perf_group.add_argument("--accurate", action="store_true",
                          help="Accurate mode: high exhaustiveness")
    perf_group.add_argument("--research-grade", action="store_true",
                          help="Research-grade: maximum accuracy and all features")
    perf_group.add_argument("--batch", action="store_true",
                          help="Batch processing mode optimizations")
    perf_group.add_argument("--parallel-conformers", action="store_true",
                          help="Process multiple conformers in parallel")
    
    # Control Flags
    control_group = parser.add_argument_group('Control Flags')
    control_group.add_argument("--skip-preprocessing", action="store_true",
                             help="Skip preprocessing if .pdbqt files exist")
    control_group.add_argument("--skip-docking", action="store_true",
                             help="Skip docking step (preprocessing only)")
    control_group.add_argument("--skip-visualization", action="store_true",
                             help="Skip all visualization steps")
    control_group.add_argument("--skip-analysis", action="store_true",
                             help="Skip post-docking analysis")
    control_group.add_argument("--keep-temp", action="store_true",
                             help="Keep temporary files for debugging")
    control_group.add_argument("--dry-run", action="store_true",
                             help="Show what would be done without executing")
    
    # Reporting Options
    report_group = parser.add_argument_group('Reporting Options')
    report_group.add_argument("--reports", action="store_true", default=True,
                            help="Generate detailed reports")
    report_group.add_argument("--csv-output", action="store_true",
                            help="Export results to CSV format")
    report_group.add_argument("--json-output", action="store_true",
                            help="Export results to JSON format")
    report_group.add_argument("--xml-output", action="store_true",
                            help="Export results to XML format")
    report_group.add_argument("--summary-only", action="store_true",
                            help="Generate only summary report")
    
    # Advanced Options
    advanced_group = parser.add_argument_group('Advanced Options')
    advanced_group.add_argument("--config", help="Configuration file (.yaml/.json)")
    advanced_group.add_argument("--verbose", "-v", action="store_true",
                              help="Enable verbose output")
    advanced_group.add_argument("--quiet", "-q", action="store_true",
                              help="Suppress non-essential output")
    advanced_group.add_argument("--debug", action="store_true",
                              help="Enable debug mode")
    advanced_group.add_argument("--log-file", help="Save logs to file")
    advanced_group.add_argument("--profile", action="store_true",
                              help="Enable performance profiling")
    advanced_group.add_argument("--version", action="version", version="Orv-Mol v1.0.0")
    
    return parser


def apply_preset_modes(args):
    """Apply preset mode configurations."""
    if args.quick:
        args.exhaustiveness = 4
        args.num_modes = 5
        args.skip_visualization = True
        args.reports = False
        print("Quick mode: Fast execution with minimal output")
        
    elif args.fast:
        args.exhaustiveness = 8
        args.num_modes = 10
        args.animation = False
        print("Fast mode: Balanced speed and accuracy")
        
    elif args.accurate:
        args.exhaustiveness = 16
        args.num_modes = 20
        args.interaction_analysis = True
        print("Accurate mode: High precision docking")
        
    elif args.research_grade:
        args.exhaustiveness = 32
        args.num_modes = 50
        args.professional_cleaning = True
        args.admet_properties = True
        args.interaction_analysis = True
        args.animation = True
        args.csv_output = True
        print("Research-grade mode: Maximum accuracy and all features")
    
    if args.standard_mode:
        args.professional_cleaning = False
        print("Standard preprocessing mode enabled")


def validate_arguments(args):
    """Validate command line arguments."""
    errors = []
    
    # Check file existence
    if not Path(args.protein).exists():
        errors.append(f"Protein file not found: {args.protein}")
    if not Path(args.ligand).exists():
        errors.append(f"Ligand file not found: {args.ligand}")
    
    # Validate parameter ranges
    if args.exhaustiveness < 1 or args.exhaustiveness > 32:
        errors.append("Exhaustiveness must be between 1 and 32")
    if args.num_modes < 1 or args.num_modes > 100:
        errors.append("Number of modes must be between 1 and 100")
    if args.energy_range < 0.1 or args.energy_range > 10.0:
        errors.append("Energy range must be between 0.1 and 10.0 kcal/mol")
    if args.ph < 0 or args.ph > 14:
        errors.append("pH must be between 0 and 14")
    
    # Check conflicting options
    if args.quick and args.research_grade:
        errors.append("Cannot use both --quick and --research-grade modes")
    if args.standard_mode and args.professional_cleaning:
        errors.append("Cannot use both --standard-mode and --professional-cleaning")
    if args.quiet and args.verbose:
        errors.append("Cannot use both --quiet and --verbose modes")
    
    if errors:
        print("❌ Argument validation errors:")
        for error in errors:
            print(f"   • {error}")
        sys.exit(1)


def main():
    """Main pipeline execution function."""
    # Parse command line arguments
    parser = create_parser()
    args = parser.parse_args()
    
    # Show header
    print("Orv-Mol: Professional Molecular Docking System")
    print("=" * 55)
    
    # Apply preset modes
    apply_preset_modes(args)
    
    # Validate arguments
    validate_arguments(args)
    
    # Handle timestamp in output directory
    if args.timestamp:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.output_dir = f"{args.output_dir}_{timestamp}"
    
    # Setup logging
    if args.debug:
        log_level = "DEBUG"
    elif args.quiet:
        log_level = "WARNING"  
    elif args.verbose:
        log_level = "INFO"
    else:
        log_level = "INFO"
    
    setup_logging(log_level, args.log_file)
    
    # Handle dry run
    if args.dry_run:
        print("DRY RUN MODE - Configuration Preview:")
        print(f"  Input protein: {args.protein}")
        print(f"  Input ligand: {args.ligand}")
        print(f"  Output directory: {args.output_dir}")
        print(f"  Professional cleaning: {args.professional_cleaning}")
        print(f"  Exhaustiveness: {args.exhaustiveness}")
        print(f"  Number of modes: {args.num_modes}")
        print(f"  Visualization: {not args.skip_visualization}")
        print(f"  Analysis: {not args.skip_analysis}")
        if args.center_x is not None:
            print(f"  Binding site: ({args.center_x:.2f}, {args.center_y:.2f}, {args.center_z:.2f})")
        else:
            print(f"  Binding site: Automatic (blind docking)")
        return
    
    print(f"Protein: {Path(args.protein).name}")
    print(f"Ligand: {Path(args.ligand).name}")
    print(f"Output: {args.output_dir}")
    print()
    
    try:
        # Validate inputs
        logger.info("Validating input files...")
        validate_inputs(args.protein, args.ligand)
        
        # Create output directory structure
        logger.info("Creating output directory structure...")
        output_paths = create_output_structure(args.output_dir, overwrite=args.overwrite)
        
        # Initialize engines
        professional_mode = args.professional_cleaning and not args.standard_mode
        preprocessor = PreprocessingEngine(
            professional_mode=professional_mode,
            ph=args.ph
        )
        docker = DockingEngine()
        parser_engine = OutputParser()
        visualizer = VisualizationEngine() if not args.skip_visualization else None
        
        if professional_mode:
            logger.info(f"Professional cleaning enabled (pH={args.ph})")
        else:
            logger.info("Using standard preprocessing mode")
        
        # Print configuration summary
        logger.info("Configuration Summary:")
        logger.info(f"  • Exhaustiveness: {args.exhaustiveness}")
        logger.info(f"  • Binding modes: {args.num_modes}")
        logger.info(f"  • Energy range: {args.energy_range} kcal/mol")
        if args.seed:
            logger.info(f"  • Random seed: {args.seed}")
        
        # Step 1: Preprocessing
        if not args.skip_preprocessing:
            logger.info("Step 1: Preprocessing input files...")
            with tqdm(total=4, desc="Preprocessing") as pbar:
                # Convert and prepare ligand
                pbar.set_description("Converting ligand to PDBQT")
                ligand_pdbqt = preprocessor.prepare_ligand(
                    args.ligand, output_paths['temp']
                )
                pbar.update(1)
                
                # Prepare protein
                pbar.set_description("Preparing protein receptor")
                protein_pdbqt = preprocessor.prepare_receptor(
                    args.protein, output_paths['temp']
                )
                pbar.update(1)
                
                # Calculate binding box
                pbar.set_description("Calculating binding box")
                if args.center_x is None:  # Blind docking
                    box_params = preprocessor.calculate_blind_docking_box(protein_pdbqt)
                    logger.info("Using blind docking approach")
                else:  # Targeted docking
                    box_params = {
                        'center_x': args.center_x, 'center_y': args.center_y, 'center_z': args.center_z,
                        'size_x': args.size_x, 'size_y': args.size_y, 'size_z': args.size_z
                    }
                    logger.info("Using targeted docking approach")
                pbar.update(1)
                
                # Log box parameters
                logger.info(f"Grid box: center=({box_params['center_x']:.2f}, {box_params['center_y']:.2f}, {box_params['center_z']:.2f})")
                logger.info(f"Grid size: ({box_params['size_x']:.2f}, {box_params['size_y']:.2f}, {box_params['size_z']:.2f})")
                pbar.update(1)
        else:
            logger.info("Skipping preprocessing (files assumed to exist)")
            ligand_pdbqt = Path(output_paths['temp']) / f"{Path(args.ligand).stem}.pdbqt"
            protein_pdbqt = Path(output_paths['temp']) / f"{Path(args.protein).stem}.pdbqt"
            # Need to define box_params for skipped preprocessing
            box_params = {
                'center_x': args.center_x or 0, 'center_y': args.center_y or 0, 'center_z': args.center_z or 0,
                'size_x': args.size_x, 'size_y': args.size_y, 'size_z': args.size_z
            }
        
        # Skip docking if requested
        if args.skip_docking:
            logger.info("Skipping docking as requested")
            return
        
        # Step 2: Docking with optimized parameters
        logger.info("Step 2: Running AutoDock Vina...")
        
        # Get optimized parameters if auto-optimization is enabled
        if args.auto_optimize_params and professional_mode:
            logger.info("Optimizing docking parameters based on molecular analysis...")
            optimized_params = preprocessor.get_optimized_docking_parameters()
            
            # Use optimized parameters (with user overrides)
            if args.exhaustiveness == 8:  # default value
                args.exhaustiveness = optimized_params.get('exhaustiveness', 8)
            if args.num_modes == 9:  # default value
                args.num_modes = optimized_params.get('num_modes', 9)
            if args.energy_range == 3.0:  # default value
                args.energy_range = optimized_params.get('energy_range', 3.0)
            
            logger.info(f"Optimized parameters: exhaustiveness={args.exhaustiveness}, modes={args.num_modes}, range={args.energy_range}")
        
        docking_results = docker.run_vina_docking(
            receptor=protein_pdbqt,
            ligand=ligand_pdbqt,
            box_params=box_params,
            output_dir=output_paths['poses'],
            exhaustiveness=int(args.exhaustiveness),
            num_modes=int(args.num_modes),
            energy_range=float(args.energy_range)
        )
        
        # Step 3: Output parsing
        logger.info("Step 3: Parsing docking results...")
        parsed_results = parser_engine.parse_vina_output(
            docking_results['log_file'],
            docking_results['output_pdbqt'],
            output_paths['poses']
        )
        
        # Log results summary
        logger.success(f"✅ Docking completed successfully!")
        logger.info(f"Best binding affinity: {parsed_results['best_affinity']:.2f} kcal/mol")
        logger.info(f"Number of poses: {len(parsed_results['poses'])}")
        
        # Step 4: Analysis
        if not args.skip_analysis:
            if professional_mode:
                logger.info("Professional cleaning summary:")
                cleaning_summary = preprocessor.get_cleaning_summary()
                
                if cleaning_summary.get('protein_analysis'):
                    protein_analysis = cleaning_summary['protein_analysis']
                    logger.info(f"  Protein quality: {protein_analysis.get('quality_score', 'N/A'):.1f}/10")
                    logger.info(f"  Protein chains: {protein_analysis.get('chain_count', 'N/A')}")
                    
                if cleaning_summary.get('ligand_properties'):
                    ligand_props = cleaning_summary['ligand_properties']
                    logger.info(f"  Drug-likeness: {ligand_props.get('druglike_score', 'N/A'):.1f}/10")
                    logger.info(f"  Molecular weight: {ligand_props.get('MW', 'N/A'):.1f} Da")
                    logger.info(f"  Conformers tested: {ligand_props.get('conformer_count', 1)}")
        
        # Step 5: Visualization
        if not args.skip_visualization and visualizer:
            logger.info("Step 4: Generating visualizations...")
            
            # Create visualization images
            visualizer.create_docking_images(
                protein_file=args.protein,
                ligand_pose=parsed_results['best_pose_pdb'],
                output_dir=output_paths['images'],
                binding_site_center=(box_params['center_x'], box_params['center_y'], box_params['center_z']),
                image_format=args.image_format,
                dpi=args.image_dpi
            )
            
            # Generate interaction plots if requested
            if args.interaction_plots:
                visualizer.create_interaction_plots(
                    protein_file=args.protein,
                    ligand_pose=parsed_results['best_pose_pdb'],
                    output_dir=output_paths['images']
                )
            
            # Create animation if requested
            if args.animation:
                logger.info("Creating 3D rotation animation...")
                visualizer.create_rotation_animation(
                    protein_file=args.protein,
                    ligand_pose=parsed_results['best_pose_pdb'],
                    output_dir=output_paths['animation'],
                    center=(box_params['center_x'], box_params['center_y'], box_params['center_z'])
                )
        
        # Step 6: Generate reports
        if args.reports:
            logger.info("Step 5: Generating summary report...")
            report_data = generate_summary_report(args, parsed_results, box_params, output_paths['root'])
            
            # Export to various formats if requested
            if args.csv_output:
                export_to_csv(parsed_results, output_paths['root'])
            if args.json_output:
                export_to_json(parsed_results, output_paths['root'])
            if args.xml_output:
                export_to_xml(parsed_results, output_paths['root'])
        
        logger.success("✅ Pipeline completed successfully!")
        logger.info(f"All results saved to: {Path(args.output_dir).absolute()}")
        
        # Print summary
        print("\nDOCKING RESULTS SUMMARY")
        print("=" * 30)
        print(f"Best binding affinity: {parsed_results['best_affinity']:.2f} kcal/mol")
        print(f"Number of poses: {len(parsed_results['poses'])}")
        print(f"\nGenerated files:")
        for key, path in output_paths.items():
            if Path(path).exists():
                print(f"  {key.capitalize()}: {path}")
        
    except Exception as e:
        logger.error(f"❌ Pipeline failed: {str(e)}")
        if args.debug:
            logger.error(f"Full traceback: {traceback.format_exc()}")
        sys.exit(1)


def generate_summary_report(args, results, box_params, output_dir):
    """Generate comprehensive summary report."""
    report_file = Path(output_dir) / "docking_summary.txt"
    
    with open(report_file, 'w') as f:
        f.write("Orv-Mol Docking Summary Report\n")
        f.write("=" * 40 + "\n\n")
        
        f.write(f"Protein: {args.protein}\n")
        f.write(f"Ligand: {args.ligand}\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("Docking Parameters:\n")
        f.write(f"  Exhaustiveness: {args.exhaustiveness}\n")
        f.write(f"  Number of modes: {args.num_modes}\n")
        f.write(f"  Energy range: {args.energy_range} kcal/mol\n\n")
        
        f.write("Grid Box:\n")
        f.write(f"  Center: ({box_params['center_x']:.2f}, {box_params['center_y']:.2f}, {box_params['center_z']:.2f})\n")
        f.write(f"  Size: ({box_params['size_x']:.2f}, {box_params['size_y']:.2f}, {box_params['size_z']:.2f})\n\n")
        
        f.write("Results:\n")
        f.write(f"  Best binding affinity: {results['best_affinity']:.2f} kcal/mol\n")
        f.write(f"  Number of poses: {len(results['poses'])}\n\n")
        
        if results['poses']:
            f.write("Top binding poses:\n")
            for i, pose in enumerate(results['poses'][:5], 1):
                f.write(f"  {i}. {pose['affinity']:.2f} kcal/mol\n")
    
    logger.info(f"Summary report saved: {report_file}")
    return report_file


def export_to_csv(results, output_dir):
    """Export results to CSV format."""
    import csv
    csv_file = Path(output_dir) / "docking_results.csv"
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Pose', 'Binding Affinity (kcal/mol)', 'RMSD Lower', 'RMSD Upper'])
        
        for i, pose in enumerate(results['poses'], 1):
            writer.writerow([
                i, 
                pose['affinity'], 
                pose.get('rmsd_lb', 'N/A'),
                pose.get('rmsd_ub', 'N/A')
            ])
    
    logger.info(f"CSV export saved: {csv_file}")


def export_to_json(results, output_dir):
    """Export results to JSON format."""
    import json
    json_file = Path(output_dir) / "docking_results.json"
    
    export_data = {
        'timestamp': datetime.now().isoformat(),
        'best_affinity': results['best_affinity'],
        'num_poses': len(results['poses']),
        'poses': results['poses']
    }
    
    with open(json_file, 'w') as f:
        json.dump(export_data, f, indent=2)
    
    logger.info(f"JSON export saved: {json_file}")


def export_to_xml(results, output_dir):
    """Export results to XML format."""
    xml_file = Path(output_dir) / "docking_results.xml"
    
    with open(xml_file, 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<docking_results>\n')
        f.write(f'  <timestamp>{datetime.now().isoformat()}</timestamp>\n')
        f.write(f'  <best_affinity>{results["best_affinity"]}</best_affinity>\n')
        f.write(f'  <num_poses>{len(results["poses"])}</num_poses>\n')
        f.write('  <poses>\n')
        
        for i, pose in enumerate(results['poses'], 1):
            f.write(f'    <pose id="{i}">\n')
            f.write(f'      <affinity>{pose["affinity"]}</affinity>\n')
            f.write(f'      <rmsd_lb>{pose.get("rmsd_lb", "N/A")}</rmsd_lb>\n')
            f.write(f'      <rmsd_ub>{pose.get("rmsd_ub", "N/A")}</rmsd_ub>\n')
            f.write('    </pose>\n')
        
        f.write('  </poses>\n')
        f.write('</docking_results>\n')
    
    logger.info(f"XML export saved: {xml_file}")


def cleanup_temp_files(temp_dir):
    """Clean up temporary files if not requested to keep them."""
    import shutil
    
    temp_path = Path(temp_dir)
    if temp_path.exists():
        try:
            shutil.rmtree(temp_path)
            logger.info("Temporary files cleaned up")
        except Exception as e:
            logger.warning(f"⚠️ Could not clean up temp files: {e}")


if __name__ == "__main__":
    main() 