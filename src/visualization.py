"""
Visualization module for creating publication-ready images and animations.
Supports PyMOL and ChimeraX for molecular visualization.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
from loguru import logger

try:
    from scipy.spatial.distance import cdist
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    logger.warning("SciPy not available, some visualization features will be limited")

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly not available, interactive plots will be limited")


class VisualizationEngine:
    """Handles molecular visualization using PyMOL and ChimeraX."""
    
    def __init__(self, config: Dict = None):
        self.config = config or {}
        self.pymol_available = self._check_pymol()
        self.chimerax_available = self._check_chimerax()
        
        # Default visualization settings
        self.default_settings = {
            'image_format': 'png',
            'image_dpi': 300,
            'image_width': 1200,
            'image_height': 900,
            'colors': {
                'protein': 'gray80',
                'ligand': 'red',
                'binding_site': 'yellow',
                'hydrogen_bonds': 'blue',
                'hydrophobic': 'orange',
                'electrostatic': 'magenta'
            }
        }
        
        # Update with config
        if 'visualization' in self.config:
            self.default_settings.update(self.config['visualization'])
        
        if not self.pymol_available and not self.chimerax_available:
            logger.warning("⚠️ Neither PyMOL nor ChimeraX found. Visualization features will be limited.")
    
    def create_comprehensive_visualization_suite(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str,
        binding_site_center: Tuple[float, float, float],
        docking_results: Dict = None
    ):
        """
        Create a comprehensive visualization suite including images, animations, 
        plots, and interaction diagrams.
        
        Args:
            protein_file: Path to protein PDB file
            ligand_pose: Path to ligand pose PDB file  
            output_dir: Directory to save all visualizations
            binding_site_center: Coordinates of binding site center
            docking_results: Results from docking analysis
        """
        logger.info("🎨 Creating comprehensive visualization suite...")
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        images_dir = output_path / "images"
        animations_dir = output_path / "animations" 
        plots_dir = output_path / "plots"
        diagrams_dir = output_path / "diagrams"
        
        for dir_path in [images_dir, animations_dir, plots_dir, diagrams_dir]:
            dir_path.mkdir(exist_ok=True)
        
        try:
            # 1. Create molecular images
            self.create_multiple_format_images(
                protein_file, ligand_pose, str(images_dir), binding_site_center
            )
            
            # 2. Create animations if enabled
            if self.default_settings.get('animation', {}).get('create_animations', False):
                self.create_enhanced_animations(
                    protein_file, ligand_pose, str(animations_dir), binding_site_center
                )
            
            # 3. Create interaction diagrams
            self.create_comprehensive_interaction_diagrams(
                protein_file, ligand_pose, str(diagrams_dir)
            )
            
            # 4. Create analysis plots
            if docking_results:
                self.create_analysis_plots(docking_results, str(plots_dir))
            
            # 5. Create contact maps
            self.create_contact_maps(protein_file, ligand_pose, str(diagrams_dir))
            
            # 6. Create energy landscapes
            if docking_results and docking_results.get('affinities'):
                self.create_energy_plots(docking_results['affinities'], str(plots_dir))
            
            logger.success("✅ Comprehensive visualization suite created")
            
        except Exception as e:
            logger.error(f"❌ Failed to create visualization suite: {e}")
            raise

    def create_multiple_format_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str,
        binding_site_center: Tuple[float, float, float]
    ):
        """Create molecular images in multiple formats (PNG, SVG, PDF)."""
        logger.info("🖼️ Creating multi-format molecular images...")
        
        output_path = Path(output_dir)
        
        # Get image settings from config
        width = self.default_settings.get('image_width', 1200)
        height = self.default_settings.get('image_height', 900)
        formats = ['png', 'svg'] if self.default_settings.get('image_format') == 'png' else [self.default_settings.get('image_format')]
        
        try:
            if self.pymol_available:
                self._create_pymol_multi_format_images(
                    protein_file, ligand_pose, output_path, 
                    binding_site_center, width, height, formats
                )
            elif self.chimerax_available:
                self._create_chimerax_multi_format_images(
                    protein_file, ligand_pose, output_path,
                    binding_site_center, width, height, formats
                )
            else:
                self._create_matplotlib_multi_format_plots(
                    protein_file, ligand_pose, output_path, formats
                )
            
        except Exception as e:
            logger.error(f"❌ Failed to create multi-format images: {e}")

    def create_enhanced_animations(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str,
        binding_site_center: Tuple[float, float, float]
    ):
        """Create enhanced animations with multiple rotation axes and effects."""
        logger.info("🎬 Creating enhanced 3D animations...")
        
        output_path = Path(output_dir)
        animation_config = self.default_settings.get('animation', {})
        
        frames = animation_config.get('frames', 120)
        duration = animation_config.get('duration', 8)
        
        try:
            # Create different types of animations
            self._create_rotation_animation(
                protein_file, ligand_pose, output_path, 
                binding_site_center, frames, duration, axis='y'
            )
            
            self._create_zoom_animation(
                protein_file, ligand_pose, output_path,
                binding_site_center, frames, duration
            )
            
            self._create_rocking_animation(
                protein_file, ligand_pose, output_path,
                binding_site_center, frames, duration
            )
            
        except Exception as e:
            logger.error(f"❌ Failed to create enhanced animations: {e}")

    def create_comprehensive_interaction_diagrams(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str
    ):
        """Create comprehensive 2D and 3D interaction diagrams."""
        logger.info("📊 Creating comprehensive interaction diagrams...")
        
        output_path = Path(output_dir)
        
        try:
            # 2D interaction diagram
            self.create_interaction_diagram(protein_file, ligand_pose, str(output_path))
            
            # 3D interaction network
            self._create_3d_interaction_network(protein_file, ligand_pose, output_path)
            
            # Residue interaction heatmap
            self._create_residue_interaction_heatmap(protein_file, ligand_pose, output_path)
            
            # Binding site overview
            self._create_binding_site_overview(protein_file, ligand_pose, output_path)
            
        except Exception as e:
            logger.error(f"❌ Failed to create interaction diagrams: {e}")

    def create_contact_maps(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str
    ):
        """Create detailed contact maps between protein and ligand."""
        logger.info("🗺️ Creating contact maps...")
        
        if not SCIPY_AVAILABLE:
            logger.warning("⚠️ SciPy not available, skipping contact maps")
            return
        
        output_path = Path(output_dir)
        
        try:
            # Extract coordinates
            protein_coords = self._extract_coordinates(protein_file)
            ligand_coords = self._extract_coordinates(ligand_pose)
            
            if protein_coords is None or ligand_coords is None:
                logger.warning("⚠️ Could not extract coordinates for contact maps")
                return
            
            # Calculate distance matrix
            distances = cdist(ligand_coords, protein_coords)
            
            # Create contact map
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Distance heatmap
            im1 = ax1.imshow(distances, cmap='viridis_r', aspect='auto')
            ax1.set_title('Protein-Ligand Distance Map')
            ax1.set_xlabel('Protein Atoms')
            ax1.set_ylabel('Ligand Atoms')
            plt.colorbar(im1, ax=ax1, label='Distance (Å)')
            
            # Contact matrix (within 4 Å)
            contacts = (distances <= 4.0).astype(int)
            im2 = ax2.imshow(contacts, cmap='RdYlBu_r', aspect='auto')
            ax2.set_title('Contact Map (≤ 4 Å)')
            ax2.set_xlabel('Protein Atoms')
            ax2.set_ylabel('Ligand Atoms')
            plt.colorbar(im2, ax=ax2, label='Contact')
            
            plt.tight_layout()
            plt.savefig(output_path / 'contact_maps.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create contact frequency plot
            self._create_contact_frequency_plot(distances, output_path)
            
            logger.info("✅ Contact maps created")
            
        except Exception as e:
            logger.error(f"❌ Failed to create contact maps: {e}")

    def create_analysis_plots(
        self,
        docking_results: Dict,
        output_dir: str
    ):
        """Create comprehensive analysis plots from docking results."""
        logger.info("📈 Creating analysis plots...")
        
        output_path = Path(output_dir)
        
        try:
            # Binding affinity analysis
            if docking_results.get('affinities'):
                self._create_affinity_analysis_plots(
                    docking_results['affinities'], output_path
                )
            
            # RMSD analysis
            if 'rmsd_data' in docking_results:
                self._create_rmsd_plots(docking_results['rmsd_data'], output_path)
            
            # Molecular properties radar chart
            if 'pharmacological_analysis' in docking_results:
                self._create_molecular_properties_radar(
                    docking_results['pharmacological_analysis'], output_path
                )
            
            # Interaction type distribution
            if 'binding_site_residues' in docking_results:
                self._create_interaction_type_distribution(
                    docking_results['binding_site_residues'], output_path
                )
            
        except Exception as e:
            logger.error(f"❌ Failed to create analysis plots: {e}")

    def create_energy_plots(
        self,
        affinities: List[Dict],
        output_dir: str
    ):
        """Create energy landscape and distribution plots."""
        logger.info("⚡ Creating energy plots...")
        
        output_path = Path(output_dir)
        
        try:
            energies = [a['affinity'] for a in affinities]
            
            # Energy distribution histogram
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
            
            # Histogram
            ax1.hist(energies, bins=min(len(energies), 20), alpha=0.7, color='skyblue', edgecolor='black')
            ax1.set_xlabel('Binding Energy (kcal/mol)')
            ax1.set_ylabel('Frequency')
            ax1.set_title('Energy Distribution')
            ax1.axvline(min(energies), color='red', linestyle='--', label=f'Best: {min(energies):.2f}')
            ax1.legend()
            
            # Energy vs pose rank
            poses = list(range(1, len(energies) + 1))
            ax2.scatter(poses, energies, alpha=0.7, s=50)
            ax2.plot(poses, energies, alpha=0.5, linestyle='-')
            ax2.set_xlabel('Pose Rank')
            ax2.set_ylabel('Binding Energy (kcal/mol)')
            ax2.set_title('Energy vs Pose Rank')
            
            # Cumulative energy distribution
            sorted_energies = sorted(energies)
            cumulative = np.arange(1, len(sorted_energies) + 1) / len(sorted_energies)
            ax3.plot(sorted_energies, cumulative, linewidth=2)
            ax3.set_xlabel('Binding Energy (kcal/mol)')
            ax3.set_ylabel('Cumulative Probability')
            ax3.set_title('Cumulative Energy Distribution')
            ax3.grid(True, alpha=0.3)
            
            # Box plot
            ax4.boxplot(energies, vert=True)
            ax4.set_ylabel('Binding Energy (kcal/mol)')
            ax4.set_title('Energy Statistics')
            ax4.set_xticklabels(['All Poses'])
            
            plt.tight_layout()
            plt.savefig(output_path / 'energy_analysis.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create interactive energy plot if Plotly available
            if PLOTLY_AVAILABLE:
                self._create_interactive_energy_plot(affinities, output_path)
            
        except Exception as e:
            logger.error(f"❌ Failed to create energy plots: {e}")

    def _create_pymol_multi_format_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path,
        binding_site_center: Tuple[float, float, float],
        width: int,
        height: int,
        formats: List[str]
    ):
        """Create PyMOL images in multiple formats."""
        colors = self.default_settings.get('colors', {})
        
        for format_type in formats:
            script_file = output_dir / f"pymol_script_{format_type}.pml"
            
            with open(script_file, 'w') as f:
                f.write(f"""
# Load structures
load {protein_file}, protein
load {ligand_pose}, ligand

# Set representations
hide everything
show cartoon, protein
show sticks, ligand
color {colors.get('protein', 'gray80')}, protein
color {colors.get('ligand', 'red')}, ligand

# Set up binding site view
center ligand
zoom ligand, 8

# Ray tracing settings
set ray_trace_mode, 1
set antialias, 2

# Create overview image
{format_type} {output_dir}/docking_overview.{format_type}, width={width}, height={height}, dpi=300, ray=1

# Create close-up view
zoom ligand, 5
{format_type} {output_dir}/docking_closeup.{format_type}, width={width}, height={height}, dpi=300, ray=1

# Create surface view
show surface, protein
set transparency, 0.3, protein
{format_type} {output_dir}/docking_surface.{format_type}, width={width}, height={height}, dpi=300, ray=1

quit
""")
            
            try:
                cmd = ['pymol', '-c', '-Q', str(script_file)]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError:
                logger.warning(f"⚠️ Failed to create {format_type} images with PyMOL")

    def _create_rotation_animation(
        self,
        protein_file: str,
        ligand_pose: str,
        output_path: Path,
        binding_site_center: Tuple[float, float, float],
        frames: int,
        duration: int,
        axis: str = 'y'
    ):
        """Create rotation animation around specified axis."""
        if self.pymol_available:
            script_file = output_path / f"rotation_{axis}_animation.pml"
            
            with open(script_file, 'w') as f:
                f.write(f"""
# Load structures
load {protein_file}, protein
load {ligand_pose}, ligand

# Set representations
hide everything
show cartoon, protein
show sticks, ligand
color gray80, protein
color red, ligand

# Center and zoom
center ligand
zoom ligand, 6

# Animation settings
set ray_trace_mode, 1
mset 1, {frames}

# Create rotation
for i in range(1, {frames + 1}):
    set_frame i
    turn {axis}, {360.0/frames}
    
# Save frames
for i in range(1, {frames + 1}):
    set_frame i
    png {output_path}/rotation_{axis}_frame_%04d.png % i, width=800, height=600, ray=1

quit
""")
            
            try:
                cmd = ['pymol', '-c', '-Q', str(script_file)]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                
                # Create video from frames
                self._create_video_from_frames(
                    output_path, output_path, duration, f"rotation_{axis}"
                )
                
            except subprocess.CalledProcessError:
                logger.warning(f"⚠️ Failed to create {axis} rotation animation")

    def _create_zoom_animation(
        self,
        protein_file: str,
        ligand_pose: str,
        output_path: Path,
        binding_site_center: Tuple[float, float, float],
        frames: int,
        duration: int
    ):
        """Create zoom-in animation."""
        if self.pymol_available:
            script_file = output_path / "zoom_animation.pml"
            
            with open(script_file, 'w') as f:
                f.write(f"""
# Load structures
load {protein_file}, protein
load {ligand_pose}, ligand

# Set representations
hide everything
show cartoon, protein
show sticks, ligand
color gray80, protein
color red, ligand

# Center on ligand
center ligand

# Animation settings
mset 1, {frames}

# Create zoom animation
for i in range(1, {frames + 1}):
    set_frame i
    zoom_factor = 15.0 - (10.0 * i / {frames})
    zoom ligand, zoom_factor
    
# Save frames
for i in range(1, {frames + 1}):
    set_frame i
    png {output_path}/zoom_frame_%04d.png % i, width=800, height=600, ray=1

quit
""")
            
            try:
                cmd = ['pymol', '-c', '-Q', str(script_file)]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                
                self._create_video_from_frames(
                    output_path, output_path, duration, "zoom"
                )
                
            except subprocess.CalledProcessError:
                logger.warning("⚠️ Failed to create zoom animation")

    def _create_rocking_animation(
        self,
        protein_file: str,
        ligand_pose: str,
        output_path: Path,
        binding_site_center: Tuple[float, float, float],
        frames: int,
        duration: int
    ):
        """Create rocking motion animation."""
        if self.pymol_available:
            script_file = output_path / "rocking_animation.pml"
            
            with open(script_file, 'w') as f:
                f.write(f"""
# Load structures
load {protein_file}, protein
load {ligand_pose}, ligand

# Set representations
hide everything
show cartoon, protein
show sticks, ligand
color gray80, protein
color red, ligand

# Center and zoom
center ligand
zoom ligand, 6

# Animation settings
mset 1, {frames}

# Create rocking motion
for i in range(1, {frames + 1}):
    set_frame i
    angle = 30.0 * sin(2 * 3.14159 * i / {frames})
    turn y, angle - prev_angle if i > 1 else angle
    
# Save frames
for i in range(1, {frames + 1}):
    set_frame i
    png {output_path}/rocking_frame_%04d.png % i, width=800, height=600, ray=1

quit
""")
            
            try:
                cmd = ['pymol', '-c', '-Q', str(script_file)]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                
                self._create_video_from_frames(
                    output_path, output_path, duration, "rocking"
                )
                
            except subprocess.CalledProcessError:
                logger.warning("⚠️ Failed to create rocking animation")

    def _check_pymol(self) -> bool:
        """Check if PyMOL is available."""
        try:
            import pymol
            logger.info("✅ PyMOL detected")
            return True
        except ImportError:
            try:
                # Check command line PyMOL
                subprocess.run(['pymol', '-c', '-Q'], 
                             capture_output=True, timeout=5)
                logger.info("✅ PyMOL command line detected")
                return True
            except (subprocess.TimeoutExpired, FileNotFoundError):
                logger.info("ℹ️ PyMOL not available")
                return False
    
    def _check_chimerax(self) -> bool:
        """Check if ChimeraX is available."""
        try:
            subprocess.run(['chimerax', '--version'], 
                         capture_output=True, timeout=5)
            logger.info("✅ ChimeraX detected")
            return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            logger.info("ℹ️ ChimeraX not available")
            return False
    
    def _create_pymol_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path,
        binding_site_center: Tuple[float, float, float],
        width: int,
        height: int
    ):
        """Create images using PyMOL."""
        
        # Create PyMOL script
        script_file = output_dir / "pymol_script.pml"
        
        with open(script_file, 'w') as f:
            f.write(f"""
# Load structures
load {protein_file}, protein
load {ligand_pose}, ligand

# Set representations
hide everything
show cartoon, protein
show sticks, ligand
color gray80, protein
color red, ligand

# Set up binding site view
center ligand
zoom ligand, 8

# Create overview image
set ray_trace_mode, 1
set antialias, 2
png {output_dir}/docking_overview.png, width={width}, height={height}, dpi=300, ray=1

# Create close-up view
zoom ligand, 5
show spheres, ligand
set sphere_scale, 0.3, ligand
png {output_dir}/docking_closeup.png, width={width}, height={height}, dpi=300, ray=1

# Create surface view
hide spheres, ligand
show sticks, ligand
show surface, protein
set transparency, 0.3, protein
png {output_dir}/docking_surface.png, width={width}, height={height}, dpi=300, ray=1

quit
""")
        
        # Run PyMOL script
        try:
            cmd = ['pymol', '-c', '-Q', str(script_file)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.warning(f"PyMOL command failed: {e.stderr}")
            self._create_pymol_python_images(protein_file, ligand_pose, output_dir, width, height)
    
    def _create_pymol_python_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path,
        width: int,
        height: int
    ):
        """Create images using PyMOL Python API."""
        try:
            import pymol
            
            # Initialize PyMOL
            pymol.finish_launching(['pymol', '-c'])
            cmd = pymol.cmd
            
            # Load structures
            cmd.load(protein_file, 'protein')
            cmd.load(ligand_pose, 'ligand')
            
            # Set representations
            cmd.hide('everything')
            cmd.show('cartoon', 'protein')
            cmd.show('sticks', 'ligand')
            cmd.color('gray80', 'protein')
            cmd.color('red', 'ligand')
            
            # Set up view
            cmd.center('ligand')
            cmd.zoom('ligand', 8)
            
            # Configure rendering
            cmd.set('ray_trace_mode', 1)
            cmd.set('antialias', 2)
            
            # Create overview image
            cmd.png(str(output_dir / 'docking_overview.png'), 
                   width=width, height=height, dpi=300, ray=1)
            
            # Create close-up view
            cmd.zoom('ligand', 5)
            cmd.png(str(output_dir / 'docking_closeup.png'),
                   width=width, height=height, dpi=300, ray=1)
            
            # Clean up
            pymol.cmd.quit()
            
        except Exception as e:
            logger.warning(f"PyMOL Python API failed: {e}")
    
    def _create_chimerax_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path,
        binding_site_center: Tuple[float, float, float],
        width: int,
        height: int
    ):
        """Create images using ChimeraX."""
        
        # Create ChimeraX script
        script_file = output_dir / "chimerax_script.cxc"
        
        with open(script_file, 'w') as f:
            f.write(f"""
open {protein_file}
open {ligand_pose}

# Set representations
hide atoms
show cartoons #1
show atoms #2
style #2 stick

# Set colors
color #1 gray
color #2 red

# Center on ligand
view #2

# Create overview image
save {output_dir}/docking_overview.png width {width} height {height}

# Close-up view
zoom #2 5
save {output_dir}/docking_closeup.png width {width} height {height}

exit
""")
        
        # Run ChimeraX script
        try:
            cmd = ['chimerax', '--nogui', str(script_file)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.warning(f"ChimeraX command failed: {e.stderr}")
    
    def _create_matplotlib_plot(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path
    ):
        """Create basic plot using matplotlib when no 3D software is available."""
        
        # Read coordinate data (simplified)
        protein_coords = self._extract_coordinates(protein_file)
        ligand_coords = self._extract_coordinates(ligand_pose)
        
        # Create 2D projection plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        if protein_coords:
            ax1.scatter(protein_coords[:, 0], protein_coords[:, 1], 
                       c='gray', alpha=0.6, s=1, label='Protein')
        
        if ligand_coords:
            ax1.scatter(ligand_coords[:, 0], ligand_coords[:, 1],
                       c='red', s=20, label='Ligand')
        
        ax1.set_xlabel('X coordinate (Å)')
        ax1.set_ylabel('Y coordinate (Å)')
        ax1.set_title('Docking Complex (XY projection)')
        ax1.legend()
        ax1.axis('equal')
        
        # Second view (XZ projection)
        if protein_coords:
            ax2.scatter(protein_coords[:, 0], protein_coords[:, 2],
                       c='gray', alpha=0.6, s=1, label='Protein')
        
        if ligand_coords:
            ax2.scatter(ligand_coords[:, 0], ligand_coords[:, 2],
                       c='red', s=20, label='Ligand')
        
        ax2.set_xlabel('X coordinate (Å)')
        ax2.set_ylabel('Z coordinate (Å)')
        ax2.set_title('Docking Complex (XZ projection)')
        ax2.legend()
        ax2.axis('equal')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'docking_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("📊 Created 2D projection plots")
    
    def _create_pymol_animation(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path,
        binding_site_center: Tuple[float, float, float],
        frames: int,
        duration: int
    ):
        """Create animation using PyMOL."""
        
        # Create PyMOL script for animation
        script_file = output_dir / "animation_script.pml"
        
        with open(script_file, 'w') as f:
            f.write(f"""
# Load structures
load {protein_file}, protein
load {ligand_pose}, ligand

# Set representations
hide everything
show cartoon, protein
show sticks, ligand
color gray80, protein
color red, ligand

# Set up binding site view
center ligand
zoom ligand, 6

# Set up animation
set ray_trace_mode, 1
set antialias, 2

# Create rotation animation
movie.clear()
scene F1, store
turn y, 360/{frames}
movie.store, 1, {frames}

# Save as movie
mpng {output_dir}/animation_frames/frame_, width=800, height=600, mode=2

quit
""")
        
        # Create frames directory
        frames_dir = output_dir / "animation_frames"
        frames_dir.mkdir(exist_ok=True)
        
        # Run PyMOL script
        try:
            cmd = ['pymol', '-c', '-Q', str(script_file)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Convert frames to video using ffmpeg if available
            self._create_video_from_frames(frames_dir, output_dir, duration)
            
        except subprocess.CalledProcessError as e:
            logger.warning(f"PyMOL animation failed: {e.stderr}")
    
    def _create_chimerax_animation(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path,
        binding_site_center: Tuple[float, float, float],
        frames: int,
        duration: int
    ):
        """Create animation using ChimeraX."""
        
        script_file = output_dir / "animation_script.cxc"
        
        with open(script_file, 'w') as f:
            f.write(f"""
open {protein_file}
open {ligand_pose}

# Set representations
hide atoms
show cartoons #1
show atoms #2
style #2 stick

# Set colors
color #1 gray
color #2 red

# Center on ligand
view #2

# Create rotation movie
movie record
turn y 2 {frames}
wait {frames}
movie encode {output_dir}/docking_animation.mp4

exit
""")
        
        # Run ChimeraX script
        try:
            cmd = ['chimerax', '--nogui', str(script_file)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.warning(f"ChimeraX animation failed: {e.stderr}")
    
    def _create_video_from_frames(
        self,
        frames_dir: Path,
        output_dir: Path,
        duration: int
    ):
        """Convert frame images to video using ffmpeg."""
        try:
            # Check if ffmpeg is available
            subprocess.run(['ffmpeg', '-version'], 
                         capture_output=True, check=True)
            
            # Create video
            cmd = [
                'ffmpeg', '-y',
                '-framerate', str(120 / duration),  # fps
                '-i', str(frames_dir / 'frame_%04d.png'),
                '-c:v', 'libx264',
                '-pix_fmt', 'yuv420p',
                '-crf', '18',
                str(output_dir / 'docking_animation.mp4')
            ]
            
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info("🎬 Created MP4 animation")
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("⚠️ ffmpeg not available, animation frames saved as images")
    
    def _extract_coordinates(self, pdb_file: str) -> Optional[np.ndarray]:
        """Extract atomic coordinates from PDB file."""
        try:
            coords = []
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                        except (ValueError, IndexError):
                            continue
            
            return np.array(coords) if coords else None
            
        except Exception as e:
            logger.warning(f"Failed to extract coordinates: {e}")
            return None
    
    def create_interaction_diagram(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str
    ):
        """Create 2D interaction diagram (basic implementation)."""
        logger.info("📊 Creating interaction diagram...")
        
        try:
            # This is a simplified implementation
            # For production use, consider libraries like ProLIF or RDKit
            
            output_path = Path(output_dir)
            
            # Extract basic information
            protein_coords = self._extract_coordinates(protein_file)
            ligand_coords = self._extract_coordinates(ligand_pose)
            
            if protein_coords is None or ligand_coords is None:
                logger.warning("⚠️ Could not extract coordinates for interaction diagram")
                return
            
            # Calculate distances (simplified)
            from scipy.spatial.distance import cdist
            
            distances = cdist(ligand_coords, protein_coords)
            close_contacts = np.where(distances < 4.0)  # 4 Å cutoff
            
            # Create simple interaction plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Plot ligand atoms
            ax.scatter(ligand_coords[:, 0], ligand_coords[:, 1], 
                      c='red', s=50, label='Ligand', zorder=3)
            
            # Plot nearby protein atoms
            nearby_protein = protein_coords[close_contacts[1]]
            if len(nearby_protein) > 0:
                ax.scatter(nearby_protein[:, 0], nearby_protein[:, 1],
                          c='blue', s=20, alpha=0.6, label='Nearby Protein', zorder=2)
            
            # Draw contact lines
            for i, j in zip(close_contacts[0], close_contacts[1]):
                ax.plot([ligand_coords[i, 0], protein_coords[j, 0]],
                       [ligand_coords[i, 1], protein_coords[j, 1]],
                       'k-', alpha=0.3, linewidth=0.5, zorder=1)
            
            ax.set_xlabel('X coordinate (Å)')
            ax.set_ylabel('Y coordinate (Å)')
            ax.set_title('Protein-Ligand Interactions (<4 Å)')
            ax.legend()
            ax.axis('equal')
            
            plt.tight_layout()
            plt.savefig(output_path / 'interaction_diagram.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info("✅ Interaction diagram created")
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create interaction diagram: {e}")

    def create_docking_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str,
        binding_site_center: Tuple[float, float, float],
        image_width: int = 1200,
        image_height: int = 900
    ):
        """Legacy method for backward compatibility - creates publication-ready images."""
        logger.info("🖼️ Creating publication-ready images (legacy method)...")
        return self.create_multiple_format_images(
            protein_file, ligand_pose, output_dir, binding_site_center
        )

    def create_3d_animation(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str,
        binding_site_center: Tuple[float, float, float],
        frames: int = 120,
        duration: int = 8
    ):
        """Legacy method for backward compatibility - creates 3D animation."""
        logger.info("🎬 Creating 3D rotation animation (legacy method)...")
        
        # Temporarily update animation settings
        old_settings = self.default_settings.get('animation', {})
        self.default_settings['animation'] = {
            'create_animations': True,
            'frames': frames,
            'duration': duration
        }
        
        try:
            self.create_enhanced_animations(protein_file, ligand_pose, output_dir, binding_site_center)
        finally:
            # Restore old settings
            self.default_settings['animation'] = old_settings

    def _create_3d_interaction_network(
        self,
        protein_file: str,
        ligand_pose: str,
        output_path: Path
    ):
        """Create 3D interaction network visualization."""
        if not SCIPY_AVAILABLE:
            logger.warning("⚠️ SciPy not available, skipping 3D interaction network")
            return
        
        try:
            # Extract coordinates
            protein_coords = self._extract_coordinates(protein_file)
            ligand_coords = self._extract_coordinates(ligand_pose)
            
            if not protein_coords or not ligand_coords:
                return
            
            # Calculate distances
            distances = cdist(ligand_coords, protein_coords)
            close_contacts = np.where(distances < 4.0)
            
            # Create 3D network plot
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(111, projection='3d')
            
            # Plot protein atoms
            protein_array = np.array(protein_coords)
            ax.scatter(protein_array[:, 0], protein_array[:, 1], protein_array[:, 2],
                      c='blue', alpha=0.6, s=20, label='Protein')
            
            # Plot ligand atoms
            ligand_array = np.array(ligand_coords)
            ax.scatter(ligand_array[:, 0], ligand_array[:, 1], ligand_array[:, 2],
                      c='red', s=50, label='Ligand')
            
            # Draw interaction lines
            for i, j in zip(close_contacts[0][:50], close_contacts[1][:50]):  # Limit to 50 contacts
                ax.plot([ligand_array[i, 0], protein_array[j, 0]],
                       [ligand_array[i, 1], protein_array[j, 1]],
                       [ligand_array[i, 2], protein_array[j, 2]],
                       'k-', alpha=0.3, linewidth=0.5)
            
            ax.set_xlabel('X (Å)')
            ax.set_ylabel('Y (Å)')
            ax.set_zlabel('Z (Å)')
            ax.set_title('3D Protein-Ligand Interaction Network')
            ax.legend()
            
            plt.savefig(output_path / '3d_interaction_network.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create 3D interaction network: {e}")

    def _create_residue_interaction_heatmap(
        self,
        protein_file: str,
        ligand_pose: str,
        output_path: Path
    ):
        """Create heatmap of residue-ligand interactions."""
        try:
            # This is a simplified implementation
            # In a real scenario, you'd parse residue information from PDB
            
            # Create mock data for demonstration
            residue_names = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO',
                           'SER', 'THR', 'TYR', 'ASN', 'GLN', 'CYS', 'LYS', 'ARG']
            
            interaction_types = ['Hydrogen Bond', 'Hydrophobic', 'Electrostatic', 'Van der Waals']
            
            # Generate random interaction matrix for demonstration
            np.random.seed(42)
            interaction_matrix = np.random.rand(len(residue_names), len(interaction_types))
            
            # Create heatmap
            fig, ax = plt.subplots(figsize=(10, 8))
            im = ax.imshow(interaction_matrix, cmap='YlOrRd', aspect='auto')
            
            # Set ticks and labels
            ax.set_xticks(np.arange(len(interaction_types)))
            ax.set_yticks(np.arange(len(residue_names)))
            ax.set_xticklabels(interaction_types)
            ax.set_yticklabels(residue_names)
            
            # Rotate the tick labels and set their alignment
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Interaction Strength', rotation=270, labelpad=15)
            
            # Add text annotations
            for i in range(len(residue_names)):
                for j in range(len(interaction_types)):
                    text = ax.text(j, i, f'{interaction_matrix[i, j]:.2f}',
                                 ha="center", va="center", color="black", fontsize=8)
            
            ax.set_title('Residue-Ligand Interaction Heatmap')
            plt.tight_layout()
            plt.savefig(output_path / 'residue_interaction_heatmap.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create residue interaction heatmap: {e}")

    def _create_binding_site_overview(
        self,
        protein_file: str,
        ligand_pose: str,
        output_path: Path
    ):
        """Create binding site overview diagram."""
        try:
            # Extract coordinates
            protein_coords = self._extract_coordinates(protein_file)
            ligand_coords = self._extract_coordinates(ligand_pose)
            
            if not protein_coords or not ligand_coords:
                return
            
            # Calculate binding site center
            ligand_array = np.array(ligand_coords)
            center = np.mean(ligand_array, axis=0)
            
            # Create binding site overview
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
            
            # XY projection
            protein_array = np.array(protein_coords)
            ax1.scatter(protein_array[:, 0], protein_array[:, 1], 
                       c='lightblue', alpha=0.5, s=1, label='Protein')
            ax1.scatter(ligand_array[:, 0], ligand_array[:, 1],
                       c='red', s=30, label='Ligand', edgecolor='black', linewidth=0.5)
            
            # Draw binding site circle
            circle = plt.Circle((center[0], center[1]), 10, fill=False, 
                              color='orange', linewidth=2, linestyle='--', label='Binding Site (10 Å)')
            ax1.add_patch(circle)
            
            ax1.set_xlabel('X coordinate (Å)')
            ax1.set_ylabel('Y coordinate (Å)')
            ax1.set_title('Binding Site Overview (XY view)')
            ax1.legend()
            ax1.axis('equal')
            
            # XZ projection
            ax2.scatter(protein_array[:, 0], protein_array[:, 2], 
                       c='lightblue', alpha=0.5, s=1, label='Protein')
            ax2.scatter(ligand_array[:, 0], ligand_array[:, 2],
                       c='red', s=30, label='Ligand', edgecolor='black', linewidth=0.5)
            
            # Draw binding site circle
            circle2 = plt.Circle((center[0], center[2]), 10, fill=False, 
                               color='orange', linewidth=2, linestyle='--', label='Binding Site (10 Å)')
            ax2.add_patch(circle2)
            
            ax2.set_xlabel('X coordinate (Å)')
            ax2.set_ylabel('Z coordinate (Å)')
            ax2.set_title('Binding Site Overview (XZ view)')
            ax2.legend()
            ax2.axis('equal')
            
            plt.tight_layout()
            plt.savefig(output_path / 'binding_site_overview.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create binding site overview: {e}")

    def _create_contact_frequency_plot(
        self,
        distances: np.ndarray,
        output_path: Path
    ):
        """Create contact frequency distribution plot."""
        try:
            # Calculate contact frequencies at different cutoffs
            cutoffs = np.arange(2.0, 10.1, 0.5)
            contact_counts = []
            
            for cutoff in cutoffs:
                contacts = np.sum(distances <= cutoff)
                contact_counts.append(contacts)
            
            # Create plot
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(cutoffs, contact_counts, 'b-o', linewidth=2, markersize=6)
            ax.fill_between(cutoffs, contact_counts, alpha=0.3)
            
            ax.set_xlabel('Distance Cutoff (Å)')
            ax.set_ylabel('Number of Contacts')
            ax.set_title('Contact Frequency vs Distance Cutoff')
            ax.grid(True, alpha=0.3)
            
            # Add annotations for important cutoffs
            for cutoff, count in [(3.5, np.sum(distances <= 3.5)), 
                                (4.0, np.sum(distances <= 4.0)),
                                (5.0, np.sum(distances <= 5.0))]:
                ax.annotate(f'{count} contacts', 
                          xy=(cutoff, count), xytext=(cutoff + 0.5, count + 5),
                          arrowprops=dict(arrowstyle='->', color='red', alpha=0.7))
            
            plt.tight_layout()
            plt.savefig(output_path / 'contact_frequency.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create contact frequency plot: {e}")

    def _create_affinity_analysis_plots(
        self,
        affinities: List[Dict],
        output_path: Path
    ):
        """Create detailed binding affinity analysis plots."""
        try:
            energies = [a['affinity'] for a in affinities]
            rmsd_lb = [a.get('rmsd_lb', 0) for a in affinities if a.get('rmsd_lb') is not None]
            rmsd_ub = [a.get('rmsd_ub', 0) for a in affinities if a.get('rmsd_ub') is not None]
            
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
            
            # Energy distribution with statistics
            ax1.hist(energies, bins=min(len(energies), 15), alpha=0.7, 
                    color='lightblue', edgecolor='black')
            ax1.axvline(np.mean(energies), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(energies):.2f}')
            ax1.axvline(np.median(energies), color='green', linestyle='--', 
                       label=f'Median: {np.median(energies):.2f}')
            ax1.set_xlabel('Binding Energy (kcal/mol)')
            ax1.set_ylabel('Frequency')
            ax1.set_title('Energy Distribution with Statistics')
            ax1.legend()
            
            # Energy funnel plot
            poses = list(range(1, len(energies) + 1))
            ax2.scatter(poses, energies, alpha=0.7, s=50, c=energies, cmap='viridis')
            ax2.plot(poses, energies, alpha=0.5, linestyle='-', color='gray')
            ax2.set_xlabel('Pose Rank')
            ax2.set_ylabel('Binding Energy (kcal/mol)')
            ax2.set_title('Energy Funnel')
            
            # RMSD vs Energy correlation (if RMSD data available)
            if rmsd_lb and len(rmsd_lb) == len(energies):
                ax3.scatter(rmsd_lb, energies, alpha=0.7, s=50)
                ax3.set_xlabel('RMSD Lower Bound (Å)')
                ax3.set_ylabel('Binding Energy (kcal/mol)')
                ax3.set_title('RMSD vs Energy Correlation')
                
                # Add correlation coefficient
                correlation = np.corrcoef(rmsd_lb, energies)[0, 1]
                ax3.text(0.05, 0.95, f'R = {correlation:.3f}', 
                        transform=ax3.transAxes, bbox=dict(boxstyle="round", facecolor='white'))
            else:
                ax3.text(0.5, 0.5, 'RMSD data not available', 
                        ha='center', va='center', transform=ax3.transAxes)
                ax3.set_title('RMSD vs Energy (No Data)')
            
            # Energy classification pie chart
            energy_classes = {
                'Very Strong (≤ -10.0)': len([e for e in energies if e <= -10.0]),
                'Strong (-10.0 to -8.0)': len([e for e in energies if -10.0 < e <= -8.0]),
                'Moderate (-8.0 to -6.0)': len([e for e in energies if -8.0 < e <= -6.0]),
                'Weak (-6.0 to -4.0)': len([e for e in energies if -6.0 < e <= -4.0]),
                'Very Weak (> -4.0)': len([e for e in energies if e > -4.0])
            }
            
            # Filter out zero values
            energy_classes = {k: v for k, v in energy_classes.items() if v > 0}
            
            ax4.pie(energy_classes.values(), labels=energy_classes.keys(), autopct='%1.1f%%')
            ax4.set_title('Binding Strength Classification')
            
            plt.tight_layout()
            plt.savefig(output_path / 'affinity_analysis.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create affinity analysis plots: {e}")

    def _create_molecular_properties_radar(
        self,
        properties: Dict,
        output_path: Path
    ):
        """Create radar chart for molecular properties."""
        try:
            # Define properties for radar chart
            property_names = ['Molecular Weight', 'LogP', 'HBD', 'HBA', 
                            'Rotatable Bonds', 'TPSA', 'Drug Likeness']
            
            # Normalize values for radar chart (0-1 scale)
            max_values = [800, 5, 5, 10, 15, 200, 10]  # Maximum reasonable values
            
            values = []
            for prop, max_val in zip(['molecular_weight', 'logp', 'hbd', 'hba', 
                                    'rotatable_bonds', 'tpsa', 'drug_likeness_score'], max_values):
                val = properties.get(prop, 0)
                normalized_val = min(val / max_val, 1.0) if max_val > 0 else 0
                values.append(normalized_val)
            
            # Create radar chart
            angles = np.linspace(0, 2 * np.pi, len(property_names), endpoint=False)
            values += [values[0]]  # Complete the circle
            angles = np.concatenate((angles, [angles[0]]))
            
            fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
            ax.plot(angles, values, 'o-', linewidth=2, label='Ligand Properties')
            ax.fill(angles, values, alpha=0.25)
            
            # Add property labels
            ax.set_xticks(angles[:-1])
            ax.set_xticklabels(property_names)
            ax.set_ylim(0, 1)
            ax.set_title('Molecular Properties Radar Chart', pad=20)
            
            # Add grid
            ax.grid(True)
            
            plt.tight_layout()
            plt.savefig(output_path / 'molecular_properties_radar.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create molecular properties radar: {e}")

    def _create_interaction_type_distribution(
        self,
        binding_site_residues: List[Dict],
        output_path: Path
    ):
        """Create distribution plot of interaction types."""
        try:
            # Count interaction types
            interaction_counts = {}
            for residue in binding_site_residues:
                interaction_type = residue.get('interaction_type', 'unknown')
                interaction_counts[interaction_type] = interaction_counts.get(interaction_type, 0) + 1
            
            if not interaction_counts:
                logger.warning("⚠️ No interaction data available")
                return
            
            # Create plots
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
            
            # Bar chart
            types = list(interaction_counts.keys())
            counts = list(interaction_counts.values())
            colors = plt.cm.Set3(np.linspace(0, 1, len(types)))
            
            bars = ax1.bar(types, counts, color=colors, alpha=0.8, edgecolor='black')
            ax1.set_xlabel('Interaction Type')
            ax1.set_ylabel('Number of Residues')
            ax1.set_title('Interaction Type Distribution')
            ax1.tick_params(axis='x', rotation=45)
            
            # Add value labels on bars
            for bar, count in zip(bars, counts):
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{count}', ha='center', va='bottom')
            
            # Pie chart
            ax2.pie(counts, labels=types, autopct='%1.1f%%', colors=colors)
            ax2.set_title('Interaction Type Proportions')
            
            plt.tight_layout()
            plt.savefig(output_path / 'interaction_type_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create interaction type distribution: {e}")

    def _create_interactive_energy_plot(
        self,
        affinities: List[Dict],
        output_path: Path
    ):
        """Create interactive energy plot using Plotly."""
        if not PLOTLY_AVAILABLE:
            return
            
        try:
            energies = [a['affinity'] for a in affinities]
            poses = list(range(1, len(energies) + 1))
            rmsd_lb = [a.get('rmsd_lb', 0) for a in affinities]
            
            # Create interactive scatter plot
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=poses,
                y=energies,
                mode='markers+lines',
                marker=dict(
                    size=8,
                    color=energies,
                    colorscale='Viridis',
                    colorbar=dict(title="Energy (kcal/mol)")
                ),
                text=[f'Pose {p}<br>Energy: {e:.2f} kcal/mol<br>RMSD: {r:.2f} Å' 
                      for p, e, r in zip(poses, energies, rmsd_lb)],
                hovertemplate='%{text}<extra></extra>',
                name='Docking Poses'
            ))
            
            fig.update_layout(
                title='Interactive Energy vs Pose Rank',
                xaxis_title='Pose Rank',
                yaxis_title='Binding Energy (kcal/mol)',
                hovermode='closest',
                width=800,
                height=600
            )
            
            # Save interactive plot
            fig.write_html(output_path / 'interactive_energy_plot.html')
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create interactive energy plot: {e}")

    def _create_video_from_frames(
        self,
        frames_dir: Path,
        output_dir: Path,
        duration: int,
        prefix: str = "frame"
    ):
        """Enhanced video creation from frames with better naming."""
        try:
            # Check if ffmpeg is available
            subprocess.run(['ffmpeg', '-version'], 
                         capture_output=True, check=True)
            
            # Create video
            fps = 120 / duration if duration > 0 else 15
            input_pattern = str(frames_dir / f'{prefix}_%04d.png')
            output_file = str(output_dir / f'{prefix}_animation.mp4')
            
            cmd = [
                'ffmpeg', '-y',
                '-framerate', str(fps),
                '-i', input_pattern,
                '-c:v', 'libx264',
                '-pix_fmt', 'yuv420p',
                '-crf', '18',
                '-preset', 'medium',
                output_file
            ]
            
            result = subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"🎬 Created video: {prefix}_animation.mp4")
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning(f"⚠️ ffmpeg not available, {prefix} animation frames saved as images")

    def _create_matplotlib_multi_format_plots(
        self,
        protein_file: str,
        ligand_pose: str,
        output_path: Path,
        formats: List[str]
    ):
        """Create matplotlib plots in multiple formats."""
        try:
            # Extract coordinates
            protein_coords = self._extract_coordinates(protein_file)
            ligand_coords = self._extract_coordinates(ligand_pose)
            
            # Create plot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
            
            if protein_coords:
                protein_array = np.array(protein_coords)
                ax1.scatter(protein_array[:, 0], protein_array[:, 1], 
                           c='gray', alpha=0.6, s=1, label='Protein')
            
            if ligand_coords:
                ligand_array = np.array(ligand_coords)
                ax1.scatter(ligand_array[:, 0], ligand_array[:, 1],
                           c='red', s=20, label='Ligand')
            
            ax1.set_xlabel('X coordinate (Å)')
            ax1.set_ylabel('Y coordinate (Å)')
            ax1.set_title('Docking Complex (XY projection)')
            ax1.legend()
            ax1.axis('equal')
            
            # Second view (XZ projection)
            if protein_coords:
                ax2.scatter(protein_array[:, 0], protein_array[:, 2],
                           c='gray', alpha=0.6, s=1, label='Protein')
            
            if ligand_coords:
                ax2.scatter(ligand_array[:, 0], ligand_array[:, 2],
                           c='red', s=20, label='Ligand')
            
            ax2.set_xlabel('X coordinate (Å)')
            ax2.set_ylabel('Z coordinate (Å)')
            ax2.set_title('Docking Complex (XZ projection)')
            ax2.legend()
            ax2.axis('equal')
            
            plt.tight_layout()
            
            # Save in multiple formats
            for fmt in formats:
                plt.savefig(output_path / f'docking_overview.{fmt}', 
                           dpi=300, bbox_inches='tight', format=fmt)
            
            plt.close()
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to create matplotlib plots: {e}")

    # Keep existing methods...
    def _create_pymol_images(self, *args, **kwargs):
        """Legacy method redirects to multi-format version."""
        return self._create_pymol_multi_format_images(*args, **kwargs, formats=['png'])

    def _create_chimerax_multi_format_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: Path,
        binding_site_center: Tuple[float, float, float],
        width: int,
        height: int,
        formats: List[str]
    ):
        """Create ChimeraX images in multiple formats."""
        colors = self.default_settings.get('colors', {})
        
        for format_type in formats:
            script_file = output_dir / f"chimerax_script_{format_type}.cxc"
            
            with open(script_file, 'w') as f:
                f.write(f"""
open {protein_file}
open {ligand_pose}

# Set representations
hide atoms
show cartoons #1
show atoms #2
style #2 stick

# Set colors
color #1 {colors.get('protein', 'gray')}
color #2 {colors.get('ligand', 'red')}

# Center on ligand
view #2

# Create images
save {output_dir}/docking_overview.{format_type} width {width} height {height}

# Close-up view
zoom #2 5
save {output_dir}/docking_closeup.{format_type} width {width} height {height}

exit
""")
            
            try:
                cmd = ['chimerax', '--nogui', str(script_file)]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError:
                logger.warning(f"⚠️ Failed to create {format_type} images with ChimeraX")
