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
import numpy as np
from loguru import logger


class VisualizationEngine:
    """Handles molecular visualization using PyMOL and ChimeraX."""
    
    def __init__(self):
        self.pymol_available = self._check_pymol()
        self.chimerax_available = self._check_chimerax()
        
        if not self.pymol_available and not self.chimerax_available:
            logger.warning("⚠️ Neither PyMOL nor ChimeraX found. Visualization features will be limited.")
    
    def create_docking_images(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str,
        binding_site_center: Tuple[float, float, float],
        image_width: int = 1200,
        image_height: int = 900
    ):
        """
        Create publication-ready images of the docking result.
        
        Args:
            protein_file: Path to protein PDB file
            ligand_pose: Path to ligand pose PDB file  
            output_dir: Directory to save images
            binding_site_center: Coordinates of binding site center
            image_width: Image width in pixels
            image_height: Image height in pixels
        """
        logger.info("🖼️ Creating publication-ready images...")
        
        output_path = Path(output_dir)
        
        try:
            if self.pymol_available:
                self._create_pymol_images(
                    protein_file, ligand_pose, output_path, 
                    binding_site_center, image_width, image_height
                )
            elif self.chimerax_available:
                self._create_chimerax_images(
                    protein_file, ligand_pose, output_path,
                    binding_site_center, image_width, image_height
                )
            else:
                logger.warning("⚠️ No visualization software available, creating basic plots")
                self._create_matplotlib_plot(protein_file, ligand_pose, output_path)
            
            logger.success("✅ Images created successfully")
            
        except Exception as e:
            logger.error(f"❌ Failed to create images: {e}")
            raise
    
    def create_3d_animation(
        self,
        protein_file: str,
        ligand_pose: str,
        output_dir: str,
        binding_site_center: Tuple[float, float, float],
        frames: int = 120,
        duration: int = 8
    ):
        """
        Create 3D rotation animation of the docking complex.
        
        Args:
            protein_file: Path to protein PDB file
            ligand_pose: Path to ligand pose PDB file
            output_dir: Directory to save animation
            binding_site_center: Coordinates of binding site center
            frames: Number of animation frames
            duration: Animation duration in seconds
        """
        logger.info("🎬 Creating 3D rotation animation...")
        
        output_path = Path(output_dir)
        
        try:
            if self.pymol_available:
                self._create_pymol_animation(
                    protein_file, ligand_pose, output_path,
                    binding_site_center, frames, duration
                )
            elif self.chimerax_available:
                self._create_chimerax_animation(
                    protein_file, ligand_pose, output_path,
                    binding_site_center, frames, duration
                )
            else:
                logger.warning("⚠️ Animation requires PyMOL or ChimeraX")
                return
            
            logger.success("✅ Animation created successfully")
            
        except Exception as e:
            logger.error(f"❌ Failed to create animation: {e}")
            raise
    
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
