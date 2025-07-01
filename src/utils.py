"""
Utility functions for the molecular docking pipeline.
"""

import os
import sys
from pathlib import Path
from typing import Dict, Any

from loguru import logger
import yaml


def setup_logging(level: str = "INFO"):
    """Setup loguru logging with custom format."""
    logger.remove()  # Remove default handler
    
    # Console handler with colors
    logger.add(
        sys.stdout,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level,
        colorize=True
    )
    
    # File handler for debugging
    logger.add(
        "logs/pipeline.log",
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} | {message}",
        level="DEBUG",
        rotation="10 MB",
        retention="7 days"
    )


def validate_inputs(protein_file: str, ligand_file: str):
    """Validate input files exist and have correct extensions."""
    protein_path = Path(protein_file)
    ligand_path = Path(ligand_file)
    
    # Check if files exist
    if not protein_path.exists():
        raise FileNotFoundError(f"Protein file not found: {protein_file}")
    
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    # Check protein file extension
    if protein_path.suffix.lower() != '.pdb':
        raise ValueError(f"Protein file must be in .pdb format, got: {protein_path.suffix}")
    
    # Check ligand file extension
    valid_ligand_extensions = {'.mol2', '.sdf', '.pdb'}
    if ligand_path.suffix.lower() not in valid_ligand_extensions:
        raise ValueError(f"Ligand file must be in {valid_ligand_extensions} format, got: {ligand_path.suffix}")
    
    logger.info(f"✅ Input validation passed")


def create_output_structure(output_dir: str, overwrite: bool = False) -> Dict[str, str]:
    """Create organized output directory structure."""
    base_path = Path(output_dir)
    
    # Check if directory exists and handle overwrite
    if base_path.exists() and not overwrite:
        import time
        timestamp = int(time.time())
        base_path = Path(f"{output_dir}_{timestamp}")
        logger.warning(f"⚠️  Output directory exists, using: {base_path}")
    
    # Define subdirectories
    subdirs = {
        'root': base_path,
        'images': base_path / 'images',
        'animation': base_path / 'animation', 
        'poses': base_path / 'poses',
        'logs': base_path / 'logs',
        'temp': base_path / 'temp',
        'visualizations': base_path / 'visualizations'
    }
    
    # Create directories
    for dir_path in subdirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)
    
    # Create logs directory at root level too
    Path('logs').mkdir(exist_ok=True)
    
    logger.info(f"📁 Created output structure in: {base_path.absolute()}")
    
    return {key: str(path) for key, path in subdirs.items()}


def load_config(config_file: str = "config.yaml") -> Dict[str, Any]:
    """Load configuration from YAML file."""
    config_path = Path(config_file)
    
    if not config_path.exists():
        # Create default config
        default_config = {
            'vina': {
                'exhaustiveness': 8,
                'num_modes': 9,
                'energy_range': 3.0
            },
            'visualization': {
                'image_width': 1200,
                'image_height': 900,
                'animation_frames': 120,
                'animation_duration': 8
            },
            'preprocessing': {
                'add_hydrogens': True,
                'remove_waters': True,
                'minimize_ligand': True
            }
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(default_config, f, default_flow_style=False)
        
        logger.info(f"📄 Created default config: {config_path}")
        return default_config
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    logger.info(f"📄 Loaded config from: {config_path}")
    return config


def check_dependencies():
    """Check if required software dependencies are available."""
    import subprocess
    import shutil
    
    dependencies = {
        'vina': 'AutoDock Vina',
        'obabel': 'Open Babel',
        'pymol': 'PyMOL'
    }
    
    missing = []
    
    for cmd, name in dependencies.items():
        if shutil.which(cmd) is None:
            missing.append(name)
    
    if missing:
        logger.warning(f"⚠️  Missing dependencies: {', '.join(missing)}")
        logger.warning("Some features may not work correctly")
    else:
        logger.info("✅ All dependencies found")
    
    return len(missing) == 0


def format_time(seconds: float) -> str:
    """Format time duration in human readable format."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


def get_file_size(file_path: str) -> str:
    """Get human readable file size."""
    size = Path(file_path).stat().st_size
    
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    
    return f"{size:.1f} TB"


def clean_temp_files(temp_dir: str):
    """Clean up temporary files."""
    import shutil
    
    temp_path = Path(temp_dir)
    if temp_path.exists():
        shutil.rmtree(temp_path)
        logger.info(f"🧹 Cleaned temporary files: {temp_dir}")


class ProgressReporter:
    """Helper class for reporting progress with timing."""
    
    def __init__(self, total_steps: int):
        self.total_steps = total_steps
        self.current_step = 0
        self.start_time = None
        
    def start(self):
        """Start timing."""
        import time
        self.start_time = time.time()
        
    def step(self, description: str):
        """Report progress for current step."""
        import time
        
        if self.start_time is None:
            self.start()
            
        self.current_step += 1
        elapsed = time.time() - self.start_time
        
        logger.info(f"[{self.current_step}/{self.total_steps}] {description} ({format_time(elapsed)} elapsed)")
        
    def finish(self):
        """Report completion."""
        import time
        
        if self.start_time is not None:
            total_time = time.time() - self.start_time
            logger.success(f"✅ Completed all steps in {format_time(total_time)}")


# Configuration defaults
DEFAULT_VINA_PARAMS = {
    'exhaustiveness': 8,
    'num_modes': 9,
    'energy_range': 3.0
}

DEFAULT_BOX_SIZE = {
    'x': 22.5,
    'y': 22.5, 
    'z': 22.5
}

# File type mappings
SUPPORTED_PROTEIN_FORMATS = {'.pdb'}
SUPPORTED_LIGAND_FORMATS = {'.mol2', '.sdf', '.pdb'}

# Visualization settings
VIZ_SETTINGS = {
    'protein_representation': 'cartoon',
    'ligand_representation': 'sticks',
    'binding_site_radius': 8.0,
    'image_ray_trace': True,
    'animation_fps': 15
} 