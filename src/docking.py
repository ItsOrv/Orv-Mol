"""
Docking module for running AutoDock Vina molecular docking.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

from .logger_config import logger


class DockingEngine:
    """Handles AutoDock Vina docking operations."""
    
    def __init__(self):
        self.vina_executable = self._find_vina_executable()
    
    def run_vina_docking(
        self,
        receptor: str,
        ligand: str,
        box_params: Dict[str, float],
        output_dir: str,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0
    ) -> Dict[str, str]:
        """
        Run AutoDock Vina docking simulation.
        
        Args:
            receptor: Path to receptor PDBQT file
            ligand: Path to ligand PDBQT file
            box_params: Dictionary with grid box parameters
            output_dir: Directory to save output files
            exhaustiveness: Vina exhaustiveness parameter
            num_modes: Number of binding modes to generate
            energy_range: Energy range for binding modes (kcal/mol)
            
        Returns:
            Dictionary with paths to output files
        """
        logger.info("⚙️ Starting AutoDock Vina docking...")
        
        output_path = Path(output_dir)
        receptor_path = Path(receptor)
        ligand_path = Path(ligand)
        
        # Generate output filenames
        output_pdbqt = output_path / f"{ligand_path.stem}_docked.pdbqt"
        config_file = output_path / "vina_config.txt"
        
        try:
            # Create Vina configuration file
            self._create_vina_config(
                config_file=config_file,
                receptor=receptor,
                ligand=ligand,
                output=output_pdbqt,
                box_params=box_params,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                energy_range=energy_range
            )
            
            # Run Vina docking
            self._run_vina_command(config_file)
            
            # Verify output files exist
            if not output_pdbqt.exists():
                raise FileNotFoundError(f"Vina output file not found: {output_pdbqt}")
            
            # Note: Vina v1.2.5+ doesn't create separate log files anymore
            # The output is embedded in the PDBQT file
            
            logger.success(f"✅ Docking completed successfully!")
            logger.info(f"📄 Output PDBQT: {output_pdbqt}")
            
            return {
                'output_pdbqt': str(output_pdbqt),
                'log_file': str(output_pdbqt),  # Use PDBQT file as log source
                'config_file': str(config_file)
            }
            
        except Exception as e:
            logger.error(f"❌ Docking failed: {e}")
            raise
    
    def run_quick_validation(
        self,
        receptor: str,
        ligand: str,
        box_params: Dict[str, float]
    ) -> bool:
        """
        Run a quick validation docking to check if parameters are valid.
        
        Args:
            receptor: Path to receptor PDBQT file
            ligand: Path to ligand PDBQT file
            box_params: Dictionary with grid box parameters
            
        Returns:
            True if validation passes, False otherwise
        """
        logger.info("🔍 Running quick validation docking...")
        
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                
                # Create temporary output files
                output_pdbqt = temp_path / "validation_output.pdbqt"
                config_file = temp_path / "validation_config.txt"
                
                # Create minimal config for validation
                self._create_vina_config(
                    config_file=config_file,
                    receptor=receptor,
                    ligand=ligand,
                    output=output_pdbqt,
                    box_params=box_params,
                    exhaustiveness=1,  # Minimal exhaustiveness for speed
                    num_modes=1,       # Only one mode
                    energy_range=1.0
                )
                
                # Run quick validation
                cmd = [self.vina_executable, '--config', str(config_file)]
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60  # 1 minute timeout
                )
                
                success = result.returncode == 0 and output_pdbqt.exists()
                
                if success:
                    logger.info("✅ Validation passed")
                else:
                    logger.warning(f"⚠️ Validation failed: {result.stderr}")
                
                return success
                
        except Exception as e:
            logger.warning(f"⚠️ Validation error: {e}")
            return False
    
    def _find_vina_executable(self) -> str:
        """Find AutoDock Vina executable."""
        import shutil
        
        # Try common names
        vina_names = ['vina', 'vina1', 'autodock_vina', 'AutoDock_Vina']
        
        for name in vina_names:
            executable = shutil.which(name)
            if executable:
                logger.info(f"🔍 Found Vina executable: {executable}")
                return executable
        
        # Check common installation paths
        common_paths = [
            '/usr/local/bin/vina',
            '/usr/bin/vina',
            '/opt/vina/bin/vina',
            '/Applications/AutoDock_Vina.app/Contents/MacOS/vina'  # macOS
        ]
        
        for path in common_paths:
            if Path(path).exists():
                logger.info(f"🔍 Found Vina at: {path}")
                return path
        
        # If not found, assume it's in PATH
        logger.warning("⚠️ Vina executable not found in common locations, assuming 'vina' is in PATH")
        return 'vina'
    
    def _create_vina_config(
        self,
        config_file: Path,
        receptor: str,
        ligand: str,
        output: Path,
        box_params: Dict[str, float],
        exhaustiveness: int,
        num_modes: int,
        energy_range: float
    ):
        """Create AutoDock Vina configuration file."""
        with open(config_file, 'w') as f:
            f.write(f"receptor = {receptor}\n")
            f.write(f"ligand = {ligand}\n")
            f.write(f"out = {output}\n")
            f.write("\n")
            f.write(f"center_x = {box_params['center_x']:.3f}\n")
            f.write(f"center_y = {box_params['center_y']:.3f}\n")
            f.write(f"center_z = {box_params['center_z']:.3f}\n")
            f.write("\n")
            f.write(f"size_x = {box_params['size_x']:.3f}\n")
            f.write(f"size_y = {box_params['size_y']:.3f}\n")
            f.write(f"size_z = {box_params['size_z']:.3f}\n")
            f.write("\n")
            f.write(f"exhaustiveness = {exhaustiveness}\n")
            f.write(f"num_modes = {num_modes}\n")
            f.write(f"energy_range = {energy_range:.1f}\n")
            f.write("\n")
            f.write("# Additional parameters\n")
            f.write("cpu = 0  # Use all available CPUs\n")
            f.write("seed = 42  # For reproducibility\n")
        
        logger.info(f"📄 Created Vina config: {config_file}")
    
    def _run_vina_command(self, config_file: Path):
        """Run the AutoDock Vina command."""
        cmd = [self.vina_executable, '--config', str(config_file)]
        
        logger.info(f"🚀 Running command: {' '.join(cmd)}")
        
        try:
            # Run Vina with real-time output
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
            )
            
            # Log output in real-time
            for line in process.stdout:
                line = line.strip()
                if line:
                    if 'ERROR' in line.upper() or 'FATAL' in line.upper():
                        logger.error(f"Vina: {line}")
                    elif 'WARNING' in line.upper():
                        logger.warning(f"Vina: {line}")
                    else:
                        logger.debug(f"Vina: {line}")
            
            # Wait for completion
            return_code = process.wait()
            
            if return_code != 0:
                raise subprocess.CalledProcessError(
                    return_code, cmd, f"Vina failed with return code {return_code}"
                )
            
            logger.info("🎯 Vina execution completed successfully")
            
        except subprocess.TimeoutExpired:
            logger.error("❌ Vina execution timed out")
            process.kill()
            raise
        except FileNotFoundError:
            logger.error(f"❌ Vina executable not found: {self.vina_executable}")
            logger.error("Please install AutoDock Vina and ensure it's in your PATH")
            raise
        except subprocess.CalledProcessError as e:
            logger.error(f"❌ Vina execution failed: {e}")
            raise
    
    def estimate_runtime(
        self,
        exhaustiveness: int,
        num_modes: int,
        ligand_size: Optional[int] = None
    ) -> str:
        """
        Estimate docking runtime based on parameters.
        
        Args:
            exhaustiveness: Vina exhaustiveness parameter
            num_modes: Number of binding modes
            ligand_size: Number of rotatable bonds (optional)
            
        Returns:
            Human-readable time estimate
        """
        # Basic estimation formula (very rough)
        base_time = 30  # seconds
        exhaustiveness_factor = exhaustiveness / 8  # Normalize to default
        modes_factor = num_modes / 9  # Normalize to default
        
        if ligand_size:
            size_factor = max(1, ligand_size / 10)  # Assume 10 bonds as baseline
        else:
            size_factor = 1
        
        estimated_seconds = base_time * exhaustiveness_factor * modes_factor * size_factor
        
        if estimated_seconds < 60:
            return f"~{estimated_seconds:.0f} seconds"
        elif estimated_seconds < 3600:
            return f"~{estimated_seconds/60:.0f} minutes"
        else:
            return f"~{estimated_seconds/3600:.1f} hours"
    
    def check_vina_version(self) -> str:
        """Check AutoDock Vina version."""
        try:
            result = subprocess.run(
                [self.vina_executable, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                version_info = result.stdout.strip()
                logger.info(f"🔍 Vina version: {version_info}")
                return version_info
            else:
                # Try alternative version command
                result = subprocess.run(
                    [self.vina_executable, '--help'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if 'AutoDock Vina' in result.stdout:
                    logger.info("🔍 Vina detected (version unknown)")
                    return "Unknown version"
            
        except Exception as e:
            logger.warning(f"⚠️ Could not determine Vina version: {e}")
        
        return "Unknown"
    
    def validate_inputs(self, receptor: str, ligand: str) -> bool:
        """
        Validate that input files are suitable for Vina.
        
        Args:
            receptor: Path to receptor PDBQT file
            ligand: Path to ligand PDBQT file
            
        Returns:
            True if inputs are valid, False otherwise
        """
        issues = []
        
        # Check if files exist
        if not Path(receptor).exists():
            issues.append(f"Receptor file not found: {receptor}")
        
        if not Path(ligand).exists():
            issues.append(f"Ligand file not found: {ligand}")
        
        # Check file extensions
        if not receptor.endswith('.pdbqt'):
            issues.append(f"Receptor must be PDBQT format: {receptor}")
        
        if not ligand.endswith('.pdbqt'):
            issues.append(f"Ligand must be PDBQT format: {ligand}")
        
        # Check file contents
        try:
            with open(receptor, 'r') as f:
                receptor_content = f.read()
                if not any(line.startswith('ATOM') for line in receptor_content.split('\n')):
                    issues.append("Receptor file contains no ATOM records")
        except Exception as e:
            issues.append(f"Could not read receptor file: {e}")
        
        try:
            with open(ligand, 'r') as f:
                ligand_content = f.read()
                if not any(line.startswith('ATOM') for line in ligand_content.split('\n')):
                    issues.append("Ligand file contains no ATOM records")
                
                # Check for ROOT/ENDROOT (flexibility information)
                if 'ROOT' not in ligand_content:
                    logger.warning("⚠️ Ligand may not have flexibility information (no ROOT found)")
        except Exception as e:
            issues.append(f"Could not read ligand file: {e}")
        
        if issues:
            logger.error("❌ Input validation failed:")
            for issue in issues:
                logger.error(f"  • {issue}")
            return False
        
        logger.info("✅ Input validation passed")
        return True 