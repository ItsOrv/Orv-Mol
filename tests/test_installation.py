#!/usr/bin/env python3
"""
Test script to verify the Automated Molecular Docking Pipeline installation.
This script checks dependencies and performs basic functionality tests.
"""

import sys
import subprocess
import importlib
from pathlib import Path

# Colors for output
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def print_success(msg):
    print(f"{Colors.GREEN}✅ {msg}{Colors.ENDC}")

def print_error(msg):
    print(f"{Colors.RED}❌ {msg}{Colors.ENDC}")

def print_warning(msg):
    print(f"{Colors.YELLOW}⚠️  {msg}{Colors.ENDC}")

def print_info(msg):
    print(f"{Colors.BLUE}ℹ️  {msg}{Colors.ENDC}")

def print_header(msg):
    print(f"\n{Colors.BOLD}{Colors.BLUE}{msg}{Colors.ENDC}")
    print("=" * len(msg))

def test_python_version():
    """Test Python version compatibility"""
    print_header("Testing Python Version")
    
    version = sys.version_info
    print_info(f"Python version: {version.major}.{version.minor}.{version.micro}")
    
    if version.major == 3 and version.minor >= 8:
        print_success("Python version is compatible")
        return True
    else:
        print_error("Python 3.8+ required")
        return False

def test_python_imports():
    """Test essential Python package imports"""
    print_header("Testing Python Package Imports")
    
    required_packages = [
        ('numpy', 'numpy'),
        ('pandas', 'pandas'),
        ('scipy', 'scipy'),
        ('pathlib', 'pathlib'),
        ('matplotlib', 'matplotlib.pyplot'),
        ('loguru', 'loguru'),
        ('tqdm', 'tqdm'),
        ('click', 'click'),
        ('yaml', 'yaml'),
    ]
    
    optional_packages = [
        ('rdkit', 'rdkit.Chem'),
        ('Bio', 'Bio.PDB'),
        ('pymol', 'pymol'),
    ]
    
    success_count = 0
    
    # Test required packages
    for package_name, import_name in required_packages:
        try:
            importlib.import_module(import_name)
            print_success(f"{package_name} imported successfully")
            success_count += 1
        except ImportError:
            print_error(f"{package_name} import failed")
    
    # Test optional packages
    for package_name, import_name in optional_packages:
        try:
            importlib.import_module(import_name)
            print_success(f"{package_name} imported successfully (optional)")
        except ImportError:
            print_warning(f"{package_name} not available (optional)")
    
    return success_count == len(required_packages)

def test_external_tools():
    """Test external software dependencies"""
    print_header("Testing External Software")
    
    tools = [
        ('vina', 'AutoDock Vina', True),
        ('obabel', 'Open Babel', True),
        ('pymol', 'PyMOL', False),
        ('chimerax', 'ChimeraX', False),
        ('ffmpeg', 'FFmpeg', False),
        ('prepare_ligand4.py', 'AutoDockTools', False),
    ]
    
    success_count = 0
    required_count = 0
    
    for tool, name, required in tools:
        if required:
            required_count += 1
            
        try:
            # Test if tool is available
            result = subprocess.run([tool, '--help'], 
                                  capture_output=True, 
                                  timeout=5)
            if result.returncode == 0 or 'usage' in result.stdout.decode().lower():
                print_success(f"{name} is available")
                if required:
                    success_count += 1
            else:
                raise subprocess.CalledProcessError(result.returncode, tool)
                
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            if required:
                print_error(f"{name} not found (required)")
            else:
                print_warning(f"{name} not found (optional)")
    
    return success_count == required_count

def test_pipeline_modules():
    """Test pipeline module imports"""
    print_header("Testing Pipeline Modules")
    
    # Add src to path
    sys.path.insert(0, str(Path(__file__).parent / 'src'))
    
    modules = [
        'utils',
        'preprocessing',
        'docking',
        'output_parser',
        'visualization'
    ]
    
    success_count = 0
    
    for module in modules:
        try:
            importlib.import_module(module)
            print_success(f"{module} module imported successfully")
            success_count += 1
        except ImportError as e:
            print_error(f"{module} module import failed: {e}")
    
    return success_count == len(modules)

def test_file_structure():
    """Test project file structure"""
    print_header("Testing File Structure")
    
    required_files = [
        'dock.py',
        'requirements.txt',
        'config.yaml',
        'README.md',
        'src/__init__.py',
        'src/utils.py',
        'src/preprocessing.py',
        'src/docking.py',
        'src/output_parser.py',
        'src/visualization.py',
    ]
    
    missing_files = []
    
    for file_path in required_files:
        if Path(file_path).exists():
            print_success(f"{file_path} exists")
        else:
            print_error(f"{file_path} missing")
            missing_files.append(file_path)
    
    return len(missing_files) == 0

def test_basic_functionality():
    """Test basic pipeline functionality"""
    print_header("Testing Basic Functionality")
    
    try:
        # Add src to path
        sys.path.insert(0, str(Path(__file__).parent / 'src'))
        
        # Test utility functions
        from utils import setup_logging, create_output_structure
        
        print_info("Testing logging setup...")
        setup_logging("INFO")
        print_success("Logging setup successful")
        
        print_info("Testing output structure creation...")
        output_paths = create_output_structure("test_output")
        print_success("Output structure creation successful")
        
        # Clean up test directory
        import shutil
        if Path("test_output").exists():
            shutil.rmtree("test_output")
        
        # Test basic class initialization
        from preprocessing import PreprocessingEngine
        from docking import DockingEngine
        from output_parser import OutputParser
        from visualization import VisualizationEngine
        
        print_info("Testing class initialization...")
        PreprocessingEngine()
        DockingEngine()
        OutputParser()
        VisualizationEngine()
        print_success("All classes initialized successfully")
        
        return True
        
    except Exception as e:
        print_error(f"Basic functionality test failed: {e}")
        return False

def test_example_files():
    """Check for example files"""
    print_header("Checking Example Files")
    
    example_files = [
        'examples/protein.pdb',
        'examples/ligand.sdf',
        'examples/ligand.mol2',
    ]
    
    found_examples = 0
    
    for file_path in example_files:
        if Path(file_path).exists():
            print_success(f"{file_path} exists")
            found_examples += 1
        else:
            print_warning(f"{file_path} not found")
    
    if found_examples == 0:
        print_info("No example files found. You'll need to provide your own protein and ligand files.")
    
    return True  # Not critical for functionality

def run_all_tests():
    """Run all tests and provide summary"""
    print(f"{Colors.BOLD}{Colors.BLUE}")
    print("🧪 Automated Molecular Docking Pipeline - Installation Test")
    print("=" * 60)
    print(f"{Colors.ENDC}")
    
    tests = [
        ("Python Version", test_python_version),
        ("Python Imports", test_python_imports),
        ("External Tools", test_external_tools),
        ("Pipeline Modules", test_pipeline_modules),
        ("File Structure", test_file_structure),
        ("Basic Functionality", test_basic_functionality),
        ("Example Files", test_example_files),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print_error(f"Test '{test_name}' failed with exception: {e}")
            results.append((test_name, False))
    
    # Print summary
    print_header("Test Summary")
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        if result:
            print_success(f"{test_name}")
        else:
            print_error(f"{test_name}")
    
    print(f"\n{Colors.BOLD}Tests passed: {passed}/{total}{Colors.ENDC}")
    
    if passed == total:
        print_success("All tests passed! Pipeline is ready to use.")
        print_info("Try running: python dock.py --help")
        return True
    else:
        print_error("Some tests failed. Please check the installation.")
        
        if passed >= total - 2:  # Allow for optional dependencies
            print_warning("Core functionality appears to work, but some optional features may be limited.")
        
        return False

def main():
    """Main test function"""
    if len(sys.argv) > 1 and sys.argv[1] == '--quick':
        # Quick test - just essential checks
        success = (test_python_version() and 
                  test_python_imports() and 
                  test_pipeline_modules())
        
        if success:
            print_success("Quick test passed!")
        else:
            print_error("Quick test failed!")
        
        return success
    else:
        # Full test suite
        return run_all_tests()

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
