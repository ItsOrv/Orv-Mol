#!/usr/bin/env python3
"""
Basic functionality test for Orv-Mol without heavy dependencies.
Tests core logic, configuration, and basic operations.
"""

import sys
import os
import tempfile
import unittest
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

def mock_heavy_imports():
    """Mock heavy imports that might not be available."""
    import types
    
    # Mock numpy
    if 'numpy' not in sys.modules:
        numpy_mock = types.ModuleType('numpy')
        numpy_mock.array = lambda x: x
        numpy_mock.mean = lambda x, axis=None: sum(x) / len(x) if hasattr(x, '__len__') else x
        numpy_mock.min = lambda x, axis=None: min(x) if hasattr(x, '__iter__') else x
        numpy_mock.max = lambda x, axis=None: max(x) if hasattr(x, '__iter__') else x
        numpy_mock.std = lambda x: 0
        numpy_mock.median = lambda x: sorted(x)[len(x)//2] if hasattr(x, '__len__') else x
        numpy_mock.arange = lambda start, stop, step=1: list(range(int(start), int(stop), int(step)))
        numpy_mock.sum = lambda x: sum(x) if hasattr(x, '__iter__') else x
        numpy_mock.prod = lambda x: 1
        numpy_mock.ndarray = list
        numpy_mock.linspace = lambda start, stop, num: [start + (stop-start)*i/(num-1) for i in range(num)]
        numpy_mock.concatenate = lambda arrays: sum(arrays, [])
        sys.modules['numpy'] = numpy_mock
    
    # Mock pandas
    if 'pandas' not in sys.modules:
        pandas_mock = types.ModuleType('pandas')
        pandas_mock.DataFrame = dict
        sys.modules['pandas'] = pandas_mock
    
    # Mock loguru
    if 'loguru' not in sys.modules:
        loguru_mock = types.ModuleType('loguru')
        class MockLogger:
            def info(self, msg): print(f"INFO: {msg}")
            def warning(self, msg): print(f"WARNING: {msg}")
            def error(self, msg): print(f"ERROR: {msg}")
            def success(self, msg): print(f"SUCCESS: {msg}")
        loguru_mock.logger = MockLogger()
        sys.modules['loguru'] = loguru_mock

# Apply mocks before importing our modules
mock_heavy_imports()


class TestBasicFunctionality(unittest.TestCase):
    """Test basic functionality without heavy dependencies."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_config_loading(self):
        """Test configuration loading functionality."""
        try:
            from src.utils import load_config
            
            # Test with default config creation
            config = load_config(str(Path(self.test_dir) / "test_config.yaml"))
            self.assertIsInstance(config, dict)
            
        except ImportError as e:
            self.skipTest(f"Config loading test skipped due to missing dependency: {e}")
    
    def test_output_structure_creation(self):
        """Test output directory structure creation."""
        try:
            from src.utils import create_output_structure
            
            output_dir = Path(self.test_dir) / "test_output"
            result = create_output_structure(str(output_dir), overwrite=True)
            
            # Check that function returns expected structure
            self.assertIsInstance(result, dict)
            self.assertIn('root', result)
            self.assertIn('visualizations', result)
            
            # Check directories are created
            self.assertTrue(Path(result['root']).exists())
            
        except ImportError as e:
            self.skipTest(f"Output structure test skipped due to missing dependency: {e}")
    
    def test_file_validation(self):
        """Test file validation logic."""
        try:
            from src.utils import validate_inputs
            
            # Create test files
            protein_file = Path(self.test_dir) / "test_protein.pdb"
            ligand_file = Path(self.test_dir) / "test_ligand.mol2"
            
            protein_file.write_text("ATOM  1  N   ALA A   1    1.000   1.000   1.000  1.00 20.00      N")
            ligand_file.write_text("@<TRIPOS>MOLECULE\ntest\n")
            
            # Test validation (should not raise exception)
            try:
                validate_inputs(str(protein_file), str(ligand_file))
                validation_passed = True
            except Exception:
                validation_passed = False
            
            self.assertTrue(validation_passed)
            
        except ImportError as e:
            self.skipTest(f"File validation test skipped due to missing dependency: {e}")
    
    def test_parser_initialization(self):
        """Test OutputParser initialization."""
        try:
            from src.output_parser import OutputParser
            
            parser = OutputParser()
            
            # Check basic attributes exist
            self.assertTrue(hasattr(parser, 'interaction_cutoffs'))
            self.assertTrue(hasattr(parser, 'amino_acid_properties'))
            
            # Check amino acid properties structure
            self.assertIn('hydrophobic', parser.amino_acid_properties)
            self.assertIn('ALA', parser.amino_acid_properties['hydrophobic'])
            
        except ImportError as e:
            self.skipTest(f"Parser initialization test skipped due to missing dependency: {e}")
    
    def test_visualization_engine_initialization(self):
        """Test VisualizationEngine initialization."""
        try:
            from src.visualization import VisualizationEngine
            
            # Test basic initialization
            viz_engine = VisualizationEngine()
            self.assertTrue(hasattr(viz_engine, 'default_settings'))
            
            # Test initialization with config
            test_config = {
                'visualization': {
                    'image_width': 1000,
                    'colors': {'protein': 'blue'}
                }
            }
            viz_engine_config = VisualizationEngine(config=test_config)
            self.assertEqual(viz_engine_config.default_settings['image_width'], 1000)
            
        except ImportError as e:
            self.skipTest(f"Visualization engine test skipped due to missing dependency: {e}")
    
    def test_dock_script_structure(self):
        """Test main dock.py script structure."""
        # Test that dock.py can be parsed
        dock_file = Path('dock.py')
        self.assertTrue(dock_file.exists())
        
        # Check for main function
        with open(dock_file, 'r') as f:
            content = f.read()
        
        self.assertIn('def main()', content)
        self.assertIn('if __name__ == "__main__"', content)
    
    def test_binding_analysis_structure(self):
        """Test binding analysis methods structure."""
        try:
            from src.output_parser import OutputParser
            
            parser = OutputParser()
            
            # Check analyze_binding_interactions method exists
            self.assertTrue(hasattr(parser, 'analyze_binding_interactions'))
            self.assertTrue(callable(getattr(parser, 'analyze_binding_interactions')))
            
            # Test with minimal data
            mock_results = {
                'affinities': [{'affinity': -8.5}, {'affinity': -7.2}]
            }
            
            analysis = parser.analyze_binding_interactions(
                detailed_results=mock_results,
                protein_file=None,
                ligand_file=None
            )
            
            # Check basic structure
            self.assertIsInstance(analysis, dict)
            self.assertIn('interaction_summary', analysis)
            
        except ImportError as e:
            self.skipTest(f"Binding analysis test skipped due to missing dependency: {e}")
    
    def test_energy_analysis(self):
        """Test energy analysis functionality."""
        try:
            from src.output_parser import OutputParser
            
            parser = OutputParser()
            
            # Test binding affinity distribution analysis
            mock_affinities = [
                {'affinity': -9.2},
                {'affinity': -8.5},
                {'affinity': -7.8},
                {'affinity': -7.1}
            ]
            
            analysis = parser._analyze_binding_affinity_distribution(mock_affinities)
            
            # Check basic analysis results
            self.assertIn('best_affinity', analysis)
            self.assertIn('binding_strength', analysis)
            self.assertEqual(analysis['best_affinity'], -9.2)
            
        except ImportError as e:
            self.skipTest(f"Energy analysis test skipped due to missing dependency: {e}")
    
    def test_config_file_validity(self):
        """Test that config.yaml is valid."""
        config_file = Path('config.yaml')
        self.assertTrue(config_file.exists())
        
        try:
            import yaml
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            
            # Check main sections exist
            expected_sections = ['docking', 'preprocessing', 'visualization', 'analysis']
            for section in expected_sections:
                self.assertIn(section, config)
                
        except ImportError:
            # YAML not available, skip detailed validation
            pass
    
    def test_requirements_file(self):
        """Test requirements.txt structure."""
        req_file = Path('requirements.txt')
        self.assertTrue(req_file.exists())
        
        with open(req_file, 'r') as f:
            content = f.read()
        
        # Check for essential dependencies
        essential_deps = ['numpy', 'matplotlib', 'yaml', 'pathlib']
        for dep in essential_deps:
            self.assertIn(dep, content.lower())

def run_basic_functionality_tests():
    """Run basic functionality tests."""
    print("🧪 Running Basic Functionality Tests")
    print("=" * 50)
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(TestBasicFunctionality)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 50)
    print("📊 Basic Functionality Test Results")
    print("=" * 50)
    
    total_tests = result.testsRun
    failures = len(result.failures)
    errors = len(result.errors)
    skipped = len(result.skipped)
    passed = total_tests - failures - errors - skipped
    
    print(f"Total tests: {total_tests}")
    print(f"Passed: {passed}")
    print(f"Failed: {failures}")
    print(f"Errors: {errors}")
    print(f"Skipped: {skipped}")
    
    if result.wasSuccessful():
        print("\n🎉 All basic functionality tests passed!")
        print("✨ Core logic and structure are working correctly")
        return True
    else:
        print(f"\n⚠️ Some tests failed or had errors")
        if result.failures:
            print("\nFailures:")
            for test, traceback in result.failures:
                print(f"- {test}: {traceback}")
        if result.errors:
            print("\nErrors:")
            for test, traceback in result.errors:
                print(f"- {test}: {traceback}")
        return False

if __name__ == "__main__":
    success = run_basic_functionality_tests()
    sys.exit(0 if success else 1) 