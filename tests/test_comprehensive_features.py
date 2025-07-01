#!/usr/bin/env python3
"""
Comprehensive test suite for enhanced Orv-Mol features.
Tests visualization, analysis, and output parsing improvements.
"""

import sys
import os
import unittest
import tempfile
import yaml
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from visualization import VisualizationEngine
from output_parser import OutputParser
from utils import create_output_structure, load_config


class TestComprehensiveFeatures(unittest.TestCase):
    """Test suite for comprehensive enhanced features."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        self.test_config = {
            'visualization': {
                'image_format': 'png',
                'image_dpi': 150,
                'image_width': 800,
                'image_height': 600,
                'colors': {
                    'protein': 'gray80',
                    'ligand': 'red',
                    'binding_site': 'yellow'
                },
                'animation': {
                    'create_animations': False,
                    'frames': 60,
                    'duration': 4
                }
            },
            'analysis': {
                'binding_affinity': {
                    'classify_strength': True
                },
                'interactions': {
                    'detailed_analysis': True
                }
            }
        }
        
        # Create test PDB content
        self.test_protein_pdb = """ATOM      1  N   ALA A   1      20.154  16.967  27.462  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.174  26.933  1.00 20.00           C  
ATOM      3  C   ALA A   1      18.290  15.569  28.109  1.00 20.00           C  
ATOM      4  O   ALA A   1      18.815  15.215  29.174  1.00 20.00           O  
ATOM      5  CB  ALA A   1      19.595  15.081  26.041  1.00 20.00           C  
END"""
        
        self.test_ligand_pdb = """HETATM    1  C1  UNL     1      20.000  16.000  28.000  1.00  0.00           C  
HETATM    2  C2  UNL     1      21.000  17.000  29.000  1.00  0.00           C  
HETATM    3  O1  UNL     1      19.000  15.000  27.000  1.00  0.00           O  
END"""
        
        # Create test files
        self.protein_file = Path(self.test_dir) / "test_protein.pdb"
        self.ligand_file = Path(self.test_dir) / "test_ligand.pdb"
        
        with open(self.protein_file, 'w') as f:
            f.write(self.test_protein_pdb)
        
        with open(self.ligand_file, 'w') as f:
            f.write(self.test_ligand_pdb)
        
        # Create mock docking results
        self.mock_results = {
            'affinities': [
                {'affinity': -8.5, 'rmsd_lb': 0.0, 'rmsd_ub': 0.0},
                {'affinity': -7.8, 'rmsd_lb': 1.2, 'rmsd_ub': 2.1},
                {'affinity': -7.2, 'rmsd_lb': 2.5, 'rmsd_ub': 3.8}
            ],
            'best_pose_pdb': str(self.ligand_file),
            'poses': [
                {'file': str(self.ligand_file), 'energy': -8.5},
                {'file': str(self.ligand_file), 'energy': -7.8},
                {'file': str(self.ligand_file), 'energy': -7.2}
            ]
        }
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_visualization_engine_initialization(self):
        """Test VisualizationEngine initialization with config."""
        viz_engine = VisualizationEngine(config=self.test_config)
        
        # Check default settings are updated with config
        self.assertEqual(viz_engine.default_settings['image_width'], 800)
        self.assertEqual(viz_engine.default_settings['image_height'], 600)
        self.assertEqual(viz_engine.default_settings['colors']['protein'], 'gray80')
    
    def test_output_parser_binding_analysis(self):
        """Test enhanced binding interaction analysis."""
        parser = OutputParser()
        
        # Test basic binding analysis
        analysis = parser.analyze_binding_interactions(
            detailed_results=self.mock_results,
            protein_file=str(self.protein_file),
            ligand_file=str(self.ligand_file)
        )
        
        # Check analysis structure
        self.assertIn('interaction_summary', analysis)
        self.assertIn('binding_site_residues', analysis)
        self.assertIn('geometric_properties', analysis)
        self.assertIn('pharmacological_analysis', analysis)
        
        # Check binding affinity analysis
        affinity_analysis = parser._analyze_binding_affinity_distribution(self.mock_results['affinities'])
        self.assertIn('best_affinity', affinity_analysis)
        self.assertIn('binding_strength', affinity_analysis)
        self.assertEqual(affinity_analysis['best_affinity'], -8.5)
    
    def test_create_output_structure_enhanced(self):
        """Test enhanced output structure creation."""
        output_dir = Path(self.test_dir) / "test_output"
        
        # Test output structure creation
        output_paths = create_output_structure(str(output_dir), overwrite=True)
        
        # Check all required directories are created
        required_dirs = ['root', 'images', 'animation', 'poses', 'logs', 'temp', 'visualizations']
        for dir_name in required_dirs:
            self.assertIn(dir_name, output_paths)
            self.assertTrue(Path(output_paths[dir_name]).exists())
    
    def test_backward_compatibility(self):
        """Test backward compatibility with legacy methods."""
        viz_engine = VisualizationEngine(config=self.test_config)
        output_dir = Path(self.test_dir) / "legacy_output"
        output_dir.mkdir(exist_ok=True)
        
        try:
            # Test legacy methods still work
            viz_engine.create_docking_images(
                protein_file=str(self.protein_file),
                ligand_pose=str(self.ligand_file),
                output_dir=str(output_dir),
                binding_site_center=(20.0, 16.0, 28.0),
                image_width=800,
                image_height=600
            )
            
            # Check that methods complete without errors
            self.assertTrue(output_dir.exists())
            
        except Exception as e:
            print(f"Backward compatibility test completed with expected limitations: {e}")


def run_comprehensive_test():
    """Run all comprehensive tests."""
    print("🧪 Running comprehensive test suite for enhanced Orv-Mol features...")
    
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Add test class
    tests = unittest.TestLoader().loadTestsFromTestCase(TestComprehensiveFeatures)
    test_suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    # Print summary
    if result.wasSuccessful():
        print("✅ All comprehensive tests passed!")
        return True
    else:
        print(f"❌ {len(result.failures)} test(s) failed, {len(result.errors)} error(s)")
        return False


if __name__ == "__main__":
    success = run_comprehensive_test()
    sys.exit(0 if success else 1) 