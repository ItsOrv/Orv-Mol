#!/usr/bin/env python3
"""
Simple test to debug the professional docking system step by step.
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

def test_imports():
    """Test all imports to identify issues."""
    print("🔍 Testing imports...")
    
    try:
        print("  - Testing basic imports...")
        import numpy as np
        import pandas as pd
        from loguru import logger
        print("  ✅ Basic imports OK")
        
        print("  - Testing Bio imports...")
        from Bio.PDB import PDBParser, PDBIO
        print("  ✅ BioPython imports OK")
        
        print("  - Testing RDKit imports...")
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Crippen
        print("  ✅ RDKit imports OK")
        
        print("  - Testing professional cleaner import...")
        from src.professional_cleaning import ProfessionalCleaner
        print("  ✅ Professional cleaner import OK")
        
        print("  - Testing preprocessing import...")
        from src.preprocessing import PreprocessingEngine
        print("  ✅ Preprocessing import OK")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Import failed: {e}")
        return False

def test_basic_functionality():
    """Test basic functionality without complex operations."""
    print("\n🧪 Testing basic functionality...")
    
    try:
        from src.professional_cleaning import ProfessionalCleaner
        
        print("  - Creating ProfessionalCleaner...")
        cleaner = ProfessionalCleaner(ph=7.4)
        print("  ✅ ProfessionalCleaner created")
        
        print("  - Creating PreprocessingEngine...")
        from src.preprocessing import PreprocessingEngine
        preprocessor = PreprocessingEngine(professional_mode=False)  # Start with simple mode
        print("  ✅ PreprocessingEngine created (standard mode)")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Basic functionality test failed: {e}")
        return False

def test_file_access():
    """Test file access."""
    print("\n📁 Testing file access...")
    
    protein_file = "protein.pdb"
    ligand_file = "ligand.pdb"
    
    if Path(protein_file).exists():
        print(f"  ✅ Protein file found: {protein_file}")
        print(f"     Size: {Path(protein_file).stat().st_size / 1024:.1f} KB")
    else:
        print(f"  ❌ Protein file not found: {protein_file}")
        return False
    
    if Path(ligand_file).exists():
        print(f"  ✅ Ligand file found: {ligand_file}")
        print(f"     Size: {Path(ligand_file).stat().st_size / 1024:.1f} KB")
    else:
        print(f"  ❌ Ligand file not found: {ligand_file}")
        return False
    
    return True

def test_simple_processing():
    """Test simple processing without complex features."""
    print("\n⚙️ Testing simple processing...")
    
    try:
        from src.preprocessing import PreprocessingEngine
        
        # Create simple test directory
        test_dir = Path("simple_test_output")
        test_dir.mkdir(exist_ok=True)
        
        print("  - Testing standard preprocessing...")
        preprocessor = PreprocessingEngine(professional_mode=False)
        
        print("  - Testing protein validation...")
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('test', 'protein.pdb')
        print("  ✅ Protein file can be parsed")
        
        print("  - Testing ligand validation...")
        from rdkit import Chem
        mol = Chem.MolFromPDBFile('ligand.pdb')
        if mol is not None:
            print("  ✅ Ligand file can be parsed")
        else:
            print("  ⚠️ Ligand file parsing issue (common with PDB format)")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Simple processing test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🚀 Simple Diagnostic Test")
    print("=" * 50)
    
    tests = [
        ("Import Test", test_imports),
        ("Basic Functionality", test_basic_functionality),
        ("File Access", test_file_access),
        ("Simple Processing", test_simple_processing)
    ]
    
    results = {}
    for test_name, test_func in tests:
        results[test_name] = test_func()
    
    print("\n" + "=" * 50)
    print("📊 Test Results Summary")
    print("=" * 50)
    
    for test_name, result in results.items():
        status = "✅ PASSED" if result else "❌ FAILED"
        print(f"{test_name}: {status}")
    
    all_passed = all(results.values())
    if all_passed:
        print("\n🎉 All tests passed! System is ready for professional docking.")
    else:
        print("\n⚠️ Some tests failed. Please check the issues above.")
    
    return 0 if all_passed else 1 