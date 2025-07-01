#!/usr/bin/env python3
"""
Syntax and structure test for Orv-Mol project.
Tests code syntax, imports structure, and basic file organization.
"""

import sys
import os
import ast
import importlib.util
from pathlib import Path

def test_python_syntax(file_path):
    """Test if a Python file has valid syntax."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            source = f.read()
        
        # Parse the AST to check syntax
        ast.parse(source)
        return True, None
    except SyntaxError as e:
        return False, f"Syntax error: {e}"
    except Exception as e:
        return False, f"Error reading file: {e}"

def test_file_structure():
    """Test project file structure."""
    print("📁 Testing project structure...")
    
    required_files = [
        'dock.py',
        'config.yaml',
        'requirements.txt',
        'README.md',
        'src/__init__.py',
        'src/docking.py',
        'src/preprocessing.py',
        'src/output_parser.py',
        'src/visualization.py',
        'src/utils.py',
        'src/professional_cleaning.py'
    ]
    
    missing_files = []
    for file_path in required_files:
        if not Path(file_path).exists():
            missing_files.append(file_path)
        else:
            print(f"  ✅ {file_path}")
    
    if missing_files:
        print(f"  ❌ Missing files: {missing_files}")
        return False
    
    print("  🎉 All required files present")
    return True

def test_python_files_syntax():
    """Test syntax of all Python files."""
    print("\n🐍 Testing Python file syntax...")
    
    python_files = []
    for pattern in ['*.py', 'src/*.py', 'tests/*.py', 'scripts/*.py', 'examples/*.py']:
        python_files.extend(Path('.').glob(pattern))
    
    failed_files = []
    for py_file in python_files:
        if py_file.name.startswith('.'):
            continue
            
        success, error = test_python_syntax(py_file)
        if success:
            print(f"  ✅ {py_file}")
        else:
            print(f"  ❌ {py_file}: {error}")
            failed_files.append((py_file, error))
    
    if failed_files:
        print(f"\n  💥 Syntax errors in {len(failed_files)} files")
        return False
    
    print(f"  🎉 All {len(python_files)} Python files have valid syntax")
    return True

def test_import_structure():
    """Test import structure without executing imports."""
    print("\n📦 Testing import structure...")
    
    src_files = {
        'src/utils.py': ['pathlib', 'yaml'],
        'src/visualization.py': ['pathlib', 'typing'],
        'src/output_parser.py': ['pathlib', 'typing'],
        'src/docking.py': ['pathlib', 'subprocess'],
        'src/preprocessing.py': ['pathlib', 'subprocess']
    }
    
    for file_path, expected_imports in src_files.items():
        if not Path(file_path).exists():
            continue
            
        try:
            with open(file_path, 'r') as f:
                source = f.read()
            
            tree = ast.parse(source)
            imports = []
            
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        imports.append(alias.name)
                elif isinstance(node, ast.ImportFrom):
                    if node.module:
                        imports.append(node.module.split('.')[0])
            
            print(f"  ✅ {file_path}: imports detected")
            
        except Exception as e:
            print(f"  ❌ {file_path}: error analyzing imports - {e}")
            return False
    
    print("  🎉 Import structure analysis completed")
    return True

def test_config_files():
    """Test configuration files."""
    print("\n⚙️ Testing configuration files...")
    
    # Test YAML config
    if Path('config.yaml').exists():
        try:
            import yaml
            with open('config.yaml', 'r') as f:
                config = yaml.safe_load(f)
            print("  ✅ config.yaml: valid YAML")
        except ImportError:
            print("  ⚠️ config.yaml: PyYAML not available, skipping validation")
        except Exception as e:
            print(f"  ❌ config.yaml: {e}")
            return False
    
    # Test requirements.txt
    if Path('requirements.txt').exists():
        try:
            with open('requirements.txt', 'r') as f:
                lines = f.readlines()
            
            valid_lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
            print(f"  ✅ requirements.txt: {len(valid_lines)} dependencies listed")
        except Exception as e:
            print(f"  ❌ requirements.txt: {e}")
            return False
    
    print("  🎉 Configuration files validated")
    return True

def test_code_organization():
    """Test code organization and class structure."""
    print("\n🏗️ Testing code organization...")
    
    main_classes = {
        'src/docking.py': ['DockingEngine'],
        'src/preprocessing.py': ['PreprocessingEngine'],
        'src/output_parser.py': ['OutputParser'],
        'src/visualization.py': ['VisualizationEngine'],
        'src/professional_cleaning.py': ['ProfessionalCleaner']
    }
    
    for file_path, expected_classes in main_classes.items():
        if not Path(file_path).exists():
            continue
            
        try:
            with open(file_path, 'r') as f:
                source = f.read()
            
            tree = ast.parse(source)
            found_classes = []
            
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    found_classes.append(node.name)
            
            for expected_class in expected_classes:
                if expected_class in found_classes:
                    print(f"  ✅ {file_path}: {expected_class} class found")
                else:
                    print(f"  ❌ {file_path}: {expected_class} class missing")
                    return False
                    
        except Exception as e:
            print(f"  ❌ {file_path}: error analyzing classes - {e}")
            return False
    
    print("  🎉 Code organization validated")
    return True

def test_documentation():
    """Test documentation presence."""
    print("\n📚 Testing documentation...")
    
    doc_files = ['README.md', 'QUICK_START.md']
    missing_docs = []
    
    for doc_file in doc_files:
        if Path(doc_file).exists():
            size = Path(doc_file).stat().st_size
            print(f"  ✅ {doc_file}: {size} bytes")
        else:
            missing_docs.append(doc_file)
    
    if missing_docs:
        print(f"  ⚠️ Missing documentation: {missing_docs}")
    
    print("  🎉 Documentation check completed")
    return True

def run_comprehensive_syntax_test():
    """Run all syntax and structure tests."""
    print("🧪 Orv-Mol Syntax & Structure Test Suite")
    print("=" * 50)
    
    tests = [
        ("File Structure", test_file_structure),
        ("Python Syntax", test_python_files_syntax),
        ("Import Structure", test_import_structure),
        ("Configuration Files", test_config_files),
        ("Code Organization", test_code_organization),
        ("Documentation", test_documentation)
    ]
    
    results = {}
    for test_name, test_func in tests:
        try:
            results[test_name] = test_func()
        except Exception as e:
            print(f"  💥 {test_name} test crashed: {e}")
            results[test_name] = False
    
    # Print summary
    print("\n" + "=" * 50)
    print("📊 Test Results Summary")
    print("=" * 50)
    
    passed = 0
    for test_name, result in results.items():
        status = "✅ PASSED" if result else "❌ FAILED"
        print(f"{test_name:.<25} {status}")
        if result:
            passed += 1
    
    print(f"\nResults: {passed}/{len(tests)} tests passed")
    
    if passed == len(tests):
        print("\n🎉 All syntax and structure tests passed!")
        print("✨ Project structure is well-organized and syntactically correct")
        return True
    else:
        print(f"\n⚠️ {len(tests) - passed} test(s) failed")
        print("🔧 Please fix the issues above before proceeding")
        return False

if __name__ == "__main__":
    success = run_comprehensive_syntax_test()
    sys.exit(0 if success else 1) 