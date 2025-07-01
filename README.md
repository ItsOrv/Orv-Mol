# Orv-Mol: Professional Molecular Docking System

**High-performance automated molecular docking platform powered by AutoDock Vina with advanced preprocessing and optimization algorithms.**

## Overview

Orv-Mol provides a complete solution for protein-ligand docking with professional-grade molecular cleaning, automated parameter optimization, and comprehensive analysis capabilities. Designed for researchers who need reliable, reproducible docking results with minimal manual intervention.

## Key Capabilities

- **Automated Parameter Optimization**: Intelligent exhaustiveness and grid box sizing based on molecular complexity
- **Professional Molecular Cleaning**: Advanced preprocessing with quality scoring and structure validation  
- **Multi-Format Support**: Handles PDB, SDF, MOL2, and PDBQT formats seamlessly
- **Comprehensive Analysis**: Drug-likeness assessment, ADMET properties, and interaction profiling
- **Batch Processing**: High-throughput screening capabilities with parallel execution
- **Publication-Ready Output**: High-quality visualizations and detailed reports

## Installation

### System Requirements
- Python 3.8+
- AutoDock Vina 1.2.5+
- 4+ GB RAM (8+ GB recommended for large proteins)
- Linux/macOS/Windows

### Quick Setup
```bash
# Clone and setup
git clone https://github.com/ItsOrv/Orv-Mol.git
cd Orv-Mol

# Install dependencies
pip install -r requirements.txt

# Install AutoDock Vina
# Ubuntu/Debian: sudo apt install autodock-vina
# macOS: brew install autodock-vina
# Windows: Download from https://autodock.scripps.edu/downloads/

# Verify installation
python dock.py --version
```

## Usage

### Basic Docking
```bash
# Automatic blind docking with professional cleaning
python dock.py protein.pdb ligand.sdf

# Quick screening mode (fast, lower accuracy)
python dock.py protein.pdb ligand.sdf --quick

# Research-grade docking (high accuracy)
python dock.py protein.pdb ligand.sdf --research-grade
```

### Targeted Docking
```bash
# Specify binding site coordinates
python dock.py protein.pdb ligand.sdf \
    --center-x 25.5 --center-y 10.2 --center-z 15.8 \
    --size-x 20 --size-y 20 --size-z 20
```

### Advanced Options
```bash
# Custom parameters with all outputs
python dock.py protein.pdb ligand.sdf \
    --exhaustiveness 24 \
    --num-modes 20 \
    --professional-cleaning \
    --csv-output --json-output \
    --animation --interaction-plots

# Batch processing mode
python dock.py protein.pdb ligand.sdf \
    --batch --parallel-conformers \
    --generate-conformers 5
```

## Performance Modes

| Mode | Exhaustiveness | Speed | Accuracy | Use Case |
|------|---------------|-------|----------|----------|
| `--quick` | 4 | Very Fast | Good | Initial screening |
| Default | 8 | Fast | High | Standard research |
| `--accurate` | 16 | Moderate | Very High | Important targets |
| `--research-grade` | 32 | Slow | Maximum | Publication work |

## Configuration

Create `config.yaml` for custom defaults:

```yaml
# Performance settings
docking:
  exhaustiveness: 16
  num_modes: 15
  energy_range: 3.5

# Quality settings  
preprocessing:
  professional_mode: true
  ph: 7.4
  remove_waters: true
  add_hydrogens: true

# Output preferences
output:
  save_logs: true
  generate_reports: true
  image_format: "png"
  image_dpi: 300
```

## Python API

### Simple Docking Workflow
```python
from src.preprocessing import PreprocessingEngine
from src.docking import DockingEngine

# Initialize with professional cleaning
processor = PreprocessingEngine(professional_mode=True, ph=7.4)
docker = DockingEngine()

# Prepare molecules
protein_pdbqt = processor.prepare_receptor("protein.pdb", "temp/")
ligand_pdbqt = processor.prepare_ligand("ligand.sdf", "temp/")

# Auto-calculate binding box and dock
box_params = processor.calculate_blind_docking_box(protein_pdbqt)
results = docker.run_vina_docking(
    receptor=protein_pdbqt,
    ligand=ligand_pdbqt, 
    box_params=box_params,
    output_dir="results/",
    exhaustiveness=16
)

print(f"Best affinity: {results['best_affinity']} kcal/mol")
```

### Advanced Analysis
```python
from src.professional_cleaning import ProfessionalCleaner
from src.output_parser import OutputParser

# Molecular analysis
cleaner = ProfessionalCleaner()
protein_quality = cleaner.analyze_protein_structure("protein.pdb")
ligand_props = cleaner.calculate_ligand_properties("ligand.sdf")

# Results processing
parser = OutputParser()
detailed_results = parser.parse_vina_output("output.log", "poses.pdbqt", "output/")
binding_analysis = parser.analyze_binding_interactions(detailed_results)
```

## Output Files

### Standard Results
- `docking_summary.txt` - Comprehensive text report
- `docking_results.csv` - Tabular data for analysis
- `poses/` - Individual binding pose structures
- `images/` - Molecular visualizations

### Advanced Results (with flags)
- `docking_results.json` - Structured data export
- `interaction_analysis.txt` - Detailed interaction profiles  
- `animation/` - 3D rotation videos
- `quality_report.txt` - Molecular quality metrics

## Validation & Testing

```bash
# Quick system check
python tests/test_installation.py

# Validate with sample data
python tests/simple_test.py

# Full test suite
python -m pytest tests/ -v
```

## Performance Optimization

### For Large Proteins (>500 residues)
```bash
python dock.py large_protein.pdb ligand.sdf \
    --exhaustiveness 12 \
    --size-x 25 --size-y 25 --size-z 25 \
    --batch
```

### For High-Throughput Screening
```bash
# Process multiple ligands efficiently
for ligand in ligands/*.sdf; do
    python dock.py protein.pdb "$ligand" --quick --skip-visualization
done
```

### Memory Management
- Large proteins: Use `--batch` mode
- Multiple conformers: Enable `--parallel-conformers`
- Low memory: Use `--quick` with `--skip-visualization`

## Troubleshooting

**Common Issues:**

| Problem | Solution |
|---------|----------|
| "Vina not found" | Install AutoDock Vina and add to PATH |
| Memory errors | Reduce exhaustiveness or use `--batch` |
| Slow performance | Use `--quick` mode or reduce grid box size |
| Import errors | Run `pip install -r requirements.txt` |
| Poor results | Enable `--professional-cleaning` |

**Performance Tips:**
- Use targeted docking when binding site is known
- Start with `--quick` for initial screening
- Enable professional cleaning for final results
- Adjust exhaustiveness based on accuracy needs


## Support

- **Documentation**: See `examples/` for working tutorials
- **Issues**: GitHub Issues for bug reports and feature requests
- **API Reference**: Full documentation in `docs/`

---

**License**: MIT | **Requirements**: Python 3.8+, AutoDock Vina 1.2.5+ | **Platform**: Linux, macOS, Windows
