# 🚀 Quick Start Guide

Get up and running with Orv-Mol in 5 minutes!

## ⚡ 1-Minute Installation

```bash
# Clone the repository
git clone <repository-url>
cd Orv-Mol

# Run the installer
chmod +x install.sh
./install.sh

# Activate virtual environment
source venv/bin/activate
```

## 🧪 First Docking Run

```bash
# Basic docking with professional cleaning
python dock.py data/input/protein.pdb data/input/ligand_small.pdb

# With custom parameters
python dock.py data/input/protein.pdb data/input/ligand_small.pdb \
               --exhaustiveness 32 --professional-cleaning

# Quick example script
python examples/quick_start.py
```

## 📊 Expected Results

- **Binding Affinity**: -4.7 to -5.9 kcal/mol (for test data)
- **Processing Time**: 2-5 minutes
- **Output**: Images, poses, logs, and analysis reports

## 🔧 Common Commands

| Task | Command |
|------|---------|
| **Standard Docking** | `python dock.py protein.pdb ligand.pdb` |
| **Professional Mode** | `python dock.py protein.pdb ligand.pdb --professional-cleaning` |
| **High Accuracy** | `python dock.py protein.pdb ligand.pdb --exhaustiveness 32` |
| **Custom Center** | `python dock.py protein.pdb ligand.pdb --center-x 168.67 --center-y 29.82 --center-z 182.29` |
| **No Visualization** | `python dock.py protein.pdb ligand.pdb --skip-visualization` |

## 📁 Output Structure

After successful docking:

```
results/
├── images/           # Publication-ready visualizations
├── poses/           # Docked ligand conformations
├── logs/            # Detailed execution logs
└── temp/            # Temporary processing files
```

## 🔍 Understanding Results

- **Binding Affinity**: Lower = stronger binding (good: < -7 kcal/mol)
- **Multiple Poses**: Different binding orientations
- **Quality Scores**: Protein (structure) and ligand (drug-likeness) quality

## 🧼 Professional Features

The professional cleaning system provides:

- **Protein Quality Assessment** (1-10 scale)
- **Drug-likeness Scoring** (ADMET properties)
- **Automatic Parameter Optimization**
- **Multi-conformer Generation**

## ⚠️ Troubleshooting

| Problem | Solution |
|---------|----------|
| `vina command not found` | Run `sudo apt install autodock-vina` |
| `Import error: rdkit` | Run `pip install rdkit` |
| `No binding affinity` | Check input file quality |
| `PDBQT format error` | Use professional cleaning mode |

## 📚 Next Steps

1. **Read the full README.md** for detailed documentation
2. **Explore examples/** directory for more use cases
3. **Check tests/** for validation scripts
4. **Modify config.yaml** for default settings

## 🆘 Getting Help

- **Issues**: Open a GitHub issue
- **Questions**: Check documentation
- **Examples**: Run `python examples/quick_start.py`

---

**Ready to dock? Run your first example now! 🧪** 