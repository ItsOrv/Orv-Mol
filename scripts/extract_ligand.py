#!/usr/bin/env python3
"""
Ligand Extraction Tool for Orv-Mol
"""

import os
from loguru import logger

def extract_ligand_from_complex(complex_file, ligand_name, output_file):
    """Extract a specific ligand from a protein-ligand complex"""
    logger.info(f"Extracting ligand '{ligand_name}' from {complex_file}")
    
    ligand_lines = []
    
    try:
        with open(complex_file, 'r') as f:
            for line in f:
                # Extract HETATM records for the specified ligand
                if line.startswith('HETATM') and ligand_name in line:
                    ligand_lines.append(line)
        
        if not ligand_lines:
            logger.error(f"❌ No ligand '{ligand_name}' found in {complex_file}")
            return False
        
        # Write extracted ligand to output file
        with open(output_file, 'w') as f:
            f.write("REMARK   Extracted ligand from complex\n")
            f.write(f"REMARK   Source: {complex_file}\n")
            f.write(f"REMARK   Ligand: {ligand_name}\n")
            for line in ligand_lines:
                f.write(line)
            f.write("END\n")
        
        logger.success(f"✅ Ligand extracted: {output_file} ({len(ligand_lines)} atoms)")
        return True
        
    except Exception as e:
        logger.error(f"❌ Error extracting ligand: {str(e)}")
        return False


def analyze_hetatm_records(complex_file):
    """Analyze HETATM records in a PDB file to identify ligands"""
    logger.info(f"Analyzing HETATM records in {complex_file}")
    
    hetatm_records = {}
    
    try:
        with open(complex_file, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    # Extract residue name (columns 18-20)
                    res_name = line[17:20].strip()
                    if res_name not in hetatm_records:
                        hetatm_records[res_name] = 0
                    hetatm_records[res_name] += 1
        
        if hetatm_records:
            logger.info("HETATM Summary:")
            for res_name, count in sorted(hetatm_records.items()):
                logger.info(f"  {res_name}: {count} atoms")
        else:
            logger.warning("⚠️ No HETATM records found")
        
        return hetatm_records
        
    except Exception as e:
        logger.error(f"❌ Error analyzing file: {str(e)}")
        return {}


if __name__ == "__main__":
    logger.info("Ligand Extraction Tool")
    logger.info("=" * 30)
    
    # Check if ligand.pdb exists to extract from
    complex_file = "data/input/ligand.pdb"
    
    # Analyze HETATM records first
    hetatm_records = analyze_hetatm_records(complex_file)
    
    if hetatm_records:
        # Extract the first ligand found (or you can specify which one)
        ligand_names = list(hetatm_records.keys())
        if ligand_names:
            ligand_name = ligand_names[0]  # Take the first one
            output_file = f"data/input/ligand_small.pdb"
            
            if extract_ligand_from_complex(complex_file, ligand_name, output_file):
                logger.info("\nNow you can test with the small molecule ligand:")
                logger.info(f"python dock.py data/input/protein.pdb {output_file}")
    else:
        if not os.path.exists(complex_file):
            logger.error("❌ ligand.pdb not found")
        else:
            logger.info("No ligands found to extract") 