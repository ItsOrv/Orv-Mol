#!/usr/bin/env python3
"""
Setup script for Automated Molecular Docking Pipeline
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="orv-mol-dock",
    version="1.0.0",
    author="Molecular Modeling Agent",
    description="Automated molecular docking pipeline using AutoDock Vina",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "rdkit>=2022.03.1",
        "biopython>=1.79",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "pillow>=8.3.0",
        "tqdm>=4.62.0",
        "click>=8.0.0",
        "colorama>=0.4.4",
        "pyyaml>=6.0",
        "loguru>=0.6.0",
    ],
    extras_require={
        "visualization": ["pymol-open-source>=2.5.0"],
        "advanced": ["meeko>=0.4.0", "vina>=1.2.3"],
    },
    entry_points={
        "console_scripts": [
            "dock=dock:main",
        ],
    },
)
