from setuptools import setup, find_packages

setup(
    name="barcode-extractor",
    version="1.0.0",
    description="Efficient barcode extraction for single-cell omics",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[],
    extras_require={
        'zstandard': ['zstandard>=0.14.0'],
    },
    entry_points={
        'console_scripts': [
            'extract-barcodes=barcode_extractor.cli:main',
        ],
    },
)
