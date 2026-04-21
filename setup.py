from setuptools import setup, find_packages

setup(
    name="barcode-corrector",
    version="1.0.0",
    description="Efficient barcode correction for single-cell genomics",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[],
    extras_require={
        'zstandard': ['zstandard>=0.14.0'],
    },
    entry_points={
        'console_scripts': [
            'correct-barcodes=barcode_corrector.cli:main',
        ],
    },
)
