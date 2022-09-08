from setuptools import setup, find_packages

ENTRY_POINTS = {
    'console_scripts': [
        'parpe_petab_to_hdf5 = parpe.hdf5_pe_input:main',
    ]
}

setup(
    name='parpe',
    version='0.0.0',
    url='https://github.com/ICB-DCM/parPE/',
    author='Daniel Weindl',
    author_email='daniel.weindl@helmholtz-muenchen.de',
    description='parpe python package',
    packages=find_packages(),
    install_requires=['numpy>=1.18.1',
                      'termcolor>=1.1.0',
                      'colorama>=0.4.3',
                      'petab>=0.1.18',
                      'amici>=0.11.15',
                      'h5py>=3.0.0',
                      'python-libsbml>=5.17.0',
                      'snakemake>=5.10.0',
                      'coloredlogs>=15.0',
                      'scipy',
                      ],
    entry_points=ENTRY_POINTS,
)
