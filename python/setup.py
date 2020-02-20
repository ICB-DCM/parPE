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
    install_requires=['numpy',
                      'termcolor',
                      'colorama',
                      'petab>=0.1.2',
                      'amici>=0.10.18',
                      'h5py',
                      'python-libsbml>=5.17.0',
                      'jinja2',
                      'snakemake'
                      ],
    entry_points=ENTRY_POINTS,
)
