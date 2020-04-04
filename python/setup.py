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
                      'petab>=0.1.7',
                      'amici>=0.10.21',
                      'h5py>=2.10.0',
                      'python-libsbml>=5.17.0',
                      'jinja2>=2.11.1',
                      'snakemake>=5.10.0'
                      ],
    entry_points=ENTRY_POINTS,
)
