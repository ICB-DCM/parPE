from setuptools import setup, find_packages

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
                      'petab',
                      'amici',
                      'h5py',
                      'python-libsbml'],
)
