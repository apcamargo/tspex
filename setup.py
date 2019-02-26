from setuptools import setup, find_packages

setup(
    name='tspex',
    version='0.1.0',
    packages=find_packages(include=['tspex']),
    license='GNU General Public License v3.0',
    description='A Python package for calculating tissue-specificity metrics for gene expression',
    long_description=open('README.md').read(),
    install_requires=['numpy', 'pandas'],
    url='https://github.com/apcamargo/tspex',
    author='Antonio Pedro Camargo',
    author_email='antoniop.camargo@gmail.com'
)