from setuptools import setup, find_packages

setup(
    name='tspex',
    version='0.1.0',
    packages=find_packages(include=['tspex']),
    license='GNU General Public License v3.0',
    description='A Python package for calculating tissue-specificity metrics for gene expression.',
    long_description=open('README.md').read(),
    install_requires=['numpy', 'pandas >= 0.23'],
    url='https://github.com/apcamargo/tspex',
    keywords=['gene expression', 'tissue-specificity', 'transcriptomics'],
    author='Antonio Pedro Camargo',
    author_email='antoniop.camargo@gmail.com',
    classifiers=[
        'Development Status :: ',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Programming Language :: Python :: 3',
    ]
)