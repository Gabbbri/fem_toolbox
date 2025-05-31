# setup.py

from setuptools import setup, find_packages

setup(
    name='fem_toolbox',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pyfiglet'
    ],
    author='Gabriele',
    author_email='gabrigirardi1@gmail.com',
    description='A Python toolkit for FEM analysis of 1D and 2D beam structures',
    long_description='Tools for assembling stiffness/mass matrices, applying BCs, postprocessing, and modal analysis.',
    long_description_content_type='text/plain',
    #url='https://github.com/yourusername/fem_tools', 
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
)
