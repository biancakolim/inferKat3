
from setuptools import setup, find_packages

setup(
    name='inferKat2',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'scanpy',
        'scipy',
        'umap-learn',
        'leidenalg',
        'infercnvpy',
        'hmmlearn', 
        'pykalman',
    ],
    author='CEKN',
    description='CSCB CEKN Final Project Package',
)
