from setuptools import setup

setup(
    name='musr_espresso',
    version='0.1',
    packages=['musr_espresso'],
    url='',
    license='',
    author='John Wilkinson',
    install_requires=["numpy", "ase", "spglib", "soprano", "pandas"],
    author_email='john.wilkinson@physics.ox.ac.uk',
    description='Collection of scripts for DFT+mu'
)
