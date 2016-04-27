from setuptools import setup, find_packages

setup(
    name='ssbio',
    version='0.1',
    author='Nathan Mih',
    author_email='nmih@ucsd.edu',
    packages=find_packages(),
    package_dir={'ssbio': 'ssbio'},
    package_data={'ssbio': ['ssbio/etc/*']},
    license='MIT license',
    scripts = ['ssbio/tools/cleanpdb', 'ssbio/tools/mutatepdb.py', 'ssbio/tools/tleap.py', 'ssbio/databases/drugbank.py'],
    long_description=open('README.md').read(),
    install_requires=['biopython','numpy','tqdm','pandas','requests']
)