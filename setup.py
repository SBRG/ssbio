from setuptools import setup, find_packages

setup(
    name='ssbio',
    version='0.1',
    author='Nathan Mih',
    author_email='nmih@ucsd.edu',
    license='MIT',
    url='http://github.com/nmih/ssbio',
    description='Various tools and functions to enable structural systems biology',
    packages=find_packages(),
    package_dir={'ssbio': 'ssbio'},
    package_data={'ssbio': ['ssbio/etc/*']},
    scripts = ['ssbio/scripts/cleanpdb', 'ssbio/scripts/aggprop', 'ssbio/scripts/thermostability', 'ssbio/structure/mutatepdb.py', 'ssbio/structure/tleap.py',
               'ssbio/databases/drugbank.py', 'ssbio/structure/properties/msmsprops.py',
               'ssbio/dock/dock.py'],
    long_description=open('README.md').read(),
    install_requires=['biopython',
                      'numpy',
                      'tqdm',
                      'pandas',
                      'requests',
                      'cachetools',
                      'bioservices',
                      'prody',
                      'xmltodict']
)