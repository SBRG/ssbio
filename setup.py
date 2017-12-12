from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
        name='ssbio',
        version='0.9.9.2',
        author='Nathan Mih',
        author_email='nmih@ucsd.edu',
        license='MIT',
        url='http://github.com/SBRG/ssbio',
        download_url = 'https://github.com/SBRG/ssbio/archive/v0.9.9.2.tar.gz',
        description='Tools to enable structural systems biology',
        packages=find_packages(),
        package_dir={'ssbio': 'ssbio'},
        # scripts=['ssbio/protein/structure/utils/cleanpdb.py',
        #          'ssbio/protein/structure/utils/mutatepdb.py',
        #          'ssbio/protein/structure/utils/tleap.py'],
        #          'ssbio/protein/structure/properties/msms.py'],
        long_description=open('README.rst').read(),
        install_requires=required
)
