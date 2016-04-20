from distutils.core import setup

setup(
    name='ssbio',
    version='0.1',
    packages=['ssbio',],
    license='MIT license',
    long_description=open('README.txt').read(), requires=['numpy','biopython']
)