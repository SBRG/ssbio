## EMBOSS programs (needle and others)
### Quick links:

- Home page: http://emboss.sourceforge.net/

### Installation for Ubuntu/Linux

1. Install the EMBOSS package which contains many programs

        sudo apt-get install emboss

2. And then once that installs, try running the "needle" program:

        needle

### Installation for Mac OSX and other Unix systems

1. Just install from source

       ./configure
       make
       sudo make install

2. On Unix you may encounter:

       /usr/local/bin/embossupdate: error while loading shared libraries: libnucleus.so.6: cannot open shared object file: No such file or directory

Just run (according to https://www.biostars.org/p/86327/):

       sudo ldconfig

And try running an EMBOSS program, such as `needle`.