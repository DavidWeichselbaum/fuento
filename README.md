# fuento
the **fu**nctional **en**richment **to**ol for bioinformatic analysis.
Read the *Manual.pdf* for detailed information.
This is part of the paper "Fuento: functional enrichment for bioinformatics" in *Bioinformatics*.

## REQUIREMENTS:
Requires curl and boost.
Install dependencies with: 
> apt-get install libcurl4-gnutls-dev libboost-all-dev

## INSTALL ON LINUX:
Install Fuento with:
> make
> sudo make install
uninstall with:
> sudo make uninstall

## INSTALL ON OSX:
To run fuento on Mac OS you should have CURL and BOOST library installed. 
If you don't have them, this can be simply done using MacPort - an open-source system for compiling, installing and upgrading.
Please, follow the instruction below.
1. Install Xcode (Apple developer tools), if you don't have it
2. Install MacPort that fits to your Mac OS version (https://www.macports.org/install.php)
3. Install BOOST libraries: 
 > sudo port install boost +universal
4. Install CURL libraries: 
 > sudo port install curl
The provided binary file of fuento is complied for OS 10.9.5

Fuento on Mac can be compiled as following:
> clang++  fuento.cpp -I [headers path]  -L[library path] -lcurl -lboost_regex-mt -lboost_filesystem-mt -std=c++0x --stdlib=libc++ -O2 -w -o fuento
default: [headers path] = /opt/local/include/ [library path]= /opt/local/lib/
Copy the binary to /usr/bin/fuento

## CREATE BACKGROUNDS:
Backgrounds are created from whitespace-separated file of UniProt IDs:
> fuento -c FILE
List backgrounds and gene ontolog slims with:
> fuento -l
Update database and backgrounds with:
> fuento -U

## ANALYSIS:
Start vanila analysis with:
> fuento BKG FILE(s)
BKG can be filepath to background file or corresponding number in 'fuento -l'.
FILE(s) are whitespace-separated files of UniProt IDs

Different usage cenarios and workflows are detailed in *Manual.pdf*

You received example data with this software, comprised of files containing coding sequences and backgrounds for a set of human
proteins, a file with prediction of the fraction of unfolded residues as calculated by IUpred for the same genes and a folder 
with subsets of these genes, divided by predicted fraction of unfolded residues. The example data corresponds to the analysis
done in the publication regarding this software.

David T. Weichselbaum (david.t.weichselbaum@gmail.com) 28.02.2017
