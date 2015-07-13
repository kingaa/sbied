#! /bin/sh

echo "downloading gfortran-4.8.2 and unpacking it"

mkdir -p ~/gfortran
cd ~/gfortran
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
tar jxvf gfortran-4.8.2-darwin13.tar.bz2 

echo "cleaning up: removing tarball"
rm -f gfortran-4.8.2-darwin13.tar.bz2

echo "setting up a Makevars file in ~/.R to tell R about the new fortran"
mkdir -p ~/.R
cd ~/.R
curl -O http://kinglab.eeb.lsa.umich.edu/SBIED/scripts/Makevars
cd ~

echo "now re-run the package installation and pomp test scripts"
