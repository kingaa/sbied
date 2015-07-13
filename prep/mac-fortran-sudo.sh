#! /bin/sh

echo "downloading gfortran-4.8.2 and unpacking it"
echo "this requires root priveleges"

curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar jxvf gfortran-4.8.2-darwin13.tar.bz2 -C /

echo "cleaning up: removing tarball"
rm -f gfortran-4.8.2-darwin13.tar.bz2

echo "setting up a Makevars file in ~/.R to tell R about the new fortran"
mkdir -p ~/.R
cd ~/.R
curl -O http://kinglab.eeb.lsa.umich.edu/SBIED/scripts/Makevars.sudo -o ~/.R/Makevars
cd ~

echo "now re-run the package installation and pomp test scripts"
