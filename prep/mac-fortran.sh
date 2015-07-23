#! /bin/sh

echo "downloading gfortran-4.8.2 and unpacking it"

mkdir -p ~/gfortran
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
tar -C ~/gfortran -jxvf gfortran-4.8.2-darwin13.tar.bz2

echo "cleaning up: removing tarball"
rm -f gfortran-4.8.2-darwin13.tar.bz2

echo "updating Makevars file in ~/.R to tell R about the new fortran"
mkdir -p ~/.R
cat >> ~/.R/Makevars <<EOF
F77 = ~/gfortran/usr/local/bin/gfortran
FC = ~/gfortran/usr/local/bin/gfortran
FLIBS = -L~/gfortran/usr/local/lib
EOF

echo "now re-run the package installation and pomp test scripts"
