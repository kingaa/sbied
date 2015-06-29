#! /bin/sh

mkdir -p ~/gfortran
cd ~/gfortran
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
tar jxvf gfortran-4.8.2-darwin13.tar.bz2 

mkdir -p ~/.R
cd ~/.R
curl -O http://kinglab.eeb.lsa.umich.edu/SBIED/scripts/Makevars
cd ~
