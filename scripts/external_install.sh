#!/bin/bash

# install primer-design R part
Rscript R/install.R

# install DCC
cd /tmp/
git clone https://github.com/dieterich-lab/DCC.git
cd DCC
python setup.py install --user

# install FUCHS
cd ..
git clone https://github.com/dieterich-lab/FUCHS.git
cd FUCHS
python setup.py install --user

# remove all temporary files
rm /tmp/FUCHS/ -rf
rm /tmp/DCC/ -rf

