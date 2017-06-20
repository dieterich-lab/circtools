#!/bin/bash

cd /tmp/
git clone https://github.com/dieterich-lab/DCC.git
cd DCC
python setup.py install --user

cd ..
git clone https://github.com/dieterich-lab/FUCHS.git
cd FUCHS
python setup.py install --user

rm /tmp/FUCHS/ -rf
rm /tmp/DCC/ -rf
