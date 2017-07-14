#!/bin/bash
wget -O amici.zip https://github.com/ICB-DCM/AMICI/archive/master.zip
unzip amici.zip
cd AMICI-master
scripts/run-build.sh && scripts/run-cpputest.sh
