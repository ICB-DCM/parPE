#!/bin/bash
tar -xzf cblas.tgz
cd CBLAS
ln -sf Makefile.LINUX Makefile.in
make -j12 alllib

