#!/bin/bash
wget http://www.feynarts.de/cuba/Cuba-4.2.tar.gz
tar xf Cuba-4.2.tar.gz
cd Cuba-4.2
export CFLAGS=-fPIC
./configure
make
cp libcuba.a $CMSSW_BASE/lib/$SCRAM_ARCH/
cp cuba.h $CMSSW_BASE/src/TTH/MEIntegratorStandalone/interface/
