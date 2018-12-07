#!/usr/bin/env bash
# Install solvers, for Linux/Ubuntu only
sudo apt-get update && sudo apt-get install -y libboost-dev
wget https://ampl.com/netlib/ampl/solvers.tgz
tar -xf solvers.tgz
( cd solvers && ./configure && make )
( export ASL_BUILD=`pwd`/solvers/sys.x86_64.Linux && \
	cd idaes/property_models/iapws95 && make )
