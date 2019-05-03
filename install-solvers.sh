#!/usr/bin/env bash
# Install solvers, for Linux/Ubuntu only
sudo apt-get update && sudo apt-get install -y libboost-dev
wget https://ampl.com/netlib/ampl/solvers.tgz
tar -xf solvers.tgz
( cd solvers && ./configure && make )
( export ASL_BUILD=`pwd`/solvers/sys.x86_64.Linux && make )
wget https://ampl.com/dl/open/ipopt/ipopt-linux64.zip
unzip ipopt-linux64.zip
sudo cp ipopt /usr/local/bin/
