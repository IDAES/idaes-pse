#!/usr/bin/env bash
# Installation of solvers for development
# and CircleCI environments
CWD=$(pwd)

sudo apt-get update
echo "America/Los_Angeles" > /etc/timezone

# Expand idaes-coinbinary contents:
mkdir -p /tmp/idaes-coinbinary && cd /tmp/idaes-coinbinary
wget https://idaes-files.s3.amazonaws.com/public/idaes-coinbinary-1.2.0.dev0.zip
unzip idaes-coinbinary-1.2.0.dev0.zip -d /usr/local/coinor-optimization-suite-1.8
cp /usr/local/coinor-optimization-suite-1.8/bin/* /usr/local/bin

cd $CWD
sudo apt-get install -y libboost-dev libgfortran3

wget https://ampl.com/netlib/ampl/solvers.tgz
tar -xf solvers.tgz
cd solvers
./configure && make

# export ASL_BUILD=/home/idaes/solvers/sys.x86_64.Linux
