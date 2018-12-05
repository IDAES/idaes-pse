# Building the IAPWS95 Library

## Download and Build The ASL (AMPL Solver Libraries)

Download The ASL from https://ampl.com/netlib/ampl/solvers.tgz.

Extract the file to a convenient location. In the ASL directory run configure then make.

```sh
./configure
make
```

## Set the ASL_BUILD environment variables

Set the ASL_BUILD environment variable to the location of the ASL build. The ASL build directory should be in the ASL directory with a name that starts with "sys." followed by something related to the system architecture.

For example (__replace with you own directory__):

```sh
export ASL_BUILD=$HOME/local/src/solvers/sys.x86_64.Linux
```
## Get boost

The memoization feature use a hash function from the boost library to hash tuples.  Install boost if necessary. If Boost is installed in a nonstandard location you can use the BOOST_HEADER environment variable to set a location to look for it.  For example:

```sh
export BOOST_HEADER=$HOME/anaconda2/include
```

## Build the IAPWS95 Library

In the idaes/property_models/iapws95 directory run make

```sh
make
```
