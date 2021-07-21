#!/bin/bash 
sudo apt-get install g++ cmake libeigen3-dev libboost-all-dev python3-dev python3-numpy pybind11-dev git doxygen python3-sphinx

git clone https://github.com/openbabel/openbabel.git
mkdir openbabel-build
cd openbabel-build
cmake ../openbabel
sudo make install
cd ..

git clone https://github.com/yesint/pteros.git pteros

mkdir pteros_build
mkdir pteros_install
cd pteros_build

cmake ../pteros -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../pteros_install -DWITH_GROMACS=OFF
sudo make install

cd ..
echo "source" $(pwd)"/pteros_install/lib/pterosrc" >> ~/.bashrc
source ~/.bashrc

rm -rf pteros
rm -rf pteros_build
rm -rf openbabel
rm -rf openbabel-build
