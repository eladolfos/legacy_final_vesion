#!/bin/bash
lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"
lsetup "cmake 3.14.0"
export LHAPDFSYS=/home/yfu/LHAPDF
export PATH=/home/yfu/LHAPDF/bin:${PATH}
export LD_LIBRARY_PATH=/home/yfu/LHAPDF/lib:${LD_LIBRARY_PATH}
