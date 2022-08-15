#!/bin/bash
#make directory
mkdir -p run_$1
cp run_batch.py run_$1/run_batch.py
cp twoLay.py run_$1/twoLay.py
cp Dmmex_R14B_4.f run_$1/Dmmex_R14B_4.f
cp Dmmex_R14B_4.cpython-39-x86_64-linux-gnu.so run_$1/Dmmex_R14B_4.cpython-39-x86_64-linux-gnu.so
#cp {run_EAP_batch.py, EAP_tn.py, Dmmex_R14B_4.f, Dmmex_R14B_4.cpython-39-x86_64-linux-gnu.so} run_$1/
cd run_$1

#set environment
source /usr/share/modules/init/bash
module use -a /swbuild/analytix/tools/modulefiles
module load miniconda3/v4
source activate tf2_8

#run code
echo "Executing run $1 on" `hostname` "in $PWD"
/swbuild/analytix/tools/miniconda3_220407/envs/tf2_8/bin/python run_batch.py $1