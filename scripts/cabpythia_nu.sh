#!/bin/bash
# This getsthe environment
. /etc/bashrc 
. /usr/products/etc/setups.sh 

# Change to the job's directory on the local CAB node 
cd /scratch/$PBS_JOBID 

# Get a new kerberos ticket 
kbatch 

# Copy stuff over from clued0 
scp /prj_root/3122/wz2_write/rclsa/PYTHIA_SAMPLES/wmass_pythia_interface/pythia.tar.gz .
tar -xzf pythia.tar.gz
rm -rf pythia.tar.gz

# Run
cd pythia_wz
source compile.tcsh
sed "s/XXX_RANSEED2/$RANDOM/" < pythia_z_nunu_RunIIb3.cards > cards.txt
./a.out < cards.txt
cd ../tupleMaker
./tupleMaker ../pythia_wz/bosons.txt bosons.root
scp bosons.root /prj_root/2673/wz2_write/rclsa/PYTHIA_SAMPLES/output3/wmass_dzero_v010259_pythia_6_409_tuneA_-z_nunu_sm.n_${job}.root
