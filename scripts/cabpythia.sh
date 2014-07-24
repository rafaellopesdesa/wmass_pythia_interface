#!/bin/bash
# This getsthe environment
. /etc/bashrc 
. /usr/products/etc/setups.sh 

# Change to the job's directory on the local CAB node 
cd /scratch/$PBS_JOBID 

# Get a new kerberos ticket 
kbatch 

# Copy stuff over from clued0 
scp /prj_root/3122/wz2_write/rclsa/PYTHIA_SAMPLES/RunIIb3/pythia.tar.gz .
tar -xzf pythia.tar.gz
rm -rf pythia.tar.gz

# Run
cd pythia_wz
source compile.tcsh
sed "s/XXX_RANSEED2/$RANDOM/" < pythia_gam-z_ee_RunIIb3.cards > cards.txt
./a.out < cards.txt
cd ../tupleMaker
./tupleMaker ../pythia_wz/bosons.txt bosons.root
scp bosons.root /prj_root/3122/wz2_write/rclsa/PYTHIA_SAMPLES/output2/wmass_dzero_v010100_pythia_gam-z_ee_sm.n_${job}.root
