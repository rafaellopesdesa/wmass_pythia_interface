                                           Jan Stark, February 24th, 2007

This directory contains a small stand-alone program that can be used to
run Pythia to generate W -> e nu and Z -> e e events, and write them
out in a format that looks like the ASCII files from ResBos+Photos. 
Cardfiles that reproduce the standard full Monte Carlo samples from the 
D0 farms are also included. The main purpose of this tool is to generate 
fast MC templates for the m(W) Monte Carlo closure test.

To build the executable, simply type:

 source compile.tcsh

This script is written for the tcsh, if you use another shell you might 
have to adapt it. No setups are required prior to execution of this
script. This script does setup the D0RunII product version p17.10.00.
This is done in order to get access to the following packages:

 - the D0 version of GCC 3.4.3,

 - Pythia 6.323,

 - Les Houches Accord PDF 4.1.1 .

After compilation, you can run the program to produce Z -> e e events:

 ./a.out < pythia_gam-z_ee_sm.n.cards

The generated events will be written to a file called bosons.txt .
Use the cardfile pythia_w_enu_sm.n.cards if you would like to produce
W -> e nu events.

Once you have built the executable, it is no longer necessary to setup
D0RunII. In a new window, just say:

 source setup.tcsh

and you can run the binary. The script setup.tcsh only setups the LHAPDF
product. We do need this one to have access to the files with the parton
density functions.


