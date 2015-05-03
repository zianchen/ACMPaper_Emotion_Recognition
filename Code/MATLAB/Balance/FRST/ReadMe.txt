This folder contains Matlab code for running symmetry detection algorithm presented in

G. Loy and A. Zelinsky, 
Fast Radial Symmetry for Detecting Points of Interest, 
IEEE Transactions on Pattern Analysis and Machine Intelligence, 
August 2003.

Type 'help FRST' at the Matlab command prompt for information on how to use the function.

The function histc_weighted is included. It is a modification of the Matlab function histc and uses a mex file histc_weighted_mex.c. Compiled versions of this function are provided for Mac PPC and Intel systems. If you are using a different OS you will need to compile the mex file before running the code. This is done by typing 'mex histc_weighted_mex.c' at the Matlab command prompt from the FRST directory.
