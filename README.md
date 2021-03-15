C++ implementation for "Identifying locations susceptible to micro-anatomical reentry using a spatial network representation of atrial fibre maps - a proof of concept"
==============================================

James Coleman, Samuel Dobson

Contact JamesAlecColeman@gmail.com

Implementation of reading/writing .npy files in C++ ("numpy.hpp") is owned by Masaki Saito

Environment used
--------
C++14

gcc 4.8.5

File information
---------------------
/tractography/ contains the program required to construct fibres, given a flat fibre orientation dataset in .npy in the form (i, j, k, v1, v2, v3)

/null/ contains the program to randomly place nodes >dSep apart over the fibre orientation dataset, without the construction of fibres

/DDM/ contains the program that completes the networks generated using tractography/null and runs a discrete diffusion model on these networks. The amount of re-entries at each node are saved as .npy.

    
