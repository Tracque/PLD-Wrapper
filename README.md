# PLD-Wrapper
A high-level wrapper, intended to make the package PLD.jl, more accessible and faster.

Due to technical difficulties, PLD-Wrapper is not being released as a Julia program for the foreseeable future.

NB: Currently, the program is unstable/unreliable, as I am still attempting to implement proper case handling for parralelisation.

There are thus, 2 versions currently avaialable:

PLDManagerV1.py, and the corresponding PLDJobV1.jl

PLDManagerV2.py, and the corresponding PLDJobV2.jl

V1 is currently the only functional version, but is slow as it does not attempt to implement parallel computation. It also has issues properly terminating (it will calculate the final face multiple times due to a failure to properly break out of recursion)

V2 should in theory now work, but I have yet to test on a difficult example. Expect an update in the next week. I also plan to integrate a basic GUI for V2 in the future. Currently, V2 may be functional, but I have yet to investigate issues with the parallel subprocesses.

In order to use this program you must have Julia installed, as well as the Arpack, Arblib, LinearAlgebra, Printf, GenericSVD, Oscar and HomotopyContinuation packages. You should follow the official installation instructions for these packages.
You do not need to download PLD separately, as it is included in this repo. 

Currently, you must manually edit the variables in the PLDManager file and then run it in your console. You should not have to touch the Julia files.


