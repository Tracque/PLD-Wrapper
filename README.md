# PLD-Wrapper.py

A high-level wrapper, intended to make the package PLD.jl, more accessible and faster.


There are 2 versions of this program available at this time, labelled V1 and V2 respectively.


V1 is NOT recommended to use, however if your machine is struggling with V2 and you cannot get access to a machine that can support the parallel processes, it remains an option. At such a time as I release the first "stable" version of this program, V1 will become unavailable unless one searches for old commits to this repo.

  

V2 is the version I am currently recommending for use. It is, from my testing, fully functional and compiles all output into 2 files, for the user's convenience. V2 will soon be supported by a simple GUI, intended to give users some guidance as to the required format for their inputs and enable slightly more verbose output.

When V2 is supported by the GUI, there will remain a version similar to its current form available for use in headless instances, which will hopefully make use of PLD-Wrapper for large calculations easier. (I assume such calculations will mostly be performed by SSH'ing into your institution's systems, so the GUI will just get in the way there)

## Installation Instructions

  
In order to use this program you must have Julia installed, as well as the Arpack, Arblib, LinearAlgebra, Printf, GenericSVD, Oscar and HomotopyContinuation packages. If you have any problems, you should consult the official installation instructions for these packages. However, here's the quick rundown:

Step 1: You MUST be running on a UNIX OS. This means Windows users should either parition their drive so that they can dual boot a Linux system, or make use of the Windows Subsystem for Linux (WSL). If you are performing a fresh install of a Linux system, make sure you run:

```
sudo apt update
sudo apt install build-essential
```

Step 2: Install Julia. Do NOT use your system package manager. It will likely fetch an old or broken version of Julia. Instead, you should use "juliaup". You may install juliaup (and hence Julia) on Ubuntu machines, by using the command:

```curl -fsSL https://install.julialang.org | sh ```

Step 3: Install the required libraries. This part should be simple. Start a Julia instance using the command ```julia``` and then you can simply make use of the native package manager:

```
using Pkg; Pkg.activate()
Pkg.add("Printf")
Pkg.add(...)
...
```

One issue you may encounter is with the package OSCAR failing to precompile. This is due to the package manager fetching bad outdated versions of some dependencies (namely, cxxmake.jl and thus Singular.jl, Polymake.jl and others). For whatever reason, simply running ```Pkg.update()``` will now fetch working versions of the dependencies, solving this issue. (Yes, you can "update" these packages which you added seconds earlier)

Note: PLD.jl is NOT a package that can be added using Pkg. You do not need to worry about downloading PLD separately, as it is included in this repo. (in fact, if you download it yourself, you WILL encounter errors, as OSCAR has undergone some slight syntax changes which required edits to PLD.jl)

Step 4: Make sure you have a working version of python3 installed, as well as the library psutil, which the wrapper uses to avoid frying your machine. You can download psutil by running:

```
pip install psutil
```

Or, 

```
sudo apt-get install python3-psutil
```

------

And with that, you should be up and running! Simply navigate to the directory where you installed PLD-Wrapper and run the command:

```
python3 PLD-Wrapper.py
```

## Usage Instructions

Currently, you must manually edit the variables in the PLDManager file and then run it in your console. You should not have to touch the Julia files.

Once the GUI integration has been completed, if you run the version with the GUI, you will need to give your inputs in the input fields provided.

Note: PLD-Wrapper.py limits its use of system resources, in an attempt not to fry your machine. If you are running PLD-Wrapper in a subsystem, or a VM, it will similarly attempt not to use all of the subsystem/VM resources. Make sure you account for this when allocating resources to such a subsystem or VM.

## Acknowledgements

I am grateful to my father (who wished not to be named) for many helpful discussions regarding this program. In particular, the central idea for this wrapper, to monitor the output of PLD.jl to automate it and avoid stalling, was his.

I also thank the authors of PLD.jl for their support in this endeavour.

This program was written to complement a summer project supervised by Einan Gardi and funded by the University of Edinburgh School of Physics and Astronomy's Summer Vacation Scholarship.
