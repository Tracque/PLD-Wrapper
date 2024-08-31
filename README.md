
  

# PLD-Wrapper.py

  

A high-level wrapper, intended to make the package PLD.jl, more accessible and faster.

  

This program comes with a number of scripts, a handful of which are intended for users to interact with/modify:

  

The main program, which should be the most intuitive (if at times a little clunky) to use, is PLD-Wrapper.py.

The GUI-less version of the program, which is more stable and is intended to be more convenient for those who are confident with programming, or are running in a headless instance, is PLDManager.py.

The script intended for use if the program does not terminate and you thus end execution early, is PLDCleanup.py. (currently this is only compatible with use of PLD for a Feynman diagram)

The script to calculate the PLD with a custom (manually input) polynomial, is PLDCustom.py.

Finally, the "Slurm Scripts" show an example of how the PLD-Wrapper program may be adapted for use with a compute cluster.

  

## System requirements

  

The program has not yet been tested on many systems, so the minimum specs are unknown. The recommendation is that your system can comfortably run PLD.jl, with at least 2/3 of the system resources unoccupied.

  

The main bottleneck for this program is the amount of RAM that is available. For reference, on a system with a 24 core Intel® Core™ i9-7920X CPU and 32 GiB of RAM, my experience is that the CPU usage has never climbed much higher than 50%, whilst the RAM usage can often reach 100% for challenging diagrams if one removes the limits on the program.

  

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

Pkg.add("LinearAlgebra")

...

```

  

One issue you may encounter is with the package OSCAR failing to precompile. This is likely due to you having failed to properly configure the julia installation, for example not having properly set the PATH variable. The package manager is therefore fetching very old versions of some dependencies (namely, cxxmake.jl and thus Singular.jl, Polymake.jl and others). If you want a quick fix without really bothering with diagnosis, running ```Pkg.update()``` will now fetch working versions of the dependencies, solving this issue.

  

Note: PLD.jl is NOT a package that can be added using Pkg. You do not need to worry about downloading PLD separately, as it is included in this repo. (in fact, if you download it yourself, you WILL encounter errors, as OSCAR has undergone some slight syntax changes which required edits to PLD.jl and some of the functionality of PLD-Wrapper also required edits to PLD.jl)

  

PLD-Wrapper v1.2.0 was last tested to be compatible with the current versions of the above packages on the 25th of August 2024. As a reference, in case something breaks in future, you may want to compare the result of the command ```Pkg.status()``` to the below:

  

```

Status `~/.julia/environments/v1.10/Project.toml`

⌅ [fb37089c] Arblib v0.8.1

[7d9fca2a] Arpack v0.5.4

[01680d73] GenericSVD v0.3.0

[f213a82b] HomotopyContinuation v2.10.0

[f1435218] Oscar v1.0.4

[37e2e46d] LinearAlgebra

[de0858da] Printf

Info Packages marked with ⌅ have new versions available but compatibility constraints restrict them from upgrading. To see why use `status --outdated`

```

  

Step 4: Make sure you have a working version of python3 installed, as well as the library psutil, which the wrapper uses to avoid frying your machine and (if you intend on using the GUI) pywebview. You can download these by running:

  

```

pip install psutil

pip install webview

```

  

Or,

  

```

sudo apt-get install python3-psutil

sudo apt-get install python3-webview

```

  

------

  

Note: pywebview, by default, will attempt to open the GUI window using GTK. The program will still run successfully even without GTK installed. If you attempt to run PLD-Wrapper.py without GTK, you will see error messages as follows:

  

```

[pywebview] GTK cannot be loaded

Traceback (most recent call last):

File "/usr/lib/python3/dist-packages/webview/guilib.py", line 16, in import_gtk

import webview.platforms.gtk as guilib

File "/usr/lib/python3/dist-packages/webview/platforms/gtk.py", line 26, in <module>

gi.require_version('Gtk', '3.0')

File "/usr/lib/python3/dist-packages/gi/__init__.py", line 126, in require_version

raise ValueError('Namespace %s not available' % namespace)

ValueError: Namespace Gtk not available

QStandardPaths: wrong permissions on runtime directory /run/user/1000/, 0755 instead of 0700

```

Despite the errors, the program will still run and give the correct output. If, however, it bothers you to have such errors occuring, you can find the installation instructions [here.](https://www.gtk.org/docs/installations/index)

  

And with that, you should be up and running! Simply navigate to the directory where you installed PLD-Wrapper and run the command:

  

```

python3 PLD-Wrapper.py

```

  

## Usage Instructions

  
### PLD Wrapper

PLD-Wrapper.py will open a local webpage which provides input fields and instructions for PLD. It also can produce visualisations of your inputs, to hopefully help you check that you are inputting what you think you are. Although self edges (edges of the form [i,i]) are supported by PLD.jl and will be calculated properly, they will not be drawn by the visualisation. Note that the visualisations are NOT supposed to be literature ready diagrams. 

PLD-Wrapper.py will assume it is allowed to use most of your machine's resources. This is an intentional design choice. If you need to protect these resources, for example because you are running on a shared compute cluster, you should use PLDManager.py instead.

  
### PLD Manager
I generally recommend that users interact primarily with this script, *especially if you are investigating a challenging diagram*.

To use PLDManager.py, you will have to adjust the inputs by directly modifying the main() function. PLDManager.py allows users to set a custom memory limit by passing an additional argument at runtime:

```
python3 PLDManager.py MEM_LIMIT
```
There is support for the memory limit to be given in bytes, KB, MB or GB. (you should denote the memory as something like 23M, or 10G, without the B) 

Note: Even if no memory limit is given, PLD-Wrapper.py and PLDManager.py limit the use of system resources, in an attempt not to fry your machine. If you are running them in a subsystem, or a VM, it will similarly attempt not to use all of the subsystem/VM resources. Make sure you account for this when allocating resources to such a subsystem or VM.

Specifically, CPU usage and RAM usage are both limited, although it is not theoretically impossible for either to become problematic if many numeric processes use much more resources than average.

  
 ### PLD Cleanup

PLDCleanup.py exists to account for situations where PLD-Wrapper.py or PLDManager.py are stopped prematurely and you do not wish to continue with the calculation. In this case, you've probably been left with a lot of junk output/input files clogging the directory you chose. Here, you should again directly modify the main() function to give the correct inputs. Running PLDCleanup.py should then compile whatever output there is and sort the output files for you.

It will also delete the extra output/input files once finished.

### PLD Custom

PLDCustom.py allows users to input a custom polynomial in place of the Graph polynomial that would be calculated by PLD.jl in the normal case. Users must be careful to format the polynomial in the same way they format the inputs for the variables and parameters of the polynomial. (in the maths sense, so for a graph polynomial, the Schwinger parameters are the variables and the kinematic variables are the parameters!)

This means if you used subscripts or brackets etc. in the polynomial, you must also do so in the input array for the varaibles/parameters.

## Advanced Usage 

### Sysimages
Those who are experienced with using julia may be familiar with the practice of creating sysimages, which can significantly decrease program load times. On my system, it provides a modest improvement, reducing the load time of the julia scripts from around 1 minute to around 20 seconds. This reduces the overall overhead of the program by around 30-40%, but your results may vary. I will not provide a step-by-step tutorial here, but will instead refer you to the PackageCompiler.jl [documentation](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html) if you wish to use this advantage.

Simply create a sysimage of a script containing the preamble (the using Oscar using PLD parts at the start of the julia scripts) and then go to the python scripts to find + replace the characters *"julia"* with *"julia", "--sysimage", "YOURSYSIMAGE.so"*. (including the quote marks)

Make sure you do this in every script to get the most out of it, so if you are using PLDManager.py, don't forget to also make this change in PLDUtils.py!

### Slurm Scripts
The provided scripts are given as an example only. They were written to be run on a compute cluster using the SLURM manager. I expect users who are attempting to run large calculations on such clusters are able to modify these scripts or write their own according to their use case.

  
### Julia Files
At no point when using these programs for their intended purpose should you have to touch the Julia files. However, feel free to tinker with things if you wish to add your own functionality, or if you spot a problem. (In the latter case, do let me know so I can fix it for everyone else too!) 

Something that hasn't yet been added is support for the much slower "high precision" mode of PLD.jl. This is relevant for numeric calculations, which may suffer rounding errors/tolerances which can result in  wonky coefficients such as 109656/73105 = 1.4999794... instead of 3/2. 

  

## Known Issues

The cleanup script is currently only compatible with normal use of PLD. If there is demand, I will make available a version for use with custom polynomials.

Whilst I have tried to ensure the input processing for this program is as robust as possible, inevitably, one must choose a format and stick to it. If you are sure you installed the program correctly, but are having problems running a specific diagram, I suggest to check the formatting of your inputs as a first step.

## Related Projects

Development of a companion tool to use PLD.jl for expansion by regions has come up as a potentially useful addition since PLD.jl in its current form only supports strict substitutions, so can only obtain strict limits and not expansions near the limit.
This tool, if it is made, will be made open source in another Github repository, the link to which will appear here when released...

  

## Acknowledgements

  

I am grateful to my father (who wished not to be named) for many helpful discussions regarding this program. In particular, the central idea for this wrapper, to monitor the output of PLD.jl to automate it and avoid stalling, was his.

  

I also thank the authors of PLD.jl for their support in this endeavour.

  

This program was written to complement a summer project supervised by Einan Gardi and funded by the University of Edinburgh School of Physics and Astronomy's Summer Vacation Scholarship.

  

## Contact

  

This repo is being maintained by Tristan Jacquel (Github: [Tracque](https://github.com/Tracque)). If you would like to raise an issue with the program, feel free to do so here on Github, but you are more likely to get a timely response if you contact me via email: s2146323@ed.ac.uk (or t.y.jacquel@sms.ed.ac.uk which also goes to me)

  

If, for some reason, you cannot reach me through email or on Github, then perhaps try my personal email: tristan@jacquel.net
