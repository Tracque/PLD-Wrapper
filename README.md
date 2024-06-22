
# PLD-Wrapper.py

A high-level wrapper, intended to make the package PLD.jl, more accessible and faster.

This program comes with a number of scripts, 3 of which are intended for users to interact with/modify:

The main program, which should be the most intuitive (if at times a little clunky) to use, is PLD-Wrapper.py. 
The GUI-less version of the program, intended to be more convenient for those who are confident with programming, or are running in a headless instance, is PLDManager.py.
The script intended for use if the program does not terminate and you thus end execution early, is PLDCleanup.py.

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

One issue you may encounter is with the package OSCAR failing to precompile. This is due to the package manager fetching bad outdated versions of some dependencies (namely, cxxmake.jl and thus Singular.jl, Polymake.jl and others). For whatever reason, simply running ```Pkg.update()``` will now fetch working versions of the dependencies, solving this issue. (Yes, you can "update" these packages which you added seconds earlier)

Note: PLD.jl is NOT a package that can be added using Pkg. You do not need to worry about downloading PLD separately, as it is included in this repo. (in fact, if you download it yourself, you WILL encounter errors, as OSCAR has undergone some slight syntax changes which required edits to PLD.jl)

PLD-Wrapper v1.0.1 was last tested to be compatible with the current versions of the above packages on 22/06/2024. As a reference, in case something breaks in future, you may want to compare the result of the command ```Pkg.status()``` to the below:

```
Status `~/.julia/environments/v1.10/Project.toml`
⌅ [fb37089c] Arblib v0.8.1
  [7d9fca2a] Arpack v0.5.4
  [01680d73] GenericSVD v0.3.0
  [f213a82b] HomotopyContinuation v2.9.3
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

PLD-Wrapper.py will open a local webpage which provides input fields and instructions for PLD. It also can produce visualisations of your inputs, to hopefully help you check that you are inputting what you think you are. Note that the visualisations are NOT supposed to be literature ready diagrams.

If you wish to use PLDManager.py, you will have to adjust the inputs by directly modifying the main() function. 

PLDCleanup.py exists simply to account for situations where PLD-Wrapper.py or PLDManager.py are stopped prematurely. In this case, you should again directly modify the main() function to give the correct inputs. Running PLDCleanup.py should then compile and sort the output files for you.

At no point when using these programs should you have to touch the Julia files.

Note: PLD-Wrapper.py and PLDManager.py limit the use of system resources, in an attempt not to fry your machine. If you are running them in a subsystem, or a VM, it will similarly attempt not to use all of the subsystem/VM resources. Make sure you account for this when allocating resources to such a subsystem or VM.

## Known Issues

The formatting of output gets messed up if you run the cleanup script multiple times on the same files. (although I am tempted to say there should virtually never be any reason to do such a thing)

## Acknowledgements

I am grateful to my father (who wished not to be named) for many helpful discussions regarding this program. In particular, the central idea for this wrapper, to monitor the output of PLD.jl to automate it and avoid stalling, was his.

I also thank the authors of PLD.jl for their support in this endeavour.

This program was written to complement a summer project supervised by Einan Gardi and funded by the University of Edinburgh School of Physics and Astronomy's Summer Vacation Scholarship.

## Contact

This repo is being maintained by Tristan Jacquel (Github: [Tracque](https://github.com/Tracque)). If you would like to raise an issue with the program, feel free to do so here on Github, but you are more likely to get a timely response if you contact me via email: s2146323@ed.ac.uk (or t.y.jacquel@sms.ed.ac.uk which also goes to me)

If, for some reason, you cannot reach me through email or on Github, then perhaps try my personal email: tristan@jacquel.net
