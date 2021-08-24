## Swimming in potential flow

This repository contains Julia code for simulating a freely deforming Joukowski foil that sheds vortices discretely in time, effectively modeling planar fish-like swimming. The dynamics are derived using a momentum-conserving model (see [this](https://www.ideals.illinois.edu/handle/2142/83903) or [this](http://scottdavidkelly.wdfiles.com/local--files/start/kpx12.pdf) for more information concerning how to derive the dynamics). The idea is that at each discrete instance in time, a [point vortex](http://web.mit.edu/fluids-modules/www/potential_flows/LecturesHTML/lec1011/node21.html) is placed into the fluid corresponding to the Kutta condition, which enforces that the local velocity field is zero at the preimage of the Joukowski foil. 

Several implementations of this have existed over the course of several graduate students from [Professor Scott Kelly's lab](http://www.kellyfish.com/) at UNC Charlotte. A Matlab version of this code has survived over the years and, to keep with the times, the present repository is my implementation of it in Julia.

### Current functionality and running the simulation

Currently the simulation works with no point vortices already present in the flow and with zero initial conditions for the Joukowski foil states. This means that ``computeInitialConditions.jl`` is still being worked on.

Assuming you already have [Julia installed](https://julialang.org/), to run the simulation, download and navigate to the ```vortexshedding``` directory and execute ```joukowskiMain.jl```in the REPL. Note that it is note yet clear to me how to vectorize this code (this is a current project of mine), so it may take several minutes to run a simulation. The default parameters are ```T = 4.0``` for the duration of the simulation and ```maxNumberOfVortices = 1000``` as a maximum number of vortices which can be shed throughout that duration. Feel free to modify these variables to achieve shorter or longer simulations.

It is also possible to modify the amplitude and frequency of flapping for the simulation. To adjust these values, change ```ua``` and ```uf``` in ```joukowskiMain.jl```, where ```ua``` corresponds to the amplitude and ```uf``` the frequency.

### Animating the resulting simulation

Simply run ```animateJoukowski.jl``` to generate a series of png images. Afterward, navigate to the directory containing those images and execute the following in the terminal. But first, check what the size of the output image is and replace ```1200x1200``` with the appropriate image size.

```
ffmpeg -r 60 -f image2 -s 1200x1200-i ./filename%04d.png -vcodec libx264 -crf 1 -pix_fmt yuv420p video_name.mp4
```

 Note that you may need to modify the ```%04d``` part of ```./filename%04d.png``` since you may not have images numbering into four digits.

### To do

1. Implement arbitrary initialization via ```computeInitialConditions.jl```.
2. Implement the PID controller detailed [here](http://scottdavidkelly.wdfiles.com/local--files/start/kpx12.pdf).
3. Replace the flexible foil model with a rigid foil with an internally actuated rotor.