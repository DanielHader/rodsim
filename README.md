# rodsim

Author: Daniel Hader

Last Updated: 23 Jun 2025

RodSim is software for performing simulations of the [abstract Slat Assembly Model](http://self-assembly.net/wiki/index.php?title=Abstract_Slat_Assembly_Model_(aSAM)) and the [kinetic Slat Assembly Model](http://self-assembly.net/wiki/index.php?title=Kinetic_Slat_Assembly_Model_(kSAM)).

## Building

There's an included CMakeList.txt file so it can be built using cmake so long as you have the dependencies installed.

```bash
mkdir build
cd build
cmake ..
make -j4
```

Once that's done, you'll have an executable called `rodsim_cli` which is the simulator.
