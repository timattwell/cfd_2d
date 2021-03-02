A C++ based CFD code written for my high performance computing university module. Currently bad. Will be better soon (TM).

Planned improvements/changes (in no particular order):
1) Optimisations to the numerical scheme used by the solver -
  a) Experiment with alternative schemes.
  b) Implement LAPACK functions into solver for increased calculation speed.
2) Implement parallel operations, splitting the domain using a cartesian grid system and solving across multiple processors using mpi.
3) Build in convenient plotting functionality
4) Expand the makefile to include test cases with courser domains and shorter times. Really should've been done in the initial project, but I was lazy.
5) Add in capability to include solid obstructions to the flow.
6) Include approximations to account for viscous effects at walls.
7) Experiment with changing mesh size or scheme near complex geometries or walls.

Things I might like to do one day maybe:
1) Develop a user friendly GUI
2) Actually comment things

~ Tim
