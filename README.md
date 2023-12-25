# 3D Simulation of Wheat Strands Dynamics with Frictional Contact in Wind-Driven Fields

1 - Rendering of a Strands as a system of rods

2 - Implementing the lagrangian dynamic of one Strands

3 - Testing the rendering and dynamic with multiple Strands

4 - forced wind : adding a sinusoidal forced wave to simulate the wind (naive approach)

6 - Critics on the simulation of Strands using pendulum and ressort -> trying something else and comparing

5 - friction and contact between Strands during the dynamic (first detect the collision THEN add a response) * -> then do 3D

7 - SoTA of wind chaotic dynamic .. (complex topic apparently)

![](video.mov)

## `Project Organization`

- `brin_blé.py`: Main source code file containing the simulation
- `geom/`: Folder containing geometry objects
- `graphics/`: Folder containing renderable objects, shader ...
- `mesh/`: Folder containing mesh obj
- `dynamics/`: Folder containing dynamics classes
- `pydfcp/`: Folder containing DFCP optimizer



## `Next step`
- Enhance collision detection speed by implementing optimized algorithms such as Sweep and Prune.
- Previously, we approximated a stream using an articulated system of multiple pendulums and springs. Now, let's explore another sophisticated approximation and compare the results.
- Expand the simulation from 2D to 3D.
