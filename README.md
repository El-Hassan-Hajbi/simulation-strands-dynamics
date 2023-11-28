# 3D Simulation of Wheat Strands Dynamics with Frictional Contact in Wind-Driven Fields

1 - Rendering of a strain as a system of rods

2 - Implementing the lagrangian dynamic of one strain

3 - Testing the rendering and dynamic with multiple strains

4 - forced wind : adding a sinusoidal forced wave to simulate the wind (naive approach)

6 - Critics on the simulation of strains using pendulum and ressort -> trying something else and comparing

5 - friction and contact between strains during the dynamic (first detect the collision THEN add a response) * -> then do 3D

7 - SoTA of wind chaotic dynamic .. (complex topic apparently)


# Questions : 
- on a vu que la dynamique d'un systeme de pendule est chaotique, comment remedier a ce pb pour notre approximation d'un brin de blé
- doit on ajouter un terme de frottement : Rayleigh
- on a besoin d'aide pour trouver une réponse adéquate à la collision pour notre systeme de brins de blé
- Complexité de la detection de collision 