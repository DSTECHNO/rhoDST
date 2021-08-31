![enter image description here](https://user-images.githubusercontent.com/30440239/129488825-eb1f5388-fe71-45bf-9f60-6a9d69466836.jpg)

# rhoDST
This repository provides a density based solver **rhoDST(c)** for steady and unsteady simulation of high speed compressible flows over aeronautical vehicles. The proposed solver **rhoDST(c)** has been developed based on  foam-extend 4.1  libraries for the fully-coupled solution of governing equations on highly skewed unstructured meshes. 
 

# Features

 - Continuity, momentum and energy equations are implicity coupled in a single matrix.  The block-matrix can be solved using block-matrix solvers such as GMRES, BICG and multigrid methods, which are availble in foam-extend libraries. 
 - Local time stepping is used for the rapid convergence to the steady-state solution. 
 - Shock-discontinuities can be captured using Kurganov, AUSM, AUSM+, AUSM+up, HLL, HLLC and HLLCP shock-capturing methods. 
 - Modifications of SpalartAllamaras turbulence model are incorporated to the turbulence libraries. Simpler Rotation-Curvature and Spalart-Shur rotation/curvature correction were implemented to account for rotation and curvature effects.  Negative Spalart-Allmaras is implemented to resolve numerical issues near the interface between turbulent and irrotational regions. 
 

# Developers

Ender Demirel, PhD

Aras Dogan


