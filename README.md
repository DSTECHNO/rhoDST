![enter image description here](https://user-images.githubusercontent.com/30440239/129488825-eb1f5388-fe71-45bf-9f60-6a9d69466836.jpg)

# rhoDST
This repository provides a density based solver **rhoDST** for steady and unsteady simulation of high speed compressible flows over aeronautical vehicles. The proposed solver **rhoDST** has been developed based on  foam-extend 4.1  libraries for the fully-coupled solution of governing equations on highly skewed unstructured meshes. 

![BG1](https://user-images.githubusercontent.com/30440239/131627275-4d343da5-bcdc-4502-8275-a7f8de4ea75a.jpg)
 

# Features

 - Continuity, momentum and energy equations are implicity coupled in a single matrix.  The corresponding matrix can be solved using block-matrix solvers such as GMRES, BiCGStab and multigrid methods. 
 - Local time stepping is used for the fast convergence to the steady-state solution. 
 - Shock-discontinuities can be captured using Kurganov, AUSM, AUSM+, AUSM+up, HLL, HLLC and HLLCP shock-capturing methods. 
 - Modifications of SpalartAllamaras turbulence model are incorporated to the turbulence libraries. Simpler Rotation-Curvature and Spalart-Shur rotation/curvature correction are implemented to account for rotation and curvature effects.  Negative Spalart-Allmaras is implemented to resolve numerical issues near the interface between turbulent and irrotational regions. 
 - Characteristic far-field boundary conditions are used at the outer boundaries for subsonic, transonic and supersonic freestream flows. 
 
# Installation Guide
The dependencies for the solver **rhoDST** are the same as for foam-extend-4.1 and by installing foam-extend-4.1 all dependencies would be satisfies. Before compiling the solver **rhoDST**, environment variables must be set for foam-extend-4.1.

```bash
cd $WM_PROJECT_USER_DIR
git clone https://github.com/DSTECHNO/rhoDST rhoDST
cd rhoDST
./Allwmake
```
# Examples
* Transient simulation of supersonic flow over a forward facing step for Ma = 3 
![bfs](https://user-images.githubusercontent.com/30440239/131734244-9f159e29-19da-4087-aebf-21270e1e7633.png)

* Transient simulation of supersonic flow over a forward facing step using unstructured mesh for Ma = 3 
![bfsUns](https://user-images.githubusercontent.com/30440239/131735815-6307e367-1bbd-4b34-a196-40ec75d1d2ab.png)

* Steady-state simulation of transonic flow over RAE2822 airfoil for Ma = 0.729
![pressure](https://user-images.githubusercontent.com/89465885/131667405-db903eff-9901-481e-b589-f4aba4e09966.png)

* Steady-state simulation of transonic flow over ONERA M6 Wing for Ma = 0.8395
![onera](https://user-images.githubusercontent.com/30440239/131879893-f28399bb-7e10-4537-a0e5-99221c805d0f.png)

* Steady-state simulation of transonic flow over a full aircraft model for Ma = 0.9

![mesh](https://user-images.githubusercontent.com/89465885/131667919-018a8a7f-f218-4d2b-a677-4a5eb141a59f.png) |  ![nut](https://user-images.githubusercontent.com/89465885/131667922-820729ee-6d01-4dc7-acb1-1374a36cfb6e.png)
--- | ---

# Developers

- Aras Dogan
- Prof.Dr. Ender Demirel

# License
Copyright 2020 Design and Simulation Technologies Inc. 

Distributed using the GNU General Public License (GPL), version 3; see the LICENSE file for details.
