# Supersonic Ramjet CFD Solver

A finite volume Euler solver for axisymmetric supersonic ramjet simulation.

1.  Open `Ramjet_Config.m` to set Mach number and Heat Release.
2.  Run `Solver_Main.m`.
    * Wait for 2500 iterations.
3.  Run `Performance_Analysis.m` to see Thrust, Isp, and Mass Conservation error.
4.  Run `Final_Visuals.m` to generate contour plots.

## Files
* `Solver_Main.m`: The core RK2 time marching engine.
* `Ramjet_Config.m`: Simulation parameters (Mach, Geometry, Heat).
* `Grid_Gen_Combined.m`: Generates the mesh.
* `Compute_Sources.m`: combustion and geometric terms.
* `Flux_Roe.m`: Computes inviscid fluxes using Roe's method.

## Notes
* **Do not modify grid dimensions** in `Grid_Gen` without updating `Config`.
* If the solver crashes, reduce `Heat_Release_Max` in `Ramjet_Config.m`.
* Visuals are cropped to `y=0.8m` to focus on the engine core.