# HTC-Based-Transient-1D-FEM-for-Cooling-Curve-Simulation-with-and-without-Phase-Transformation
A 1D transient FEM code for cooling curve simulation using the Heat Transfer Coefficient (HTC), incorporating temperature-dependent thermal properties and metallurgical phase transformations. Supports simulations with and without phase transformation, including latent heat effects.

Overview

This code provides a numerical framework for transient thermal analysis of steel subjected to convective cooling. The temperature field is calculated using a 1D finite element formulation, while temperature-dependent material properties and metallurgical transformations are incorporated into the thermal solution.

Key Features

Transient 1D Finite Element Method (FEM)

Heat Transfer Coefficient (HTC)-based convective boundary condition

Cooling curve simulation

Temperature-dependent thermal properties

Simulations with and without metallurgical phase transformation

Ferrite, pearlite, bainite, and martensite transformation modeling where applicable

Schiel-sum based kinetics for diffusional transformations

Latent heat contribution during phase transformation

Hypo- and hypereutectoid steel treatment according to model inputs

Export of temperature and phase-fraction histories

Numerical Approach

The transient thermal problem is solved using a 1D FEM formulation with linear finite elements and time integration using a theta-method.

For simulations including phase transformation, thermal properties are updated according to the instantaneous phase fractions. Latent heat released during transformation is incorporated as a heat source term.

Separate Schiel sums are maintained for ferrite, pearlite, and bainite transformations. Martensitic transformation is treated separately as a temperature-dependent transformation.

Phase Transformation Model

The model uses critical transformation temperatures and TTT-based kinetic data supplied through the input data.

Hypereutectoid

Austenite

Pearlite

Bainite

Martensite

Hypoeutectoid

Austenite

Ferrite

Pearlite

Bainite

Martensite

Inputs

The model requires preprocessed data including:

Mesh connectivity and nodal coordinates

Initial austenite fraction

Critical transformation temperatures

TTT transformation data

Temperature-dependent density

Thermal conductivity

Specific heat

Latent heat data

Initial temperature

Ambient temperature

Heat Transfer Coefficient (HTC)

Time-step size

Outputs

The simulation generates time-dependent results including:

Temperature field and cooling curves

Austenite fraction

Ferrite fraction

Pearlite fraction

Bainite fraction

Martensite fraction

Schiel-sum histories

The results can be used to compare thermal behavior with and without phase transformation.

Applications

The code is intended for numerical investigation of:

Transient cooling behavior of steel

HTC-based cooling curve analysis

Effects of phase transformation on cooling curves

Latent heat effects during metallurgical transformation

FEM-based thermal and metallurgical simulations

Notes

Material-specific results depend on the accuracy and range of the supplied thermophysical and transformation data.

Author

Augustine Samuel
