# ArterialWaveSolutions
### An Arterial Parameter Estimation Technique using a 1D Frequency Domain Model

This project computes a set of known parameters given 3D Numerical
simulation measurements at the outlets of an arterial system. These parameters
generally include the length and inlet radius of a vessel, as well as
a stiffness factor, and 3 Windkessel parameters. The 3D numerical simulation
measurements include pressure and volume flow rates at various locations
within the system in question.

## Directory Structure

The folders in this repository are formatted to contain different
kinds of scripts in each directory. 

- The root directory contains scripts that run the parameter estimation
technique on each of the five cases of arterial networks.
- The `Algorithms/` directory contains optimisation algorithms used in
the estimation procedure, as well as a framework to apply them seamlessly
in code in the `OptimisationProblem.m` file. The preprocessing methods
used are also contained in this folder, which includes the `MinMaxScaler.m`
file.
- The `Data/` directory contains output data from the full aorta
parameter estimation procedure.
- The `Figures/` directory contains plots of the pressure and
volume flow rate at various locations of three kinds of arterial networks.
These include the upper thoracic aorta (UTA), single bifurcation and
the full aorta.
- The `Models/` directory contains the model used in each test case present
in `PreviousCode/`. It is essentially the same as scripts in `PreviousCode/`
but in an object-oriented (class) format instead of individual scripts.
- The `NewCode/` directory contains scripts similar to `PreviousCode/`
except using the new implementation of the model from `VesselModels/`.
- The `PreviousCode/` directory contains code previously written by
Georgios Efstathiou that simulates blood flows in 1D flexible tapered
vessels in the frequency domain.
- The `PreviousCode/automeris/` directory contains the 3D numerical simulation
measurements used in the parameter estimation technique.
- The `VesselModels/` directory contains a revised implementation of the model
previously implemented in the `PreviousCode/` scripts in an object-oriented
format.

## Running Instructions

To run any of the parameter estimation scripts install Matlab R2022a then
run any of the scripts in the root directory. Set the variable
`measurements` to `false` for determining the known parameters using the
corresponding 1D outputs or `true` to estimate the same parameters
by fitting the 1D model to the 3D numerical simulation results. 