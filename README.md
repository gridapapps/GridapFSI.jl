# GridapFSI

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridapapps.github.io/GridapFSI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridapapps.github.io/GridapFSI.jl/dev)
[![Build Status](https://travis-ci.com/gridapapps/GridapFSI.jl.svg?branch=master)](https://travis-ci.com/gridapapps/GridapFSI.jl)
[![Codecov](https://codecov.io/gh/gridapapps/GridapFSI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridapapps/GridapFSI.jl)

This is an application repository with a collection of drivers for the simulation of Fluid-Structure Interaction (FSI) problems. It is based on [Gridap](https://github.com/gridap/Gridap.jl), a package for grid-based approximation of PDEs with Finite Element, and [GridapODEs](https://github.com/gridap/GridapODEs.jl), a package that provides time integration tools for Gridap.

## Installation
This is a Julia package that can be easily installed using the `Pkg` manager. Open the Julia REPL, type `]` to enter into the package mode, and install as follows
```julia
pkg> add https://github.com/gridapapps/GridapFSI.jl
```

## Usage
First, include the main GridapFSI module:
```julia
using GridapFSI
```
Then, execute the `main()` function especifying the driver name and the optional arguments. By default the `elasticFlag` driver is called, with the default parameters corresponding to the FSI2 test case. Other drivers can be called with the following sintax:
```julia
main(problemName="<driverName>",kwargs...)
```
Currently the following drivers are implemented:
  - `problemName="analytical"`: A monolithic FSI simulation using a manufactured analytical solution.
  - `problemName="elasticFlag"`: A monolighic FSI simulation of the elastic flag behind a cylinder benchmark.

## Examples
### Elastic flag after a cylinder: FSI-1
The FSI1 test case for the elastic flag benchmark proposed by S. Turek and J. Horn can be run with the following command:
```julia
output = main(
  Um=0.2,
  Re=20,
  rho_s=1.0e-3,
  strategy="biharmonic",
  alpha_m=1.0e-6,
  alpha_m_weight="volume",
  model="models/elasticFlagFine.json",
  dt=1.0,
  tf=25.0,
  theta=0.6
  )
```
This call will run the FSI1 test case using a θ-method with a time step size of `dt=1.0`, using the Biharmonic mesh motion strategy with a cell volume weighted constant, and a discrete model defined in `models/elasticFlagFine.json`. In this test the displacement, velocity and pressure fields are outputed to `.vtu` files located in the automatically generated folder `fsi-results`. The result `output` consists on a tuple with the values of time, drag and lift forces. For this specific case we have the following values at the final time:

  \# elements | \# DOFs | x-displ of A | y-displ of A | Drag force | Lift force 
  --- | --- | --- | --- | --- | --- 
  7342 | 152,745 | 2.257e-5 | 8.728e-4 | 14.215 | 0.6803
  Reference | | 2.27e-5 | 8.209e-4 | 14.295 | 0.7638

  ![](/models/velFSI1.png)

## Contributing
Contributions with the definition of new drivers and additional FSI formulations are welcome. The repository is organized as follows:
  - `Gridap.jl`: module with the main function and inclusion of submodules.
  - `FSIDrivers.jl`: module with a list of user-defined drivers. Each driver must implement the `execute` function with the corresponding problem name (`<driverName>`). 
    ```julia
    function execute(problem::Problem{:<driverName>}; kwargs...)
      # user-defined driver
    end
    ```
  - `WeakForms.jl`: module with the definition of the weak forms for all required residuals and Jacobians appearing in the FSI formulation.

Please, follow the [Gridap contributing guides](https://github.com/gridap/Gridap.jl/blob/master/CONTRIBUTING.md) when developing new features in the repository.

## How to cite
If you have used these drivers in a scientific publication, please cite Gridap library and the stabilization method as follows:

```
@software{gridap_project,
  author       = {Francesc Verdugo and
                  Santiago Badia and
                  Víctor Sande and
                  Alberto F. Martin and
                  Oriol Colomés and
                  Jesús Bonilla},
  title        = {Gridap.jl},
  year         = 2020,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.3934468},
  url          = {https://doi.org/10.5281/zenodo.3934468}
}
```
