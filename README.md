# SnowFlake

## Parameters

| Parameter  | Description                                 |
|------------|---------------------------------------------|
| size       | Dimension of the lattice                    |
| iterations | Number of iterations                        |
| rho        | Initial uniform vapor density on the domain |
| beta       | Attachment parameter                        |
| alpha      | Knife-edge instability parameter (liquid)   |
| theta      | Knife-edge instability parameter (vapor)    |
| kappa      | Freezing parameter                          |
| mu         | Evaporation parameter                       |
| gamma      | Melting parameter                           |
| sigma      | Random perturbation parameter               |

## Building

On Uniy-like systems simply run `./compile.sh`. This will create a dist directory with the compiled binary.

Build dependencies:

- CMake
- Make

## Usage

1. Define parameters in `main.cpp`
2. Build
3. `cd dist`
4. `./snowflake`

The program saves the snowflake into four binary matrices: vapor mass (double), solid mass (double), liquid mass (double) and the snowflake's interior (boolean).

### Generating figures

The Python file `hex_grid.py` contains a function to plot/export the different masses on the domain. A hexagonal grid is used to represent the underlying crystalline structure. A snapshot of the snowflake can also be exported.
