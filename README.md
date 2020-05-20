# ChocoboPy
Chocobo CFD Code - Python/Cython Version

## Dependencies
*ChocoboPy is written for Python 3.X, and will not run under Python 2.X*

ChocoboPy depends on the following Python packages (all of which can be installed via e.g. Anaconda if required).

* Cython
* Numpy
* PyInstaller

In addition, a C compiler is needed. This should be the same C compiler used to build Python on your system (usually gcc on Linux).

## Building ChocoboPy
ChocoboPy is built in two stages:

1. Compilation of Lagrangian hydrodynamics kernal using Cython.
2. Compilation of `chocobopy` executable using PyInstaller.

To make this process as easy as possible, a makefile is provided which carries out the build tasks, along with optional installation tasks. The usual way of using this is:

    make; make install

By default, the installation path is `/prod/ChocoboPy/bin`, with a simlink placed in `/prod/bin`, although these can be changed in the makefile.

## Running ChocoboPy
To run ChocoboPy, assuming its installation location is in your `$PATH`, simply type:

    chocobopy

This will assume that the following input files are located in the current directory:

| File        | Purpose                          |
|-------------|----------------------------------|
| input.in    | Main input file                  |
| mesh.key    | Mesh file (in LSDYNA KEY format) |
| eos_data.in | EoS data file                    |

These default names can be changed using the following command-line input options

| Option            | Purpose                                                                                                         |
|-------------------|-----------------------------------------------------------------------------------------------------------------|
| `--control`, `-c` | Sets name of main input file                                                                                    |
| `--mesh`, `-m`    | Sets name of mesh file                                                                                          |
| `--eos`, `-e`     | Sets name of EoS file                                                                                           |
| `--problem`, `-p` | Sets problem name. If specified, ChocoboPy will expect to find input files inside a directory with the same name |

Note that use of the `-p` switch is optional. For example, say you had the following input files in a subfolder as follows:

```
.
+-- sod
    +-- sod.in
    +-- sod.key
    +-- sod.eos.in
```

You could run this using either of the following commands:

    chocobopy -c sod/sod.in -m sod/sod.key -e sod/sod.eos.jl

or

    chocobopy -p sod -c sod.in -m sod.key -e sod.eos.jl

These is one important difference between the commands, however: The first command will put output inside a `results` subdirectory in the current directory, whereas the second will place them in a `results` subdirectory of the `sod` directory. The output location can be changed, however, using the `--output` (`-o`) option:

| Option            | Purpose                                                                                                         |
|-------------------|-----------------------------------------------------------------------------------------------------------------|
| `--output`, `-o`  | Sets location to place output files                                                                             |

## Control (input) file
The main control file takes the form of an INI-style configuration file. That is, inputs are of the form `variable = value`, and sit inside blocks demarked by a `[block name]` style header. Comments may be included in the file using the `#` character, and may appear inline, or as lines by themselves. The following example details all the current input options, with the values listed here also being the default values for each option.

    # Core problem control parameters
    [Control]
    t0 = 0.0            # Initial problem time, usually set to zero
    tend = 0.205        # Problem end time
    dtinit = 0.0001     # Initial timestep
    dtmax = 0.001       # Maximum permitted timestep
    growth = 1.02       # Maximum growth of dt in a single step

    # Options controlling artificial viscosity parameters (CL and CQ)
    [Q]
    CL = 0.5
    CQ = 0.75

    # This block is supported for reading, but not currently used
    [SphSod Mesh]
    X0 = 0.0
    X1 = 1.0
    Y0 = 0.0
    Y1 = 1.0
    nregions = 1
    meshx = 50
    meshy = 50

    # Debugging options
    [Debug]
    debug_step_count = 0    # If set > 0, run only this many steps

    # Material/Region linkage
    [Material]
    material_numbers = 1, 2     # Sets material number for each region

## Mesh file
Meshes should be generated using a mesh generator which supports the following:

* Generation of **quad** meshes
* Output in LSDYNA Key format

The GMSH generator (https://gmsh.info/) has been tested and demonstrated to produce meshes of the form required. Note that by default it generates triangular meshes, and some care is needed to generate a good quad mesh. The following GMSH script will generate a mesh suitable for running the Sod shock tube problem:

    // Mesh physical size
    lc = 0.02;
    l = 1.0;
    h = 0.1;

    // Mesh resolution
    xc = 20;
    yc = 10;

    // Create mesh geometry from Points and Lines

    // Create geometry points
    Point(1) = {0, 0, 0, lc};
    Point(2) = {l/2, 0,  0, lc};
    Point(3) = {l/2, h, 0, lc};
    Point(4) = {0, h, 0, lc};
    Point(5) = {l, 0, 0, lc};
    Point(6) = {l, h, 0, lc};

    // Create geometry lines
    Line(1) = {1, 2};
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 1};

    Line(5) = {2, 5};
    Line(6) = {5, 6};
    Line(7) = {6, 3};

    // Create surfaces. These produce regions in the mesh. Defined via loops of lines
    Curve Loop(1) = {4, 1, 2, 3};
    Curve Loop(2) = {-2, 5, 6, 7};
    Plane Surface(1) = {1};
    Plane Surface(2) = {2};

    // xc and yc refer to number of *nodes*, so add one to them
    xc++;
    yc++;

    // Create regions, each spanning half the X space
    Transfinite Surface {1} = {1, 2, 3, 4};
    Transfinite Line {1, 3} = xc/2.0;
    Transfinite Line {4, 2} = yc;

    Transfinite Surface {2} = {2, 5, 6, 3};
    Transfinite Line {5, 7} = xc/2.0;
    Transfinite Line {6} = yc;

    // Boundaries
    Physical Curve("YLOW") = {1, 5};
    Physical Curve("YHIGH") = {3, 7};
    Physical Curve("XLOW") = {4};
    Physical Curve("XHIGH") = {6};

    // Set meshing options
    Mesh.Algorithm = 8;
    Mesh.RecombinationAlgorithm = 3;
    Mesh.RecombineAll = 1;

    Mesh.SaveAll = 84;
    Mesh.SaveElementTagType = 1;
    Mesh.SaveTopology = 0;
    Mesh.SaveParametric = 0;
    Mesh.SaveGroupsOfElements = 1;
    Mesh.SaveGroupsOfNodes = 3;

    // Generate mesh
    Mesh 2;

    // Save mesh file
    Save "sod.key";

    Exit;

To correctly set the boundary conditions, the edges of the mesh **must** be included on `Physical Curve`s with the names *XLOW*, *XHIGH*, *YLOW* and *YHIGH*, as shown in the example.

## EoS file
The EoS file should consist of one line of data per material, with each line being a comma-separated list of the following parameters:

* Material Number, Initial Density, Initial Pressure, Gamma (for perfect gas EoS)

Lines beginning with a `#` are ignored. For example:

    # Format
    # Mat number, rho0, P0, gamma
    1, 1, 1, 1.4
    2, 0.125, 0.1, 1.4

