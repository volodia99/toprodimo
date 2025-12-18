# toprodimo

***Command-line interface to interpolate and convert 2D (r, theta) datasets from simulations outputs (readable with [nonos](https://github.com/la-niche/nonos)) to ProDiMo models. Based on the [ProDiMo 2D interface notebook](https://prodimopy.readthedocs.io/en/stable/notebooks/interface2D.html) by Christian Rab***

## Development status

To be implemented : 
- Take care of the 2D distributions of dust fluids described by their size, when it will be fully implemented in ProDiMo
- Add flexibility on how the temperature is computed (for now the temperature in the simulation is directly given to ProDiMo, using a fixed solar mean molecular weight of 1.37)
- For now, `toprodimo` takes only the upper half of a 2D (`r`, `theta`) disk, if `theta` is symmetric compared to the midplane. Could add some flexibility if the user wants to focus on the lower half of the disk, or to average azimuthally a 3D (`r`, `theta`, `phi`) disk, or even add a prescription to extend vertically a 2D (`r`, `phi`) disk (e.g., with vertical hydrostatic equilibrium).  

## Installation

We recommend to install this repo using the package and project manager `uv`. See the [documentation](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer) to install `uv` on your system, then clone this repository.

## Run the CLI

After `cd lick_template`, you can run the CLI inside the project's virtual environment, via the following:

```shell
uv run toprodimo PATH_TO_SIMULATION_FILE -UNIT_LENGTH float -UNIT_MASS float -from PATH_TO_INIT_PRODIMO_MODEL_DIR -to PATH_TO_PRODIMO_MODEL_DIR
```

To get help, run
```shell
uv run toprodimo -help
```

```shell
usage: toprodimo [-h] [-dir DIRECTORY] [-UNIT_LENGTH UNIT_LENGTH] [-UNIT_MASS UNIT_MASS]
                 [-mask_inside MASK_INSIDE] [-from_pmdir INIT_PRODIMO_MODEL_DIRECTORY]
                 [-to_pmdir PRODIMO_MODEL_DIRECTORY] [-plot]
                 filename

positional arguments:
  filename              name of the output file to be read by ProDiMo.

options:
  -h, --help            show this help message and exit
  -dir DIRECTORY        location of output files and param files (default: '.').
  -UNIT_LENGTH UNIT_LENGTH
                        code unit of length [au].
  -UNIT_MASS UNIT_MASS  code unit of mass [solMass].
  -mask_inside MASK_INSIDE
                        mask the velocities inside given radius [inner edge unit]. put 0 for no
                        masking. (default: 1.2).
  -from_pmdir, -from INIT_PRODIMO_MODEL_DIRECTORY
                        location of the initialized prodimo model from which to extract ProDiMo.out.
                        Works only with -to_pmdir.
  -to_pmdir, -to PRODIMO_MODEL_DIRECTORY
                        location of the prodimo model on which the simulation model is then
                        interpolated. Works only with -from_pmdir.
  -plot, -p             rough plotting procedure to check the validity of toprodimo results.
```

## Additional arguments

- `mask_inside` removes the contribution of the radial and vertical velocities close to the grid's inner edge, to avoid some spurious effects in ProDiMo. By default we cancel these velocity components in a band that extends from r_inner to 1.2*r_inner.
- `plot` performs creates 3 .pdf files with some plots of the density/velocities/temperature
    - simulation.pdf: from the simulation file, with some post-processing (e.g., removing for all the fields the region inside the cylindrical radius corresponding to the inner edge).
    - prodimo.pdf: from the ProDiMo model, ready to be run with ProDiMo.
    - compare_simulation_prodimo.pdf: look at the 1D density in the midplane and vertically at R=UNIT_LENGTH. 

## Remarks

In order for the procedure to work, you need to keep in mind that:
- `toprodimo` needs the typical UNIT_LENGTH and UNIT_MASS of the simulated model.
- `toprodimo` works on top of an initialized ProDiMo model that has to be run beforehand with the typical parameters of the simulated model (disk, star, ...). The corresponding ProDiMo.out parameter is then copied to a new prodimo directory to perform the interpolation of the simulated data to this new ProDiMo model.

See also the [ProDiMo documentation](https://prodimowiki.readthedocs.io/en/latest/userguide/interface2D.html).
