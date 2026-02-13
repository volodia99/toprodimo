# toprodimo
[![PyPI](https://img.shields.io/pypi/v/toprodimo.svg?logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/toprodimo/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)

***Interface to interpolate and convert 2D (r, theta) datasets from simulation outputs (readable with [nonos](https://github.com/la-niche/nonos)) to ProDiMo models. Based on the [ProDiMo 2D interface notebook](https://prodimopy.readthedocs.io/en/stable/notebooks/interface2D.html) by Christian Rab***

## Development status

*Word of caution*: `toprodimo` has still to be tested, in particular to be sure that everything works smoothly when the ProDiMo model created by `toprodimo` is used in a ProDiMo computation.

To be implemented : 
- For now, `toprodimo` takes only the upper half of a 2D (`r`, `theta`) disk, if `theta` is symmetric compared to the midplane. Could add some flexibility if the user wants to focus on the lower half of the disk, or to average azimuthally a 3D (`r`, `theta`, `phi`) disk, or even add a prescription to extend vertically a 2D (`r`, `phi`) disk (e.g., with vertical hydrostatic equilibrium).  

## Installation

We recommend to install `toprodimo` using the package and project manager `uv`. See the [documentation](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer) to install `uv` on your system. Run the following:

```shell
uv tool install toprodimo
```

## Use the interface

You can use the interface inside the project's virtual environment using a parameter .toml file: 

```shell
toprodimo toprodimo.toml
```

See the [TOML documentation](https://toml.io/en) to know more about this config file format.

## Configuration file

### Example

You can find an example for the parameter file in `toprodimo/toprodimo.toml`.

### 1. Section `[simulation]`

***Mandatory parameters:***
- `on` : simulation output number (`int`)
- `input_dir` : directory of the simulated output (`str`)
- `unit_length_au` : code unit of length \[au\] (`float`)
- `unit_mass_msun` : code unit of mass \[solMass\] (`float`)
- `component` : which component is included (`"dust"` and/or `"gas"`) (`str|list[str]`)
- `internal_rho` : if the `"dust"` component is included, the internal density used in the simulation \[g/cm3\] (`float`)

***Optional parameters:***
- `mask_inside` : removes the contribution of the radial and vertical velocities close to the grid's inner edge `r_inner`, at `mask_inside*r_inner`, to avoid some spurious effects in ProDiMo (`float`). By default we do not cancel these velocity components (`mask_inside=0.0`).
- `tgas` : compute the gas temperature from the simulations. Implemented: `tgas = {eos="ideal", mu_star=1.37}`, using an ideal equation of state (tgas = pressure/density) and a user-defined mean molecular weight. For the moment, we recommand to let ProDiMo handle the gas temperature, by not providing any tgas.

### 2. Section `[prodimo]`

***Mandatory parameters:***
- `from` : directory of the initialized ProDiMo model from which to extract ProDiMo.out and *.in files (`str`)
- `to` : directory of the ProDiMo model on which the simulation grid and fields are interpolated (`str`)

***Optional parameters:***
- `plot` : creates 3 .pdf files in the ProDiMo model directory (given by the `to` parameter), with plots of the fields (`bool`).
    - `simulation.pdf`: from the simulation output file, with some post-processing (e.g., removing for all the fields the region inside the cylindrical radius corresponding to the inner edge).
    - `prodimo.pdf`: from the ProDiMo model, ready to be run with ProDiMo.
    - `compare_simulation_prodimo.pdf`: look at the 1D density in the midplane and vertically at R=`unit_length_au`. 

## Remarks

In order for the procedure to work, you need to keep in mind that:
- `toprodimo` needs the typical `unit_length_au` and `unit_mass_msun` of the simulated model.
- `toprodimo` works on top of an initialized ProDiMo model that has to be run beforehand with the typical parameters of the simulated model (disk, star, ...). The corresponding (ProDiMo.out, *.in) files are then copied to a new prodimo directory to perform the interpolation of the simulated data to this new ProDiMo model.

See also the [ProDiMo documentation](https://prodimowiki.readthedocs.io/en/latest/userguide/interface2D.html).






