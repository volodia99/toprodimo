import os
import argparse
import shutil
import glob

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import astropy.units as u
import astropy.constants as uc
import tomli
import tomli_w

import inifix
from deep_chainmap import DeepChainMap
from nonos.api import GasDataSet

import prodimopy.interface2D.infile as pin2D
import prodimopy.read as pread
import prodimopy.plot as pplot
import prodimopy.plot_models as pplotm

from toprodimo._typing import F, FArray1D, FArray2D
from toprodimo._parsing import is_set, list_of_middle_keys
from toprodimo.default import DEFAULT_LAYER, MANDATORY_SET

plt.style.use("nonos.default")

def computeSizeMM(betai:FArray1D, *, internal_rho:u.Quantity, unit_length_au:u.Quantity, unit_mass_msun:u.Quantity) -> u.Quantity:
    sizeMM_k = betai/(np.sqrt(np.pi/8.0)*internal_rho.to(u.g/u.cm/u.cm/u.cm)*(unit_length_au.to(u.cm))**2/unit_mass_msun.to(u.g))
    return sizeMM_k.to(u.mm)

def iterative_mean(*, data:FArray2D[F], count:int=0, mean:FArray2D[F]|None=None) -> FArray2D[F]:
    """
    Compute the iterative mean of data arrays.
    """
    if mean is None:
        mean = np.zeros_like(data)
    mean += 1.0 / (count + 1) * (data - mean)
    return mean

def pol2cart(*, r:FArray2D[F], theta:FArray2D[F]) -> tuple[FArray2D[F], FArray2D[F]]:
    """
    Conversion of the coordinate system to cartesian.
    """
    x = r*np.sin(theta)
    z = r*np.cos(theta)
    return (x, z)

def vpol2cart(*, vr:FArray2D[F], vtheta:FArray2D[F], r:FArray2D[F], theta:FArray2D[F]) -> tuple[FArray2D[F], FArray2D[F]]:
    """
    Convert velocity components from polar to cartesian coordinates.
    """
    sint = np.sin(theta)
    cost = np.cos(theta)
    vx = vr * sint + vtheta * cost
    vz = vr * cost - vtheta * sint
    return (vx, vz)

# TODO: take care of 2D dust fluids when prodimo can handle it
def load_model(
    ds:GasDataSet, 
    unit_length_au:float, 
    unit_mass_msun:float, 
    component:list,
    config:dict,
    ):
    """
    Load the simulation model from a simulation file
    - ds (GasDataSet): nonos.api.analysis.GasDataSet
    - unit_length_au (float): typical length in au
    - unit_mass_msun (float): typical mass in solMass
    - component (list[str]): what components are included
    - config (dict): ChainMap containing all the parameters 

    Returns a prodimopy Interface2Din object.
    """
    # For converting code units to real units
    unit_length_au = unit_length_au * u.au
    unit_mass_msun = unit_mass_msun * u.M_sun
    UNIT_DENSITY = unit_mass_msun / unit_length_au**3
    UNIT_VELOCITY = np.sqrt(uc.G*unit_mass_msun/unit_length_au).to(u.m/u.s)
    # TODO: for now temperature conversion with fixed mustar
    # TODO: check if prodimo parameter how temperature is computed
    MUSTAR = 1.37
    UNIT_TEMPERATURE = ((MUSTAR*uc.m_p*uc.G/uc.k_B)*unit_mass_msun/unit_length_au).to(u.K)
    if "dust" in component:
        internal_rho = config["simulation"]["internal_rho"] * (u.g/u.cm/u.cm/u.cm)

    density = None
    velocity_r = None
    velocity_theta = None
    velocity_phi = None
    temperature = None
    dust_density = None
    dust_size_distribution = None

    if ds.native_geometry!="spherical":
        raise ValueError(f"native_geometry='{ds.native_geometry}' should be 'spherical'")
    nr, ntheta, nphi = list(ds.values())[0].data.shape
    if (nphi!=1) | (nr<2) | (ntheta<2):
        raise ValueError(f"{(nr, ntheta, nphi) = } should be (nr, nt, 1): works for 2D r-theta simulations")
    if density is None:
        r, theta, phi = (ds.coords.get_axis_array_med(kk) for kk in ("r", "theta", "phi"))
        # TODO: for now we take only the upper half if theta is symmetric compared to the midplane 
        theta_edge = ds.coords.get_axis_array("theta")
        min_half_dtheta_edge = np.diff(theta_edge).min()/2
        if not np.isclose(np.mean(theta_edge)-np.pi/2, np.float32(0.0), atol=min_half_dtheta_edge):
            raise NotImplementedError("Should be symmetric compared to the midplane (pi/2), in order to extract upper half")
        theta = theta[0:ntheta//2+1]
    if "gas" in component:
        density = iterative_mean(
            data=ds["RHO"].data[:,0:ntheta//2+1,0], 
            count=0, 
            mean=density
        )
        density = ((density * UNIT_DENSITY).to(u.g / u.cm**3)).value

        velocity_r = iterative_mean(
            data=ds["VX1"].data[:,0:ntheta//2+1,0], 
            count=0, 
            mean=velocity_r
        )
        velocity_theta = iterative_mean(
            data=ds["VX2"].data[:,0:ntheta//2+1,0], 
            count=0, 
            mean=velocity_theta
        )
        velocity_phi = iterative_mean(
            data=ds["VX3"].data[:,0:ntheta//2+1,0], 
            count=0, 
            mean=velocity_phi
        )
        velocity = ((np.dstack((velocity_r, velocity_theta, velocity_phi))*UNIT_VELOCITY).to(u.cm / u.s)).value

        temperature = iterative_mean(
            data=ds["PRS"].data[:,0:ntheta//2+1,0]/ds["RHO"].data[:,0:ntheta//2+1,0], 
            count=0, 
            mean=temperature
        )
        temperature = (temperature * UNIT_TEMPERATURE).value

    if "dust" in component:
        print(f"WARN: 'dust' not implemented in a general way with nonos. Implementation specific to IDEFIX.")
        directory = ds._parameters_input["directory"]
        inifile = inifix.load(os.path.join(directory, "idefix.ini"))
        dragType = inifile["Dust"]["drag"][0]
        if dragType!="size":
            raise ValueError(f"{dragType=} should be 'size'.")
        dustBeta = np.array(inifile["Dust"]["drag"][1:])
        dustSize_cm = computeSizeMM(
            dustBeta, 
            internal_rho=internal_rho,
            unit_length_au=unit_length_au,
            unit_mass_msun=unit_mass_msun,
            ).to(u.cm)
        dustSize_cm_ascending_indices = np.argsort(dustSize_cm)
        dust_density = np.empty((len(dustSize_cm),)+density.shape)
        for kk in dustSize_cm_ascending_indices:
            dust_density[kk, ...] = iterative_mean(
                    data=ds[f"DUST{kk}_RHO"].data[:,0:ntheta//2+1,0], 
                    count=0, 
                    mean=None,
                )

    rr, tt = np.meshgrid(r, theta, indexing="ij")
    xx, zz = pol2cart(
        r=rr, 
        theta=tt
    )
    # flip so that z=0 has zidx=0
    xx = np.flip(xx, -1)
    zz = np.flip(zz, -1)

    # z can be < 0 in the midplane: set it to zero
    zz[zz[:, 0] < 0, 0] = 0.0

    if "gas" in component:
        # velocity is in vr,vtheta,vphi
        # get vx and vz from vr and vtheta
        vy = velocity[:, :, 2]
        vx, vz = vpol2cart(
            vr=velocity[:, :, 0], 
            vtheta=velocity[:, :, 1], 
            r=rr, 
            theta=tt
        )

        velocity[:, :, 0] = vx
        velocity[:, :, 1] = vy
        velocity[:, :, 2] = vz

        # flip so that z=0 has zidx=0
        density = np.flip(density, -1)
        velocity = np.flip(velocity, 1)
        temperature = np.flip(temperature, -1)
    if "dust" in component:
        dust_density = np.flip(dust_density, -1)
        dust_size_distribution = pin2D.DustSizeDistribution(
            asize=dustSize_cm.value,
            fsize_rho=((dust_density * UNIT_DENSITY).to(u.g / u.cm**3)).value,
        )

    # Use the prodimopy tools to generate an object for further processing
    return pin2D.Interface2Din(
        ((xx * unit_length_au).to(u.cm)).value,
        ((zz * unit_length_au).to(u.cm)).value,
        rhoGas=density,
        velocity=velocity,
        tgas=temperature,
        dustSDF=dust_size_distribution,
    )

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="toprodimo",
        description=__doc__,
    )
    parser.suggest_on_error = True  # type: ignore [attr-defined]

    parser.add_argument(
        "parameter_file", 
        type=str, 
        help="work on toprodimo parameter file"
    )

    return parser

def plot_model(model, pdf_name:str=""):
    """
    Plot various quantities of the model using prodimopy plotting tools.
    """
    xmin = model.x.min()
    xmax = model.x.max()
    ymax = model.z.max()
    rhogmin = np.nanmin(model.rhog[np.nonzero(model.rhog)])
    rhogmax = np.nanmax(model.rhog)

    xlim = [xmin, xmax]
    ylim = [0, ymax]

    constyle = {"zr": False, "xlog": False, "ylog": False, "axequal": True}

    with PdfPages(f"{pdf_name}.pdf") as pdf:
        pp = pplot.Plot(pdf)
        pp.plot_grid(model, xlim=xlim, ylim=ylim, **constyle)

        pp.plot_cont(
            model,
            "rhog",
            **constyle,
            zlim=[rhogmin, rhogmax],
            extend="both",
            xlim=xlim,
            ylim=ylim,
        )
        pp.plot_cont(
            model,
            "rhog",
            **constyle,
            zlim=[rhogmin, rhogmax],
            extend="both",
            xlim=[0.75*xmin, 3*xmin],
            ylim=[0, 3*xmin],
        )

        pp.plot_cont(model, "tg", **constyle, xlim=xlim, ylim=ylim)
        pp.plot_cont(model, "tg", **constyle, xlim=[0.75*xmin, 3*xmin], ylim=[0, 3*xmin])

        for i, label in enumerate(["vx", "vy", "vz"]):
            if label=="vy":
                zlim = [0, np.nanmax(model.velocity[:,2])]
                cmap = "plasma"
            else:
                zlim = [-0.01*np.nanmax(model.velocity[:,2]), 0.01*np.nanmax(model.velocity[:,2])]
                cmap = "berlin"
            pp.plot_cont(
                model,
                model.velocity[:, :, i],
                label,
                **constyle,
                zlog=False,
                extend="both",
                zlim=zlim,
                cmap=cmap,
                xlim=xlim,
                ylim=ylim,
            )
            pp.plot_cont(
                model,
                model.velocity[:, :, i],
                label,
                **constyle,
                zlog=False,
                extend="both",
                zlim=zlim,
                xlim=[0.75*xmin, 3*xmin],
                ylim=[0, 3*xmin],
                cmap=cmap,
            )

        if os.path.basename(pdf_name)=="prodimo":
            rhodmin = np.nanmin(model.rhod[np.nonzero(model.rhod)])
            rhodmax = np.nanmax(model.rhod)
            pp.plot_cont(
                model,
                "rhod",
                **constyle,
                zlim=[rhogmin, rhogmax],
                extend="both",
                xlim=xlim,
                ylim=ylim,
            )
            pp.plot_cont(
                model,
                "rhod",
                **constyle,
                zlim=[rhogmin, rhogmax],
                extend="both",
                xlim=[0.75*xmin, 3*xmin],
                ylim=[0, 3*xmin],
            )

        # Check the Keplerian velocity
        pp.plot_radial(
            model,
            model.velocity[:, :, 1],
            "vy",
            zidx=0,
            ylog=False,
            xlog=True,
            ylim=[np.nanmin(model.velocity[:, 0, 1]), np.nanmax(model.velocity[:, 0, 1])],
            title="Keplerian Velocity Vy at z=0 au",
        )

def main(argv: list[str] | None = None) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)
    with open(args.parameter_file, "rb") as f:
        config_file_layer = tomli.load(f)

    print("")

    expected_length_mandatory = 7
    if len(MANDATORY_SET)!=expected_length_mandatory:
        raise ValueError(
            f"expecting {expected_length_mandatory} mandatory parameters, not {len(MANDATORY_SET)}."
        )
    if "dust" in config_file_layer["simulation"]["component"]:
        MANDATORY_SET.update(["internal_rho"])
    if MANDATORY_SET.difference(set(list_of_middle_keys(config_file_layer))):
        raise ValueError(
            f"at least one mandatory parameter is missing: {MANDATORY_SET.difference(set(list_of_middle_keys(config_file_layer)))}"
        )
    if "dust" not in config_file_layer["simulation"]["component"] and is_set(config_file_layer["simulation"]["internal_rho"]):
        print("WARN: unused 'internal_rho' when 'dust' not included.")
        config_file_layer["simulation"]["internal_rho"] = "unset"

    config = DeepChainMap(config_file_layer, DEFAULT_LAYER)
    if not is_set(config["simulation"]["mask_inside"]):
        config["simulation"]["mask_inside"] = 0.0

    component = config["simulation"]["component"]
    component = np.atleast_1d(component).tolist()
    if not (set(component) & set(["gas","dust"])):
        raise ValueError(
            f"{component=} should be 'dust', 'gas' or ['dust', 'gas']."
        )

    if not all((config["prodimo"]["from"], config["prodimo"]["to"])):
        raise ValueError(
            f"init_prodimo_model_directory={config["prodimo"]["from"]}, prodimo_model_directory={config["prodimo"]["to"]}. "\
            "Both the init_prodimo_model_directory and the prodimo_model_directory have to be specified in .toml file, "\
            "using from = 'init_prodimo_model_directory' and to = 'prodimo_model_directory'"\
        )

    if (not(os.path.isdir(config["prodimo"]["from"]))) | (not(os.path.exists(os.path.join(os.getcwd(), config["prodimo"]["from"], "ProDiMo.out")))):
        raise ValueError(
            f"{os.path.join(config["prodimo"]["from"], 'ProDiMo.out')} must exist. "\
            "In order to do the conversion from simulation outputs to ProDiMo, an init ProDiMo model has to be run beforehand. "\
            "The resulting proper grid and parameters will then be used by ProDiMo. "\
            "One needs to run it at least once with stop_after_init=.true. in order to create the ProDiMo.out file."
        )

    file = config["simulation"]["on"]
    input_dir = config["simulation"]["input_dir"]
    ds = GasDataSet(file, directory=input_dir)

    # load the model
    model = load_model(
        ds=ds,
        unit_length_au=config["simulation"]["unit_length_au"], 
        unit_mass_msun=config["simulation"]["unit_mass_msun"], 
        component=component,
        config=config,
    )

    # Make some manipulations
    # remove everything inside the inner boundary of the model
    xxnew = (model.x * u.cm).to(u.au).value
    zznew = (model.z * u.cm).to(u.au).value
    rrnew = np.sqrt(xxnew**2 + zznew**2)
    rincut = rrnew.min()*0.95  # use this as the cut ...

    mask_inner_edge = np.zeros_like(xxnew, dtype=bool)
    mask_inner_edge = (xxnew < rincut)

    if config["simulation"]["mask_inside"]:
        print(f"INFO: canceling (vx, vz) inside {config["simulation"]["mask_inside"]:.2f} r_inner...\n")
        model.velocity[(xxnew < rincut * config["simulation"]["mask_inside"]), 0] = 0.0
        model.velocity[(xxnew < rincut * config["simulation"]["mask_inside"]), 2] = 0.0

    model.rhoGas[mask_inner_edge] = np.nan
    model.tgas[mask_inner_edge] = np.nan
    for i in range(3):
        model.velocity[mask_inner_edge, i] = np.nan

    if not(os.path.isdir(config["prodimo"]["to"])):
        print(f"WARN: '{config["prodimo"]["to"]}' must exist. Creating it...\n")
        os.makedirs(config["prodimo"]["to"])

    if os.path.isfile(os.path.join(config["prodimo"]["to"], "ProDiMo.out")):
        raise ValueError(
            f"{os.path.join(config["prodimo"]["to"], 'ProDiMo.out')} already exists. "\
            "Check if you are sure of what you are doing. If you do, delete the preexisting 'ProDiMo.out' file yourself."
        )
    shutil.copy2(os.path.join(config["prodimo"]["from"], "ProDiMo.out"), config["prodimo"]["to"])
    for file in glob.glob(os.path.join(config["prodimo"]["from"], "*.in")):
        shutil.copy2(file, config["prodimo"]["to"])

    prodimo_model = pread.read_prodimo(config["prodimo"]["to"], name="ProDiMo model")

    # The new 2D input files are written to the prodimo_model directory
    interpolated_prodimo_model = model.toProDiMo(prodimo_model, outdir=prodimo_model.directory, fixmethod=0)

    if config["prodimo"]["plot"]:
        pmodel = model.get_pmodel()
        # Plot directly the data from the simulation
        plot_model(pmodel, pdf_name=os.path.join(config["prodimo"]["to"], "simulation"))

        # Plot the new ProDiMo model
        plot_model(interpolated_prodimo_model, pdf_name=os.path.join(config["prodimo"]["to"], "prodimo"))

        models = [pmodel, interpolated_prodimo_model]

        with PdfPages(os.path.join(config["prodimo"]["to"], "compare_simulation_prodimo.pdf")) as pdf:
            ppm = pplotm.PlotModels(pdf, styles=["-", "--"])
            # We can directly compare the new intepolated numbers to the original MDH model
            ppm.plot_midplane(models, "rhog", "rhog", ylim=[np.nanmin(models[0].rhog[np.nonzero(models[0].rhog)]), np.nanmax(models[0].rhog)])#[1e-24, 1e-10])
            # This plot shows some deviations, the reason is that we have to few vertical points, and most of them are concentrated towards the midplane
            # where the interpolation is still accurate
            ppm.plot_vertical(models, config["simulation"]["unit_length_au"], "rhog", "rhog", xlim=[1, 0])

    # Writing config to a TOML file
    with open(os.path.join(config["prodimo"]["to"], "toprodimo.full.toml"), "wb") as outfile:
        tomli_w.dump(config, outfile)

    return 0