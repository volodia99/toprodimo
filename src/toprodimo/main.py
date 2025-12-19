import os
import argparse
import shutil
import json

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import astropy.units as u
import astropy.constants as uc

from nonos.api import GasDataSet

import prodimopy.interface2D.infile as pin2D
import prodimopy.read as pread
import prodimopy.plot as pplot
import prodimopy.plot_models as pplotm

from toprodimo._typing import F, FArray2D

plt.style.use("nonos.default")

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
def load_model(file:str|int, *, directory:None|str=None, UNIT_LENGTH:None|float=None, UNIT_MASS:None|float=None):
    """
    Load the simulation model from a simulation file
    - file (str|int): absolute path of the simulation output file, or output number
    - directory (str): location of the simulated output file, if described by its output number
    - UNIT_LENGTH (float): typical length in au
    - UNIT_MASS (float): typical mass in solMass

    Returns a prodimopy Interface2Din object.
    """

    if isinstance(file, str):
        if directory is not None:
            raise ValueError(
                f"if {file=} is an absolute path, {directory=} should not be defined"
            )
        ds = GasDataSet(file)
    elif isinstance(file, int):
        if directory is None:
            raise ValueError(
                f"if {file=} is an integer output number, {directory=} should be defined"
            )
        ds = GasDataSet(file, directory=directory)

    # For converting code units to real units
    UNIT_LENGTH = UNIT_LENGTH * u.au
    UNIT_MASS = UNIT_MASS * u.M_sun
    UNIT_DENSITY = UNIT_MASS / UNIT_LENGTH**3
    UNIT_VELOCITY = np.sqrt(uc.G*UNIT_MASS/UNIT_LENGTH).to(u.m/u.s)
    # TODO: for now temperature conversion with fixed mustar
    # TODO: check if prodimo parameter how temperature is computed
    MUSTAR = 1.37
    UNIT_TEMPERATURE = ((MUSTAR*uc.m_p*uc.G/uc.k_B)*UNIT_MASS/UNIT_LENGTH).to(u.K)

    density = None
    velocity_r = None
    velocity_theta = None
    velocity_phi = None
    temperature = None

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
        if not np.isclose(np.ptp(theta_edge)-np.pi/2, np.float32(0.0), atol=min_half_dtheta_edge):
            raise NotImplementedError("Should be symmetric compared to the midplane (pi/2), in order to extract upper half")
        theta = theta[0:ntheta//2+1]
    density = iterative_mean(
        data=ds["RHO"].data[:,0:ntheta//2+1,0], 
        count=0, 
        mean=density
    )
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
    temperature = iterative_mean(
        data=ds["PRS"].data[:,0:ntheta//2+1,0]/ds["RHO"].data[:,0:ntheta//2+1,0], 
        count=0, 
        mean=temperature
    )

    rr, tt = np.meshgrid(r, theta, indexing="ij")
    xx, zz = pol2cart(
        r=rr, 
        theta=tt
    )

    velocity = np.dstack(
        (
            velocity_r,
            velocity_theta,
            velocity_phi,
        )
    )
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
    xx = np.flip(xx, 1)
    zz = np.flip(zz, 1)
    density = np.flip(density, 1)
    velocity = np.flip(velocity, 1)
    temperature = np.flip(temperature, 1)

    # z can be < 0 in the midplane: set it to zero
    zz[zz[:, 0] < 0, 0] = 0.0

    # Use the prodimopy tools to generate an object for further processing
    return pin2D.Interface2Din(
        ((xx * UNIT_LENGTH).to(u.cm)).value,
        ((zz * UNIT_LENGTH).to(u.cm)).value,
        rhoGas=((density * UNIT_DENSITY).to(u.g / u.cm**3)).value,
        velocity=((velocity * UNIT_VELOCITY).to(u.cm / u.s)).value,
        tgas=(temperature * UNIT_TEMPERATURE).value,
    )

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="toprodimo",
        description=__doc__,
    )
    parser.suggest_on_error = True  # type: ignore [attr-defined]

    subparsers = parser.add_subparsers(
        help="Choice between from_on (from output number) and from_path (from absolute path to file name)", 
        dest="input_processing"
    )

    ###
    ## Create the parser for the "on" command
    parser_on = subparsers.add_parser(
        "from_on", 
        help="Work with output number approach. Try 'toprodimo from_on -h' for more info."
    )
    parser_on.add_argument(
        "file_on", 
        type=int, 
        help="Get simulation output from its output number"
    )

    parser_on.add_argument(
        "-dir",
        type=str,
        required=True,
        dest="directory",
        help="required: location of the simulated output, if described by its output number.",
    )

    ## Create the parser for the "name" command
    parser_file = subparsers.add_parser(
        "from_path", 
        help="Work with file name approach. Try 'toprodimo from_path -h' for more info."
    )
    parser_file.add_argument(
        "file_path", 
        type=str, 
        help="Get simulation output from its name"
    )

    for subparser in [parser_on, parser_file]:
        subparser.add_argument(
            "-UNIT_LENGTH",
            type=float,
            default=None,
            required=True,
            help="required: code unit of length [au].",
        )

        subparser.add_argument(
            "-UNIT_MASS",
            type=float,
            default=None,
            required=True,
            help="required: code unit of mass [solMass].",
        )

        subparser.add_argument(
            "-mask_inside",
            type=float,
            default=1.2,
            help="mask the velocities inside given radius [inner edge unit]. put 0 for no masking. (default: 1.2).",
        )

        subparser.add_argument(
            "-from_pmdir",
            "-from",
            type=str,
            default=None,
            dest="init_prodimo_model_directory",
            help="location of the initialized prodimo model from which to extract ProDiMo.out. Works only with -to_pmdir.",
        )

        subparser.add_argument(
            "-to_pmdir",
            "-to",
            type=str,
            default=None,
            dest="prodimo_model_directory",
            help="location of the prodimo model on which the simulation model is then interpolated. Works only with -from_pmdir.",
        )

        subparser.add_argument(
            "-plot",
            "-p",
            action="store_true",
            help="rough plotting procedure to check the validity of toprodimo results.",
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

    print(rhogmin, rhogmax)

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
    print("")

    if not all((args.init_prodimo_model_directory, args.prodimo_model_directory)):
        raise ValueError(
            f"init_prodimo_model_directory={args.init_prodimo_model_directory}, prodimo_model_directory={args.prodimo_model_directory}. "\
            "Both the init_prodimo_model_directory and the prodimo_model_directory have to be specified, "\
            "using: -from 'init_prodimo_model_directory' -to 'prodimo_model_directory'"\
        )

    if (not(os.path.isdir(args.init_prodimo_model_directory))) | (not(os.path.exists(os.path.join(os.getcwd(), args.init_prodimo_model_directory, "ProDiMo.out")))):
        raise ValueError(
            f"{os.path.join(args.init_prodimo_model_directory, 'ProDiMo.out')} must exist. "\
            "In order to do the conversion from simulation outputs to ProDiMo, an init ProDiMo model has to be run beforehand. "\
            "The resulting proper grid and parameters will then be used by ProDiMo. "\
            "One needs to run it at least once with stop_after_init=.true. in order to create the ProDiMo.out file."
        )

    if args.input_processing=="from_path":
        file = args.file_path
        directory = None
        if not(os.path.isfile(file)):
            raise FileNotFoundError(f"absolute path of the file '{file}' must exist.")
    elif args.input_processing=="from_on":
        file = args.file_on
        directory = args.directory

    # load the model
    model = load_model(
        file=file,
        directory=directory, 
        UNIT_LENGTH=args.UNIT_LENGTH, 
        UNIT_MASS=args.UNIT_MASS, 
    )

    # Make some manipulations
    # remove everything inside the inner boundary of the model
    xxnew = (model.x * u.cm).to(u.au).value
    zznew = (model.z * u.cm).to(u.au).value
    rrnew = np.sqrt(xxnew**2 + zznew**2)
    rincut = rrnew.min()*0.95  # use this as the cut ...

    vzcart = model.velocity[:, :, 2]
    mask_inner_edge = np.zeros_like(xxnew, dtype=bool)
    mask_inner_edge = (xxnew < rincut)

    if args.mask_inside:
        print(f"INFO: canceling (vx, vz) inside {args.mask_inside:.2f} r_inner...\n")
        model.velocity[(xxnew < rincut * args.mask_inside), 0] = 0.0
        model.velocity[(xxnew < rincut * args.mask_inside), 2] = 0.0

    model.rhoGas[mask_inner_edge] = np.nan
    model.tgas[mask_inner_edge] = np.nan
    for i in range(3):
        model.velocity[mask_inner_edge, i] = np.nan

    if not(os.path.isdir(args.prodimo_model_directory)):
        print(f"WARN: '{args.prodimo_model_directory}' must exist. Creating it...\n")
        os.makedirs(args.prodimo_model_directory)

    if os.path.isfile(os.path.join(args.prodimo_model_directory, "ProDiMo.out")):
        raise ValueError(
            f"{os.path.join(args.prodimo_model_directory, 'ProDiMo.out')} already exists. "\
            "Check if you are sure of what you are doing. If you do, delete the preexisting 'ProDiMo.out' file yourself."
        )
    shutil.copy2(os.path.join(args.init_prodimo_model_directory, "ProDiMo.out"), args.prodimo_model_directory)

    prodimo_model = pread.read_prodimo(args.prodimo_model_directory, name="ProDiMo model")

    # The new 2D input files are written to the prodimo_model directory
    interpolated_prodimo_model = model.toProDiMo(prodimo_model, outdir=prodimo_model.directory, fixmethod=0)

    if args.plot:
        pmodel = model.get_pmodel()
        # Plot directly the data from the simulation
        plot_model(pmodel, pdf_name=os.path.join(args.prodimo_model_directory, "simulation"))

        # Plot the new ProDiMo model
        plot_model(interpolated_prodimo_model, pdf_name=os.path.join(args.prodimo_model_directory, "prodimo"))

        models = [pmodel, interpolated_prodimo_model]

        with PdfPages(os.path.join(args.prodimo_model_directory, "compare_simulation_prodimo.pdf")) as pdf:
            ppm = pplotm.PlotModels(pdf, styles=["-", "--"])
            # We can directly compare the new intepolated numbers to the original MDH model
            ppm.plot_midplane(models, "rhog", "rhog", ylim=[np.nanmin(models[0].rhog[np.nonzero(models[0].rhog)]), np.nanmax(models[0].rhog)])#[1e-24, 1e-10])
            # This plot shows some deviations, the reason is that we have to few vertical points, and most of them are concentrated towards the midplane
            # where the interpolation is still accurate
            ppm.plot_vertical(models, args.UNIT_LENGTH, "rhog", "rhog", xlim=[1, 0])

    dargs = vars(args)
    # Writing CLI args to a JSON file
    with open(os.path.join(args.prodimo_model_directory, "toprodimo_args.json"), "w") as outfile:
        json.dump(dargs, outfile, indent=1)

    return 0