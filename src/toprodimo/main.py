import os
import argparse
import shutil

import numpy as np
import astropy.units as u
import astropy.constants as uc

from nonos.api import GasDataSet

import prodimopy.interface2D.infile as pin2D
import prodimopy.read as pread

def iterative_mean(*, data, count:int=0, mean=None):
    """
    Compute the iterative mean of data arrays.
    """
    if mean is None:
        mean = np.zeros(data.shape)
    mean += 1.0 / (count + 1) * (data - mean)

    return mean

def pol2cart(*, r, theta):
    """
    Conversion of the coordinate system to cartesian.
    """
    x = r*np.sin(theta)
    z = r*np.cos(theta)
    return (x, z)

def vpol2cart(*, vr, vtheta, r, theta):
    """
    Convert velocity components from polar to cartesian coordinates.
    """
    sint = np.sin(theta)
    cost = np.cos(theta)
    vx = vr * sint + vtheta * cost
    vz = vr * cost - vtheta * sint
    return (vx, vz)

# TODO: make iterative_mean optional when called
# For now compute time average of all .vtk contained in a given directory
def load_model(filename:str, *, directory:str=".", UNIT_LENGTH:None|float=None, UNIT_MASS:None|float=None):
    """
    Load the simulation model from a directory containing multiple files
    - filename (str): name of the simulation output file
    - directory (str): location where to find 'filename'
    - UNIT_LENGTH (float): typical length in au
    - UNIT_MASS (float): typical mass in solMass

    Returns a prodimopy Interface2Din object.
    """

    if UNIT_LENGTH is None:
        raise ValueError("UNIT_LENGTH should be specified with corresponding command-line argument.")
    if UNIT_MASS is None:
        raise ValueError("UNIT_MASS should be specified with corresponding command-line argument.")
    # For converting code units to real units
    UNIT_LENGTH = UNIT_LENGTH * u.au
    UNIT_MASS = UNIT_MASS * u.M_sun
    UNIT_DENSITY = UNIT_MASS / UNIT_LENGTH**3
    UNIT_VELOCITY = np.sqrt(uc.G*UNIT_MASS/UNIT_LENGTH).to(u.m/u.s)
    # TODO: for now temperature conversion with fixed mustar
    # TODO: check if prodimo parameter how temperature is computed
    MUSTAR = 1.37
    UNIT_TEMPERATURE = ((MUSTAR*uc.m_p*uc.G/uc.k_B)*UNIT_MASS/UNIT_LENGTH).to(u.K)

    if directory.startswith("~"):
        directory = os.path.expanduser(directory)

    density = None
    velocity_r = None
    velocity_theta = None
    velocity_phi = None
    temperature = None
    ds = GasDataSet(os.path.join(directory, filename))
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

    parser.add_argument(
        "filename",
        type=str,
        dest="filename",
        help="name of the file location of output files and param files.",
    )

    parser.add_argument(
        "-dir",
        type=str,
        default=".",
        dest="directory",
        help="location of output files and param files (default: '.').",
    )

    parser.add_argument(
        "-UNIT_LENGTH",
        type=float,
        default=None,
        dest="UNIT_LENGTH",
        help="code unit of length [au].",
    )

    parser.add_argument(
        "-UNIT_MASS",
        type=float,
        default=None,
        dest="UNIT_MASS",
        help="code unit of mass [solMass].",
    )

    parser.add_argument(
        "-clean",
        type=float,
        default=None,
        dest="clean",
        help="mask the velocities inside given cleaning radius [inner edge unit] (default: None).",
    )

    parser.add_argument(
        "-from_pmdir",
        "-from",
        type=str,
        default=None,
        dest="init_prodimo_model_directory",
        help="location of the initialized prodimo model from which to extract ProDiMo.out. Works with -to_pmdir.",
    )

    parser.add_argument(
        "-to_pmdir",
        "-to",
        type=str,
        default=None,
        dest="prodimo_model_directory",
        help="location of the prodimo model on which the simulation model is then interpolated. Works with -from_pmdir.",
    )

    # parser.add_argument(
    #     "-prodimo_model_dir",
    #     "-pmdir",
    #     type=str,
    #     default=None,
    #     dest="prodimo_model_directory",
    #     help="location of initialized prodimo model on which the simulation model is interpolated.",
    # )

def main(argv: list[str] | None = None) -> int:
    parser = get_parser()
    args = vars(parser)

    # load the model
    model = load_model(
        filename=args.filename,
        directory=args.directory, 
        UNIT_LENGTH=args.UNIT_LENGTH, 
        UNIT_MASS=args.UNIT_MASS, 
    )

    # Make some manipulations
    # remove everything inside the inner boundary of the model
    xxnew = (model.x * u.cm).to(u.au).value
    zznew = (model.z * u.cm).to(u.au).value
    rrnew = np.sqrt(xxnew**2 + zznew**2)
    rincut = rrnew.min()  # use this as the cut ...

    vzcart = model.velocity[:, :, 2]
    mask_inner_edge = np.zeros_like(xxnew, dtype=bool)
    mask_inner_edge = (xxnew < rincut)

    if args.clean:
        model.velocity[(xxnew < rincut * args.clean), 0] = 0.0
        model.velocity[(xxnew < rincut * args.clean), 2] = 0.0

    model.rhoGas[mask_inner_edge] = np.nan
    model.tgas[mask_inner_edge] = np.nan
    for i in range(3):
        model.velocity[mask_inner_edge, i] = np.nan

    if not all(args.init_prodimo_model_directory, args.prodimo_model_directory):
        raise ValueError(
            f"init_prodimo_model_directory={args.init_prodimo_model_directory}, prodimo_model_directory={args.prodimo_model_directory}"\
            "Both the init_prodimo_model_directory and the prodimo_model_directory have be specified,"\
            "using: -from 'init_prodimo_model_directory' -to 'prodimo_model_directory'"\
        )

    if not(os.path.isdir(args.init_prodimo_model_directory)) | not(os.path.exists(os.path.join(os.getcwd(), args.init_prodimo_model_directory, "ProDiMo.out"))):
        raise ValueError(
            f"Both '{args.init_prodimo_model_directory}' and {os.path.join(args.init_prodimo_model_directory, 'ProDiMo.out')} must exist."\
            "In order to do the conversion from simulation outputs to ProDiMo, an init ProDiMo model has to be run beforehand."\
            "The resulting proper grid and parameters will then be used by ProDiMo."\
            "One needs to run it at least once with stop_after_init=.true. in order for the toProDiMo routine to read the ProDiMo.out file."
        )

    if not(os.path.isdir(args.prodimo_model_directory)):
        print(f"'{args.prodimo_model_directory}' must exist. Creating it...")
        os.makedirs(args.prodimo_model_directory)

    if os.path.exists(os.path.join(os.getcwd(), args.prodimo_model_directory, "ProDiMo.out")):
        raise ValueError(
            f"{os.path.join(args.init_prodimo_model_directory, 'ProDiMo.out')} already exists."\
            "Check if you are sure of what you are doing. If you do, delete the preexisting 'ProDiMo.out' file yourself."
        )
    shutil.copy2(os.path.join(args.init_prodimo_model_directory, "ProDiMo.out"), args.prodimo_model_directory)

    prodimo_model = pread.read_prodimo(args.prodimo_model_directory, name="ProDiMo model")

    # The new 2D input files are written to the init_prodimo_model directory
    interpolated_prodimo_model = model.toProDiMo(init_prodimo_model, outdir=init_prodimo_model.directory, fixmethod=0)

    return 0