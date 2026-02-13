import os
from pathlib import Path
import tarfile
import numpy as np
import shutil
import pytest

from toprodimo.main import main

def partial_extraction(*, target_tarfile:str, target_subdirectory:str, path:str):
    with tarfile.open(target_tarfile) as tar:
        subdir_and_files = [
            tarinfo for tarinfo in tar.getmembers()
            if tarinfo.name.startswith(target_subdirectory)
        ]
        tar.extractall(path=path, members=subdir_and_files, filter="data")

class TestFileWrite:
    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    def test_all_files(self, test_data_dir, data_dir):
        if os.path.isdir(test_data_dir / "idefix_1_dust_fluid"):
            shutil.rmtree(test_data_dir / "idefix_1_dust_fluid")
        partial_extraction(
            target_tarfile=data_dir / "idefix_1_dust_fluid.tar.gz", 
            target_subdirectory="idefix_1_dust_fluid/prodimo_from_ref/",
            path=test_data_dir,
        )
        partial_extraction(
            target_tarfile=data_dir / "idefix_1_dust_fluid.tar.gz", 
            target_subdirectory="idefix_1_dust_fluid/prodimo_to_ref/",
            path=test_data_dir,
        )
        partial_extraction(
            target_tarfile=data_dir / "idefix_1_dust_fluid.tar.gz", 
            target_subdirectory="idefix_1_dust_fluid/idefix_ref/",
            path=test_data_dir,
        )
        data_dir_from_ref = test_data_dir / "idefix_1_dust_fluid" / "prodimo_from_ref"
        data_dir_to_ref = test_data_dir / "idefix_1_dust_fluid" / "prodimo_to_ref"
        data_dir_to = test_data_dir / "idefix_1_dust_fluid" / "prodimo_to"

        if os.path.isdir(data_dir_to):
            shutil.rmtree(data_dir_to)

        main([str(test_data_dir / "toprodimo.toml")])

        pluto_rho_ref = np.loadtxt(os.path.join(data_dir_to_ref, "pluto_rho.dat"))
        pluto_rho = np.loadtxt(os.path.join(data_dir_to, "pluto_rho.dat"))
        np.testing.assert_array_equal(pluto_rho_ref, pluto_rho)

        pluto_tgas_ref = np.loadtxt(os.path.join(data_dir_to_ref, "pluto_tgas.dat"))
        pluto_tgas = np.loadtxt(os.path.join(data_dir_to, "pluto_tgas.dat"))
        np.testing.assert_array_equal(pluto_tgas_ref, pluto_tgas)

        pluto_vx_ref = np.loadtxt(os.path.join(data_dir_to_ref, "pluto_vx.dat"))
        pluto_vx = np.loadtxt(os.path.join(data_dir_to, "pluto_vx.dat"))
        np.testing.assert_array_equal(pluto_vx_ref, pluto_vx)

        pluto_vy_ref = np.loadtxt(os.path.join(data_dir_to_ref, "pluto_vy.dat"))
        pluto_vy = np.loadtxt(os.path.join(data_dir_to, "pluto_vy.dat"))
        np.testing.assert_array_equal(pluto_vy_ref, pluto_vy)

        pluto_vz_ref = np.loadtxt(os.path.join(data_dir_to_ref, "pluto_vz.dat"))
        pluto_vz = np.loadtxt(os.path.join(data_dir_to, "pluto_vz.dat"))
        np.testing.assert_array_equal(pluto_vz_ref, pluto_vz)

        pluto_dustsizedist_ref = np.loadtxt(os.path.join(data_dir_to_ref, "pluto_dustsizedist.dat"), skiprows=2)
        pluto_dustsizedist = np.loadtxt(os.path.join(data_dir_to, "pluto_dustsizedist.dat"), skiprows=2)
        np.testing.assert_array_equal(pluto_dustsizedist_ref, pluto_dustsizedist)

        shutil.rmtree(data_dir_to)
        shutil.rmtree(test_data_dir / "idefix_1_dust_fluid")