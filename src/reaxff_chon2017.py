import numpy as np
import subprocess
import shlex
import lammps_logfile
from copy import deepcopy
from typing import List
import os
import shutil


# constants
_HARTREE_2_KCAL = 627.509
_ANGSTROM_2_BOHR = 1.8897259886


# data template
_DATA_TEMPLATE = """# chon2017

NUM_ATOMS atoms
4 atom types

 -100.000000000  100.000000000    xlo    xhi
 -100.000000000  100.000000000    ylo    yhi
 -100.000000000  100.000000000    zlo    zhi

Masses

       1   12.011000000  # C
       2    1.008000000  # H
       3   15.999000000  # O
       4   14.007000000  # N

Atoms

"""

_ATOM_MAP = {
    "C": 1,
    "H": 2,
    "O": 3,
    "N": 4,
}


class CHON2017:
    def __init__(self) -> None:
        pass

    def calc_energy(
        self,
        symbols: List[str],
        coords: np.ndarray,
    ) -> float:
        self._exec(symbols, coords)
        log = lammps_logfile.File("log.lammps")

        return log.data_dict["PotEng"][0] / _HARTREE_2_KCAL

    def calc_gradients(
        self,
        symbols: List[str],
        coords: np.ndarray,
    ) -> np.ndarray:
        forces = []
        with open("dump.force") as f:
            while line := f.readline():
                if "ITEM: ATOMS" in line:
                    for _ in range(len(coords)):
                        line = f.readline()
                        forces.append([float(x) for x in line.split()])

        # (kcal/mol) / Angstrom -> Hartree / Bohr
        return -1.0 * np.array(forces) / _HARTREE_2_KCAL / _ANGSTROM_2_BOHR

    def _exec(
        self,
        symbols: List[str],
        coords: np.ndarray,
    ) -> None:
        # append atoms to data
        data_contants = deepcopy(_DATA_TEMPLATE)
        data_contants = data_contants.replace("NUM_ATOMS", f"{len(symbols)}")

        for (atom_idx, (atom_symbol, atom_xyz)) in enumerate(zip(symbols, coords)):
            data_contants += \
                f"{(atom_idx + 1):8d}{_ATOM_MAP[atom_symbol]:8d}{0.0:15.9f}" \
                f"{atom_xyz[0]:15.9f}{atom_xyz[1]:15.9f}{atom_xyz[2]:15.9f}\n"

        with open("data.gau", "w") as f:
            print(data_contants, file=f)

        # run lammps
        subprocess.run(
            args=shlex.split("lmp -in in.force"),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )



if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver.from_stdio()

    if driver.derivs == 2:
        raise NotImplementedError(
            "ReaxFF calculator only supports energy and gradients."
        )

    pes = CHON2017()

    driver.write(
        energy=pes.calc_energy(driver.symbols, driver.coords),
        gradients=(pes.calc_gradients(driver.symbols, driver.coords)
                   if driver.derivs == 1 else None),
        force_constants=None,
    )
