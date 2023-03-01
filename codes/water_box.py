"""make water box with packmmol
Preparing and running input file for PACKMOL
It needs the Radius of the NP, also the size of the box, and number of
water molecules,
The output file has such a structure:

    tolerance TOLERANCE
    structure INPUT.FILE
       number NUMBER__OF_WATER_MOLECULES
       inside box X_MIN Y_MIN Z_MIN X_MAX Y_MAX Z_MAX
       outside sphere 0. 0. 0. NP_RADIUS+1
    end structure

output OUTPUT.FILE
INPUT.file:

    HEADER    water
    COMPND
    SOURCE
    HETATM    1  H   HOH     1       9.626   6.787  12.673
    HETATM    2  H   HOH     1       9.626   8.420  12.673
    HETATM    3  O   HOH     1      10.203   7.604  12.673
    CONECT    1    3
    CONECT    2    3
    CONECT    3    1    2
    END

The script kept the input file format as the standard PACKMOL PDB
input file for water.
The script first makes a PDB file, then converts it to the standard PDB
and data file structure.


The maximum radius is from silanized data.
TOLERANCE, X, Y, and Z limitations are from the "statinfo" module.
The number of molecules should be calculated from the final volume of
the water box based on the limitations.
"""

import sys
import subprocess
import pandas as pd
import static_info as stinfo
from colors_text import TextColor as bcolors


class InFile:
    """preparing input file for the"""
    def __init__(self,
                 radius: float  # Radius of the NP after silanization
                 ) -> None:
        self.__prepare_input(radius)

    def __prepare_input(self,
                        radius: float  # Radius of the NP after silanization
                        ) -> None:
        """prepare the input file for the PACKMOL"""
        self.__get_mols_num(radius)

    def __get_mols_num(self,
                       radius: float  # Radius of the NP after silanization
                      ) -> int:
        """get numbers of molecules based on the valume of the water
        box"""
        box_volume: float = self._get_box_valume()

    def _get_box_valume(self) -> float:
        """calculate the valume of the box including sphere's area"""
        x_lim: float = stinfo.Hydration.X_MAX - stinfo.Hydration.X_MIN
        y_lim: float = stinfo.Hydration.Y_MAX - stinfo.Hydration.Y_MIN
        z_lim: float = stinfo.Hydration.Z_MAX - stinfo.Hydration.Z_MIN
        box_valume: float = x_lim*y_lim*z_lim
        if box_valume <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                     f'\tZero valume, there in problem in setting box '
                     f'limitaion, box_volume is "{box_valume:.3f}"'
                     f'{bcolors.ENDC}')
