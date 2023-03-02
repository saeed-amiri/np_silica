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
import numpy as np
import static_info as stinfo
from colors_text import TextColor as bcolors


class InFile:
    """preparing input file for the"""
    def __init__(self,
                 radius: float  # Radius of the NP after silanization
                 ) -> None:
        self.prepare_input(radius)
        self.print_info()

    def prepare_input(self,
                      radius: float  # Radius of the NP after silanization
                      ) -> None:
        """prepare the input file for the PACKMOL"""
        self.__get_mols_num(radius)

    def __get_mols_num(self,
                       radius: float  # Radius of the NP after silanization
                      ) -> int:
        """get numbers of molecules based on the volume of the water
        box"""
        box_volume: float = self.__get_box_volume()
        sphere_volume: float = self.__get_sphere_volume(radius)
        net_volume: float = self.__check_volumes(box_volume, sphere_volume)
        print(net_volume)

    def __get_box_volume(self) -> float:
        """calculate the volume of the box including sphere's area"""
        x_lim: float = stinfo.Hydration.X_MAX - stinfo.Hydration.X_MIN
        y_lim: float = stinfo.Hydration.Y_MAX - stinfo.Hydration.Y_MIN
        z_lim: float = stinfo.Hydration.Z_MAX - stinfo.Hydration.Z_MIN
        box_volume: float = x_lim*y_lim*z_lim
        if box_volume <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                     f'\tZero volume, there in problem in setting box '
                     f'limitaion, box_volume is "{box_volume:.3f}"'
                     f'{bcolors.ENDC}')
        return box_volume

    def __get_sphere_volume(self,
                            radius: float  # Radius of the NP after silanizatio
                            ) -> float:
        """calculate the volume of the sphere for the NP"""
        sphere_volume: float = 4*np.pi*(radius**3)/3
        if sphere_volume <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                     f'\tZero volume, there in problem in setting sphere '
                     f'valume, it is "{sphere_volume:.3f}"'
                     f'{bcolors.ENDC}')
        return sphere_volume

    def __check_volumes(self,
                        box: float,  # Volume of the box
                        sphere: float  # Volume of the sphere
                        ) -> float:
        """check wether the box volume is bigger then the sphere volume"""
        if box - sphere <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}: '
                     f'({self.__module__})\n'
                     f'\tVolume of the Sphere is less or equal to the '
                     f'Volume of the Box'
                     f'{bcolors.ENDC}\n')
        return box - sphere

    def print_info(self) -> None:
        """print infos"""
        print("TEST")

if __name__ == "__main__":
    infile = InFile(10)
