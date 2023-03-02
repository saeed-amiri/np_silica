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
import numpy as np
import static_info as stinfo
from colors_text import TextColor as bcolors


class NumberMols:
    """Getting the number of water molecules in the volume"""
    def __init__(self,
                 radius: float  # Radius of the NP after silanization
                 ) -> None:
        self.prepare_input(radius)
        self.print_info()

    def prepare_input(self,
                      radius: float  # Radius of the NP after silanization
                      ) -> None:
        """prepare the input file for the PACKMOL"""
        self.edge_cube: float  # Edge of the box that inscibed the NP
        self.number_mols: int  # Number of the molecules in the box
        self.number_mols, self.edge_cube = self.__get_mols_num(radius)

    def __get_mols_num(self,
                       radius: float  # Radius of the NP after silanization
                       ) -> tuple[int, float]:
        """get numbers of molecules based on the volume of the water
        box"""
        sphere_volume: float = self.__get_sphere_volume(radius)
        box_volume: float  # volume of the box
        edge_cube: float  # Edge of the box that inscibed the NP
        box_volume, edge_cube = self.__get_box_volume(sphere_volume)
        net_volume: float = self.__check_volumes(box_volume, sphere_volume)
        return self.__calc_mols_num(net_volume), edge_cube

    def __calc_mols_num(self,
                        volume: float  # Net volume of the water box
                        ) -> int:
        """return number of water molecules to get the density of
        water"""
        lit_m3: float = 1e-24  # convert units
        m_water: float  # Mass of the water in the volume
        m_water = volume * stinfo.Hydration.WATER_DENSITY * lit_m3
        num_moles: float
        num_moles = int(m_water * stinfo.Hydration.AVOGADRO /
                        stinfo.Hydration.WATER_MOLAR_MASS) + 1
        return int(num_moles*2.5)

    def __get_box_volume(self,
                         sphere_volume: float,  # Volume of the sphere
                         ) -> tuple[float, float]:
        """calculate the volume of the box including sphere's area
        For the largest possible sphere is inscribed in cube, the ratio
        of volumes is: V_sphere/V_cube = pi/6"""
        v_inscribed_box: float = 6*sphere_volume/np.pi
        edge_cube: float  # Edge of the cube that inscribed the sphere
        edge_cube = v_inscribed_box**(1/3)
        x_lim: float = (stinfo.Hydration.X_MAX -
                        stinfo.Hydration.X_MIN) + edge_cube
        y_lim: float = (stinfo.Hydration.Y_MAX -
                        stinfo.Hydration.Y_MIN) + edge_cube
        z_lim: float = (stinfo.Hydration.Z_MAX -
                        stinfo.Hydration.Z_MIN) + edge_cube
        box_volume: float = x_lim*y_lim*z_lim
        if box_volume <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                     f'\tZero volume, there in problem in setting box '
                     f'limitaion, box_volume is "{box_volume:.3f}"'
                     f'{bcolors.ENDC}')
        return box_volume, edge_cube

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
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}: ('
              f'{self.__module__}):\n'
              f'\tThe number of water molecules is set to '
              f'"{self.number_mols}"'
              f'{bcolors.ENDC}')


class InFile:
    """preparing input file for the"""
    def __init__(self,
                 radius: float,  # Radius of the NP
                 num_mols: int,  # Number of the molecules in the water volume
                 edge: float  # Edge of the inscribed cube
                 ) -> None:
        self.write_file(radius, num_mols, edge)
        self.print_info()

    def write_file(self,
                   radius: float,  # Radius of the NP
                   num_mols: int,  # Number of the molecules in the water volum
                   edge: float  # Edge of the inscribed cube
                   ) -> None:
        """write the input file for the PACKMOL"""
        tolerence: float = stinfo.Hydration.TOLERANCE
        out_file: str = 'water_box.pdb'
        with open(stinfo.Hydration.INP_FILE, 'w', encoding="utf8") as f_out:
            f_out.write('# Input file for PACKMOL, Water box for a NP ')
            f_out.write(f'with the radius of {radius}\n\n')
            f_out.write(f'tolerance {stinfo.Hydration.TOLERANCE}\n\n')
            f_out.write(f'structure {stinfo.Hydration.WATER_PDB}\n')
            f_out.write(f'\tnumber {num_mols}\n')
            f_out.write('\tinside box ')
            f_out.write(f'{stinfo.Hydration.X_MIN - edge - tolerence: .2f} ')
            f_out.write(f'{stinfo.Hydration.Y_MIN - edge - tolerence: .2f} ')
            f_out.write(f'{stinfo.Hydration.Z_MIN - edge - tolerence: .2f} ')
            f_out.write(f'{stinfo.Hydration.X_MAX + edge + tolerence: .2f} ')
            f_out.write(f'{stinfo.Hydration.Y_MAX + edge + tolerence: .2f} ')
            f_out.write(f'{stinfo.Hydration.Z_MAX + edge + tolerence: .2f}\n')
            f_out.write(f'\t outside sphere 0. 0. 0. {radius: .2f}\n')
            f_out.write('end structure\n\n')
            f_out.write(f'output {out_file}\n\n')

    def print_info(self) -> None:
        """print infos"""
        print(f'\n{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__}):\n'
              '\tThe input file for PACKMOL is written in: '
              f'"{stinfo.Hydration.INP_FILE}"\n'
              f'{bcolors.ENDC}')


class RunPackMol:
    """call PACKMOL and run the input script prepared for it"""
    def __init__(self) -> None:
        """run the input with subprocess"""
        pack_mol: str = stinfo.Hydration.PACKMOL  # Compiler of PACKMOL
        inp_file: str = stinfo.Hydration.INP_FILE  # Input file for packmol
        pack_flag: int  # If PACKMOL executed successfully
        pack_flag = self.make_water(pack_mol, inp_file)
        self.print_info(pack_flag)

    def make_water(self,
                   pack_mol: str,  # Compiler of PACKMOL
                   inp_file: str  # Input file for packmol
                   ) -> int:
        """call the subprocess and run the input file"""
        pack_flag: int = subprocess.call(f'{pack_mol} < {inp_file}>/dev/null',
                                         shell=True, cwd='./')
        return pack_flag

    def print_info(self,
                   pack_flag: int  # If PACKMOL executed successfully
                   ) -> None:
        """print infos"""
        if pack_flag == 0:
            print(f'\n{bcolors.OKCYAN}{self.__class__.__name__}:'
                  f' ({self.__module__})\n'
                  '\tPACKMOL executed successfully, output is: '
                  f'{stinfo.Hydration.OUT_FILE}'
                  f'{bcolors.ENDC}')
        else:
            sys.exit(f'{bcolors.FAIL}\tError! in executing PACKMOL'
                     f'{bcolors.ENDC}')


if __name__ == "__main__":
    moles = NumberMols(radius=50)
    in_file = InFile(radius=50,
                     num_mols=moles.number_mols,
                     edge=moles.edge_cube)
    water_box = RunPackMol()
