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
import typing
import subprocess
import static_info as stinfo
import box_dimensions as boxd
from colors_text import TextColor as bcolors


class InFile:
    """preparing input file for the PACKMOL execution"""
    def __init__(self,
                 radius: float,  # Radius of the NP
                 dimensions: boxd.BoxEdges  # Number of moles & dims of the box
                 ) -> None:
        self.radius: float = radius  # Radius of the silanized NP
        self.num_water: int = dimensions.num_mols['sol']
        self.write_file(dimensions)
        self.print_info()

    def write_file(self,
                   dimensions: boxd.BoxEdges  # Number of moles & dimensions
                   ) -> None:
        """write the input file for the PACKMOL, Subtract the number of
        water atoms so can fit the number of ions into box"""
        # For now, the scripts take the +1 as total charge of ODAp in itp
        # to be true, later it should it be checked by the scripts
        out_file: str = 'water_box.pdb'
        with open(stinfo.Hydration.INP_FILE, 'w', encoding="utf8") as f_out:
            f_out.write('# Input file for PACKMOL, Water box for a NP ')
            f_out.write(f'with the radius of {self.radius}\n\n')
            f_out.write('filetype pdb\n')
            f_out.write(f'tolerance {stinfo.Hydration.TOLERANCE}\n')
            f_out.write(f'output {out_file}\n\n')
            self.__water_section(f_out, dimensions)
            self.__oil_section(f_out, dimensions)

    def __oil_section(self,
                      f_out: typing.IO,  # The file to write into it
                      dimensions: boxd.BoxEdges  # Num_moles, dims of box
                      ) -> None:
        """set the data for mols in water section"""
        self.__write_inp_sections(f_out,
                                  dimensions.oil_axis,
                                  stinfo.Hydration.OIL_PDB,
                                  dimensions.num_mols['oil'])
        self.__check_oda(f_out, dimensions, 'odn')

    def __water_section(self,
                        f_out: typing.IO,  # The file to write into it
                        dimensions: boxd.BoxEdges  # Num_moles, dims of box
                        ) -> None:
        """set the data for mols in water section"""
        self.__write_inp_sections(f_out,
                                  dimensions.water_axis,
                                  stinfo.Hydration.WATER_PDB,
                                  dimensions.num_mols['sol'])
        self.__check_ions(f_out, dimensions)
        self.__check_oda(f_out, dimensions, 'oda')

    def __check_ions(self,
                     f_out: typing.IO,  # The file to write into it
                     dimensions: boxd.BoxEdges  # Num_moles, dims of box
                     ) -> None:
        """check ions charges and write it part"""
        pdb_file: str  # Name of the pdbfile based on the charge of the system
        num_ions: int = dimensions.num_mols['ion']
        if num_ions > 0:
            pdb_file = stinfo.Hydration.CL_PDB
        else:
            pdb_file = stinfo.Hydration.NA_PDB
        self.__write_inp_sections(f_out,
                                  dimensions.water_axis,
                                  pdb_file,
                                  dimensions.num_mols['ion'])

    def __check_oda(self,
                    f_out: typing.IO,  # The file to write into it
                    dimensions: boxd.BoxEdges,  # Num_moles, dims of box
                    style: str  # Protonated or not, ODA, ODAN
                    ) -> None:
        """check odap and write them if needed"""
        num_oda: int = dimensions.num_mols[style]
        if num_oda == 0:
            pass
        else:
            if num_oda < 0:
                sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}: '
                         f'({self.__module__})\n'
                         f'\tWrong number is set for the {style} molecules!\n'
                         f'{bcolors.ENDC}')
            else:
                if style == 'oda':
                    pdb_file = stinfo.Hydration.ODAP_PDB
                    dimens = dimensions.water_axis
                elif style == 'odn':
                    pdb_file = stinfo.Hydration.ODAN_PDB
                    dimens = dimensions.oil_axis
                self.__write_inp_sections(f_out,
                                          dimens,
                                          pdb_file,
                                          dimensions.num_mols[style])

    def __write_inp_sections(self,
                             f_out: typing.IO,  # The file to write into it
                             dimens: dict[str, float],  # Section dimensions
                             pdb_file: str,  # Name of the pdb file to write
                             num_mol: int  # Number of the molecules in secti
                             ) -> None:
        """write water section: which include the water, ions, and
        protonated ODA"""
        tlr: float = stinfo.Hydration.TOLERANCE
        if num_mol == 0:
            pass
        else:
            if num_mol < 0:
                sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}: '
                         f'({self.__module__})\n'
                         f'\tWrong number is set for the {num_mol} molecules!'
                         f'{bcolors.ENDC}\n')
            else:
                f_out.write(f'structure {pdb_file}\n')
                f_out.write(f'\tnumber {num_mol}\n')
                f_out.write('\tinside box ')
                f_out.write(f'{dimens["x_lo"] - tlr : .2f} ')
                f_out.write(f'{dimens["y_lo"] - tlr : .2f} ')
                f_out.write(f'{dimens["z_lo"] : .2f} ')
                f_out.write(f'{dimens["x_hi"] + tlr : .2f} ')
                f_out.write(f'{dimens["y_hi"] + tlr : .2f} ')
                f_out.write(f'{dimens["z_hi"] : .2f}\n')
                f_out.write(
                    f'\toutside sphere 0. 0. 0. {self.radius + tlr: .2f}\n')
                f_out.write('end structure\n\n')

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
            print(f'{bcolors.OKCYAN}{self.__class__.__name__}:'
                  f' ({self.__module__})\n'
                  '\tPACKMOL executed successfully, output is: '
                  f'"{stinfo.Hydration.OUT_FILE}"'
                  f'{bcolors.ENDC}')
        else:
            sys.exit(f'{bcolors.FAIL}\tError! in executing PACKMOL'
                     f'{bcolors.ENDC}')


if __name__ == "__main__":
    dims = boxd.BoxEdges(radius=20, net_charge=10)
    in_file = InFile(radius=20, dimensions=dims)

    water_box = RunPackMol()
