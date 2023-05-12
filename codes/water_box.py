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
import my_tools
import numpy as np
import run_packmol as pakml
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
        out_file: str = stinfo.Hydration.OUT_FILE
        inp_file: str = stinfo.Hydration.INP_FILE
        if stinfo.Hydration.CHECK_WATER_PDB:
            my_tools.check_file(out_file, delete=True)
            my_tools.check_file(inp_file, delete=True)
        with open(inp_file, 'w', encoding="utf8") as f_out:
            f_out.write('# Input file for PACKMOL, Water box for a NP ')
            f_out.write(f'with the radius of {self.radius}\n')
            f_out.write(f'# Contact angle: {stinfo.Hydration.CONATCT_ANGLE}\n')
            f_out.write('\nfiletype pdb\n')
            f_out.write(f'tolerance {stinfo.Hydration.TOLERANCE}\n')
            f_out.write(f'output {out_file}\n\n')
            self.__water_section(f_out, dimensions)
            self.__oil_section(f_out, dimensions)

    def __oil_section(self,
                      f_out: typing.IO,  # The file to write into it
                      dimensions: boxd.BoxEdges  # Num_moles, dims of box
                      ) -> None:
        """set the data for mols in water section"""
        if stinfo.Hydration.CONATCT_ANGLE > 0:
            self.__print_header(f_out, ' Oil ')
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
        self.__print_header(f_out, ' Water ')
        self.__write_inp_sections(f_out,
                                  dimensions.water_axis,
                                  stinfo.Hydration.WATER_PDB,
                                  dimensions.num_mols['sol'])
        self.__print_header(f_out, ' Counter ions ')
        self.__check_ions(f_out, dimensions)
        self.__check_oda(f_out, dimensions, 'oda')
        self.__check_nacl(f_out, dimensions)

    def __check_nacl(self,
                     f_out: typing.IO,  # The file to write into it
                     dimensions: boxd.BoxEdges  # Num_moles, dims of box
                     ) -> None:
        """check nacl charges and write its part"""
        if dimensions.num_mols['sal'] > 0:
            self.__print_header(f_out, ' NaCl ')
            for pdb in [stinfo.Hydration.CL_PDB, stinfo.Hydration.NA_PDB]:
                self.__write_inp_sections(f_out,
                                          dimensions.water_axis,
                                          pdb,
                                          dimensions.num_mols['sal'])

    def __check_ions(self,
                     f_out: typing.IO,  # The file to write into it
                     dimensions: boxd.BoxEdges  # Num_moles, dims of box
                     ) -> None:
        """check ions charges and write its part"""
        pdb_file: str  # Name of the pdbfile based on the charge of the system
        num_ions: int = dimensions.num_mols['ion']
        if num_ions > 0:
            pdb_file = stinfo.Hydration.CL_PDB
        else:
            pdb_file = stinfo.Hydration.NA_PDB
        self.__write_inp_sections(f_out,
                                  dimensions.water_axis,
                                  pdb_file,
                                  np.abs(num_ions))

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
                    if stinfo.Hydration.ODAP_PROTONATION:
                        self.__print_header(f_out,
                            ' Protonated Octadecylamine ')
                        pdb_file = stinfo.Hydration.ODAP_PDB
                    else:
                        self.__print_header(f_out,
                            ' Unprotonated Octadecylamine in Water ')
                        pdb_file = stinfo.Hydration.ODAN_PDB
                    dimens = dimensions.water_axis
                    self.__odap_interface(dimensions)
                elif style == 'odn':
                    self.__print_header(f_out, ' Unprotonated Octadecylamine ')
                    pdb_file = stinfo.Hydration.ODAN_PDB
                    dimens = dimensions.oil_axis
                self.__write_inp_sections(f_out,
                                          dimens,
                                          pdb_file,
                                          dimensions.num_mols[style])

    def __odap_interface(self,
                         dimensions: boxd.BoxEdges  # Num_moles, dims of box
                         ) -> dict[str, float]:
        """check if the ODAP should be at the interface"""
        dimes: dict[str, float]  # The dimension of the box for ODAP
        # if there is an interface:
        if stinfo.Hydration.CONATCT_ANGLE < 0:
            dimes = dimensions.water_axis
        else:
            # If they should be at the interface
            if not stinfo.Hydration.ODAP_INTERFACE:
                dimes = dimensions.water_axis
            else:
                dimes = self.__get_odap_area(dimensions)
        return dimes

    def __get_odap_area(self,
                        dimensions: boxd.BoxEdges  # Num_moles, dims of box
                        ) -> dict[str, float]:
        """calculate the area for the ODAP in the interface.
        Empirically the ODAP tail is entirely in the oil phase, and its
        head (NH3) is in the water phase. However, here, they will give
        an area bigger than their length in the interface since it is
        faster for PACKMOL to run.
        The area for the ODAP is split equally between the oil and water
        phase, but if the length of oil, water, or system cannot take
        that, it should try to solve it with a warning or exit with an
        error.
        """
        oda_length: float = self.__check_oda_box(dimensions)
        oda_box: dict[str, float]  # Edges of the box for the interface
        oda_box = dimensions.water_axis
        oda_box['z_lo'] = dimensions.water_axis['z_hi'] - oda_length / 2
        oda_box['z_hi'] = dimensions.oil_axis['z_lo'] + oda_length / 2
        return oda_box

    def __check_oda_box(self,
                        dimensions: boxd.BoxEdges  # Num_moles, dims of box
                        ) -> float:
        """check if the length of the selected area for the ODAP can
        be put in the system box"""
        oda_length: float = stinfo.Constants.ODA_length + 1
        box_z: float = dimensions.box_edges['box']['z_lim']/2
        oil_z: float = dimensions.box_edges['oil']['z_lim']
        sol_z: float = dimensions.box_edges['sol']['z_lim']
        for item, box in zip([box_z, oil_z, sol_z], ['box', 'oil', 'water']):
            if oda_length/2 > item:
                sys.exit(f'{bcolors.WARNING}{self.__class__.__name__}: '
                         f'({self.__module__})\n\tThe are of the interface '
                         f' cannot suit into the system {box}: {oda_length} >'
                         f' {box_z} {bcolors.ENDC}')
        return oda_length

    def __write_inp_sections(self,
                             f_out: typing.IO,  # The file to write into it
                             dimens: dict[str, float],  # Section dimensions
                             pdb_file: str,  # Name of the pdb file to write
                             num_mol: int  # Number of the molecules in secti
                             ) -> None:
        """write water section: which include the water, ions, and
        protonated ODA"""
        tlr: float = stinfo.Hydration.TOLERANCE
        expend_edg: float = 1 * tlr  # To expand edges, fasting the test run
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
                f_out.write(f'{dimens["x_lo"] - expend_edg : .2f} ')
                f_out.write(f'{dimens["y_lo"] - expend_edg : .2f} ')
                f_out.write(f'{dimens["z_lo"] : .2f} ')
                f_out.write(f'{dimens["x_hi"] + expend_edg : .2f} ')
                f_out.write(f'{dimens["y_hi"] + expend_edg : .2f} ')
                f_out.write(f'{dimens["z_hi"] : .2f}\n')
                f_out.write(
                    f'\toutside sphere 0. 0. 0. {self.radius + tlr: .2f}\n')
                f_out.write('end structure\n\n')

    def __print_header(self,
                       f_out: typing.IO,  # The file to write into it
                       section: str  # Name of the section
                       ) -> None:
        """write the name of section in the inp file"""
        row_length: int = 66
        f_out.write(f'{section.center(row_length, "#")}\n')

    def print_info(self) -> None:
        """print infos"""
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__}):\n'
              '\tThe input file for PACKMOL is written in: '
              f'"{stinfo.Hydration.INP_FILE}"\n'
              f'{bcolors.ENDC}')


if __name__ == "__main__":
    dims = boxd.BoxEdges(radius=10, net_charge=10)
    in_file = InFile(radius=10, dimensions=dims)
    water_b = pakml.RunPackMol(inp_file=stinfo.Hydration.INP_FILE,
                               out_file=stinfo.Hydration.OUT_FILE)
