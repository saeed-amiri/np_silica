"""combine water_box.pdb with silanized nanoparticles
For this, the script calls PACKMOL and runs it.
Here the script only needs the name of the input files path of PACKMOL,
which all should be defined in the static_info and must be made by
previous tools such as silanization.py and hydration.py
"""


import sys
import typing
import run_packmol as pakml
import static_info as stinfo
from colors_text import TextColor as bcolors


class InFile:
    """preparing input file for PACKMOL"""
    def __init__(self,
                 silaniz_pdb: str  # Name of the file to put into the box
                 ) -> None:
        self.write_file(silaniz_pdb)
        self.print_info()

    def write_file(self,
                   silaniz_pdb: str  # Name of the file to put into the box
                   ) -> None:
        """write the input file for the PACKMOL, Subtract the number of
        water atoms so can fit the number of ions into box"""
        out_file: str = 'silica_water.pdb'
        with open(stinfo.Hydration.WS_INP, 'w', encoding="utf8") as f_out:
            f_out.write('# Input file for PACKMOL, Silanized and '
                        'Water box for a NP\n\n')
            f_out.write('filetype pdb\n')
            f_out.write(f'tolerance {stinfo.Hydration.TOLERANCE}\n')
            f_out.write(f'output {out_file}\n\n')
            self.__write_water(f_out)
            self.__write_nano_p(f_out, silaniz_pdb)

    def __write_nano_p(self,
                       f_out: typing.IO,  # The file to write into it
                       silaniz_pdb: str  # Name of the file to put into the box
                       ) -> None:
        """write the ions section in the box"""
        x_o: float = 0.0  # Position of the center for the nanoparticle
        y_o: float = 0.0  # Position of the center for the nanoparticle
        z_o: float = 0.0  # Position of the center for the nanoparticle
        a_o: float = 0.0  # Angle for rotation in radians
        b_o: float = 0.0  # Angle for rotation in radians
        g_o: float = 0.0  # Angle for rotation in radians
        f_out.write(f'structure {silaniz_pdb}\n')
        f_out.write('\tnumber 1\n')
        f_out.write('\tcenter\n')
        f_out.write('\tfixed ')
        f_out.write(f'{x_o} {y_o} {z_o} {a_o} {b_o} {g_o}\n')
        f_out.write('end structure\n\n')

    def __write_water(self,
                      f_out: typing.IO  # The file to write into it
                      ) -> None:
        """write the water box section in the packmol inputfile"""
        x_o: float = 0.0  # Position of the center for the nanoparticle
        y_o: float = 0.0  # Position of the center for the nanoparticle
        z_o: float = 0.0  # Position of the center for the nanoparticle
        a_o: float = 0.0  # Angle for rotation in radians
        b_o: float = 0.0  # Angle for rotation in radians
        g_o: float = 0.0  # Angle for rotation in radians
        f_out.write(f'structure {stinfo.Hydration.OUT_FILE}\n')
        f_out.write('\tnumber 1\n')
        f_out.write('\tcenter\n')
        f_out.write('\tfixed ')
        f_out.write(f'{x_o} {y_o} {z_o} {a_o} {b_o} {g_o}\n')
        f_out.write('end structure\n\n')

    def print_info(self) -> None:
        """print infos"""
        print(f'\n{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__}):\n'
              '\tThe input file for PACKMOL is written in: '
              f'"{stinfo.Hydration.WS_INP}"\n'
              f'{bcolors.ENDC}')


if __name__ == "__main__":
    in_file = InFile(silaniz_pdb=sys.argv[1])
    pakml.RunPackMol(inp_file=stinfo.Hydration.WS_INP,
                     out_file=stinfo.Hydration.GRO_PDB)
