"""run PACKMOL package to combine system and make oil, water, ions,
surfactants"""


import os
import sys
import subprocess
import static_info as stinfo
from colors_text import TextColor as bcolors


class RunPackMol:
    """call PACKMOL and run the input script prepared for it"""
    def __init__(self,
                 inp_file: str,  # Name the input file for PACKMOL
                 out_file: str  # Name of the out put of PACKMOL
                 ) -> None:
        """run the input with subprocess
        inputs:
            inp_file name from stinfo to write PACKMOL file in it
            out_file name from stinfo to check if the PACKMOL ran
        """
        pack_mol: str = stinfo.Hydration.PACKMOL  # Compiler of PACKMOL
        pack_flag: int  # If PACKMOL executed successfully
        pack_flag = self.make_water(pack_mol, inp_file, out_file)
        self.print_info(pack_flag)

    def make_water(self,
                   pack_mol: str,  # Compiler of PACKMOL
                   inp_file: str,  # Input file for packmol
                   out_file: str  # Name of the out put of PACKMOL
                   ) -> int:
        """call the subprocess and run the input file"""
        self.__check_file(out_file, delete=True)
        pack_flag: int  # Check if PACKMOL executed successfully
        subprocess.call(f'{pack_mol} < {inp_file}>/dev/null',
                        shell=True, cwd='./')
        pack_flag = self.__check_file(out_file, delete=False)
        return pack_flag

    def __check_file(self,
                     out_file: str,  # Name of the out put of PACKMOL
                     delete: bool = False  # Keep the file or not
                     ) -> int:
        """check if water box exist, if delete"""
        pack_flag: int = -1  # Check if PACKMOL executed successfully
        if delete:
            if os.path.isfile(out_file):
                print(f'{bcolors.CAUTION}{self.__class__.__name__} '
                      f'({self.__module__})\n'
                      f'\tAn old "{out_file}" exists, it will be deleted'
                      f'{bcolors.ENDC}')
                os.remove(out_file)
        else:
            if os.path.isfile(out_file):
                pack_flag = 0
        return pack_flag

    def print_info(self,
                   pack_flag: int  # If PACKMOL executed successfully
                   ) -> None:
        """print infos"""
        if pack_flag == 0:
            print(f'{bcolors.OKCYAN}{self.__class__.__name__}: '
                  f'({self.__module__})\n'
                  '\tPACKMOL executed successfully, output is: '
                  f'"{stinfo.Hydration.GRO_PDB}"'
                  f'{bcolors.ENDC}')
        else:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}: '
                     f'({self.__module__})\n'
                     f'\tError! in executing PACKMOL\n'
                     f'{bcolors.ENDC}')


if __name__ == '__main__':
    pckml = RunPackMol(inp_file='water_box.inp', out_file='water_box.pdb')
