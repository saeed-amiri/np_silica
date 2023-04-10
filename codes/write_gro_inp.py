"""write input files for simulation inputs for GROMACS
For topology file it needs:
    - itp of the forcefield (stinfo)
    - itp of ions, water, oil, nanoparticle if needed
    - the number of the molecules in [ molecules ] section
"""

import typing
import static_info as stinfo
import box_dimensions as boxd
from colors_text import TextColor as bcolors


class WriteTop:
    """write topology file"""
    def __init__(self,
                 mol_nums: dict[str, int]  # Number of each molecule in system
                 ) -> None:
        self.write_topo(mol_nums)

    def write_topo(self,
                   mol_nums: dict[str, int]  # Number of each molecule
                   ) -> None:
        """writ the data"""
        with open(stinfo.GroInp.TOPFILE,  'w', encoding='utf8') as f_out:
            f_out.write('; Topology file, radius XX and contact angle\n')
            self.__write_include(stinfo.GroInp.FORCEFIELD, 'forcefiled', f_out)
            if mol_nums['sol'] > 0:
                self.__write_include(stinfo.GroInp.WATERITP, 'water', f_out)

    def __write_include(self,
                        itp_file: str,  # Path of the forcefield
                        molecule: str,  # Name of the molecule
                        f_out: typing.IO  # Topology file
                        ) -> None:
        """write the include section of the forcefiled"""
        f_out.write(f'\n; Include {molecule} topology\n')
        f_out.write(f'#include "{itp_file}"\n')


if __name__ == '__main__':
    axis_limits = boxd.BoxEdges(radius=25, net_charge=10)
    topp = WriteTop(axis_limits.num_mols)
