"""write input files for simulation inputs for GROMACS
For topology file it needs:
    - itp of the forcefield (stinfo)
    - itp of ions, water, oil, nanoparticle if needed
    - the number of the molecules in [ molecules ] section
"""

import os
import typing
import static_info as stinfo
import box_dimensions as boxd
from colors_text import TextColor as bcolors


class WriteTop:
    """write topology file"""
    def __init__(self,
                 mol_nums: dict[str, int],  # Number of each molecule in system
                 net_charge: int,  # Total charge of silica system with sign
                 ) -> None:
        self.write_topo(mol_nums, net_charge)

    def write_topo(self,
                   mol_nums: dict[str, int],  # Number of each molecule
                   net_charge: int,  # Total charge of silica system with sign
                   ) -> None:
        """writ the data"""
        with open(stinfo.GroInp.TOPFILE,  'w', encoding='utf8') as f_out:
            f_out.write('; Topology file, radius XX and contact angle\n')
            self.__write_include(stinfo.GroInp.FORCEFIELD, 'forcefiled', f_out)
            if mol_nums['sol'] > 0:
                self.__write_include(stinfo.GroInp.WATERITP, 'water', f_out)
            if mol_nums['ion'] > 0 or mol_nums['sal']:
                self.__write_include(stinfo.GroInp.IONITP, 'ions', f_out)
            if mol_nums['oil'] > 0:
                itp_file = os.path.basename(stinfo.Hydration.OIL_ITP)
                self.__write_include(f'./{itp_file}', 'oil', f_out)
            if mol_nums['oda'] > 0:
                itp_file = os.path.basename(stinfo.Hydration.ODAP_ITP)
                self.__write_include(f'./{itp_file}', 'unprotonat ODA', f_out)
            if mol_nums['odn'] > 0:
                itp_file = os.path.basename(stinfo.Hydration.ODAN_ITP)
                self.__write_include(f'./{itp_file}', 'protonate ODA', f_out)
            self.__write_molecules(mol_nums, f_out, net_charge)

    def __write_include(self,
                        itp_file: str,  # Path of the forcefield
                        molecule: str,  # Name of the molecule
                        f_out: typing.IO  # Topology file
                        ) -> None:
        """write the include section of the forcefiled"""
        f_out.write(f'\n; Include {molecule} topology\n')
        f_out.write(f'#include "{itp_file}"\n')

    def __write_molecules(self,
                          mol_nums: dict[str, int],  # Number of each molecule
                          f_out: typing.IO,  # Topology file
                           net_charge: int,  # Charge of silica with sign
                          ) -> None:
        """write the molecule section"""
        f_out.write('\n\n[ molecules ]\n')
        f_out.write('; Compound\t\t\t#mols\n')
        ordered_list: list[str]  # Order the list in order to appear in topol
        ordered_list = ('sol', 'ion', 'oil', 'oda', 'odn', 'sal')
        for key in ordered_list:
            num = mol_nums[key]
            if num > 0:
                if key == 'sol':
                    mol_name = stinfo.PdbMass.water_residue
                    self.__write_mol_line(mol_name, num, f_out)
                if key == 'ion':
                    if net_charge > 0:
                        mol_name = stinfo.PdbMass.cl_residue
                        if mol_nums['sal'] > 0:
                            num += mol_nums['sal']
                    else:
                        mol_name = stinfo.PdbMass.na_residue
                    self.__write_mol_line(mol_name, num, f_out)
                if key == 'oil':
                    mol_name = stinfo.PdbMass.oil_residue
                    self.__write_mol_line(mol_name, num, f_out)
                if key == 'oda':
                    mol_name = stinfo.PdbMass.odap_residue
                    self.__write_mol_line(mol_name, num, f_out)
                if key == 'odn':
                    mol_name = stinfo.PdbMass.odan_residue
                    self.__write_mol_line(mol_name, num, f_out)
                if key == 'sal':
                    mol_name = stinfo.PdbMass.na_residue
                    self.__write_mol_line(mol_name, num, f_out)
                    if mol_nums['ion'] == 0:
                        mol_name = stinfo.PdbMass.cl_residue
                        self.__write_mol_line(mol_name, num, f_out)

    def __write_mol_line(self,
                         mol_name: str,  # Name of the molecule in itp and pdb
                         mol_num: int,  # Number of the molecules
                         f_out: typing.Any  # Output file
                         ) -> None:
        """write the line for molecule in topology file"""
        f_out.write(f'{mol_name:<15}{mol_num:>10}\n')


if __name__ == '__main__':
    axis_limits = boxd.BoxEdges(radius=25, net_charge=10)
    topp = WriteTop(axis_limits.num_mols, net_charge=10)
