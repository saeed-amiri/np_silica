"""silanization of the nanoparticles"""


import re
import sys
import numpy as np
import pandas as pd
import set_masses_name
import write_lmp as wrlmp
import bond_check as bchek
import static_info as stinfo
import make_chains as mkchin
import update_coords as upcord
from colors_text import TextColor as bcolors


class ConcatAll:
    """append all the aminopropyle to the silicon data file"""
    def __init__(self,
                 silica: upcord.UpdateCoords,  # Silica updated
                 aminos: mkchin.PrepareAmino  # Update aminos
                 ) -> None:
        self.Nmols: int  # Number in the data file
        self.NAtoms: int  # Number in the data file
        self.NBonds: int  # Number in the data file
        self.NAngles: int  # Number in the data file
        self.NDihedrals: int  # Number in the data file
        self.NAtomTyp: int  # Number in the data file
        self.NBondTyp: int  # Number in the data file
        self.NAngleTyp: int  # Number in the data file
        self.NDihedralTyp: int  # Number in the data file
        self.concate_all(silica, aminos)
        self.print_infos()

    def concate_all(self,
                    silica: upcord.UpdateCoords,  # Silica updated
                    aminos: mkchin.PrepareAmino  # Update aminos
                    ) -> None:
        """concate the all atoms, bonds, angles, diedrlas"""
        self.Atoms_df: pd.DataFrame  # DF in write_lmp format
        self.Atoms_df = self.__concate_atoms(silica.Atoms_df,
                                             aminos.All_amino_atoms)
        self.Bonds_df = self.__concate_bonds(silica.Bonds_df,
                                             aminos.All_amino_bonds)
        self.Angles_df = self.__concate_angles(silica.Angles_df,
                                               aminos.All_amino_angles)
        self.Dihedrals_df = self.__concate_dihedrals(silica.Dihedrals_df,
                                                     aminos.All_amino_dihedrals
                                                     )
        self.Masses_df = self.__concate_masses(silica.Masses_df,
                                               aminos.Masses_df)
        self.max_radius: float = self.__get_max_radius()  # Max radius of NP
        print(f'{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\tConcatenating all the atoms, max radius: '
              f'"{self.max_radius:.3f}"'
              f'{bcolors.ENDC}')
        self.__set_attrs()

    def __concate_atoms(self,
                        silica_atoms: pd.DataFrame,  # Silica atoms
                        aminos_atoms: pd.DataFrame  # Aminos atoms
                        ) -> pd.DataFrame:
        """concate all the atoms, make sure they all have same columns"""
        columns: list[str]  # Name of the wanted columns
        columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z',
                   'nx', 'ny', 'nz', 'cmt', 'name']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_atoms[col]
            amino_df[col] = aminos_atoms[col]
        df_i: pd.DataFrame  # All atoms dataframe
        df_i = pd.concat([si_df, amino_df], ignore_index=True)
        df_i.index += 1
        del si_df
        del amino_df
        return df_i

    def __concate_bonds(self,
                        silica_bonds: pd.DataFrame,  # Silica bonds
                        aminos_bonds: pd.DataFrame  # Aminos bonds
                        ) -> pd.DataFrame:
        """concate all the bonds in the write_lmp format"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'cmt', 'name']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_bonds[col]
            amino_df[col] = aminos_bonds[col]
        df_i: pd.DataFrame  # All atoms dataframe
        df_i = pd.concat([si_df, amino_df], ignore_index=True)
        df_i.index += 1
        del si_df
        del amino_df
        return df_i

    def __concate_angles(self,
                         silica_angles: pd.DataFrame,  # Silica bonds
                         aminos_angles: pd.DataFrame  # Aminos bonds
                         ) -> pd.DataFrame:
        """concate all the angles in the write_lmp format"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'ak']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_angles[col]
            amino_df[col] = aminos_angles[col]
        df_i: pd.DataFrame  # All atoms dataframe
        df_i = pd.concat([si_df, amino_df], ignore_index=True)
        df_i.index += 1
        del si_df
        del amino_df
        return df_i

    def __concate_dihedrals(self,
                            silica_dihedrals: pd.DataFrame,  # Silica bonds
                            aminos_dihedrals: pd.DataFrame  # Aminos bonds
                            ) -> pd.DataFrame:
        """concate all the angles in the write_lmp format"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'ak', 'ah']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_dihedrals[col]
            amino_df[col] = aminos_dihedrals[col]
        df_i: pd.DataFrame  # All atoms dataframe
        df_i = pd.concat([si_df, amino_df], ignore_index=True)
        df_i.index += 1
        del si_df
        del amino_df
        return df_i

    def __concate_masses(self,
                         silica_masses: pd.DataFrame,  # Silica masses
                         amino_masses: pd.DataFrame  # Aminos masses
                         ) -> pd.DataFrame:
        """update and make total masses df"""
        columns: list[str] = ['typ', 'mass', 'cmt', 'name']
        df_silica = pd.DataFrame(columns=columns)
        df_amino = pd.DataFrame(columns=columns)
        for col in columns:
            df_silica[col] = silica_masses[col]
            df_amino[col] = amino_masses[col]

        df_c = pd.DataFrame(columns=columns)
        df_c = pd.concat([df_silica, df_amino], ignore_index=True)
        df_c.index += 1
        df_c = set_masses_name.SetMasses(df_c).df_masses
        return df_c

    def __get_max_radius(self) -> float:
        """find the maximum radius of the NP
        Doing it in a simple way, find the maximum of the absolute
        values in each direction and return the maximum one."""
        axis: list[str] = ['x', 'y', 'z']
        max_list: list[float] = []  # All the max along each direction
        for a_x in axis:
            i_max: float = np.max([np.abs(np.max(self.Atoms_df[a_x])),
                                   np.abs(np.min(self.Atoms_df[a_x]))])
            max_list.append(i_max)
        return np.max(max_list)

    def __set_attrs(self) -> None:
        """set attributes to object(self)"""
        self.Nmols = np.max(self.Atoms_df['mol'])
        self.NAtoms = len(self.Atoms_df)
        self.NBonds = len(self.Bonds_df)
        self.NAngles = len(self.Angles_df)
        self.NDihedrals = len(self.Dihedrals_df)
        self.NAtomTyp = np.max(self.Masses_df['typ'])
        self.NBondTyp = np.max(self.Bonds_df['typ'])
        self.NAngleTyp = np.max(self.Angles_df['typ'])
        self.NDihedralTyp = np.max(self.Dihedrals_df['typ'])

    def print_infos(self) -> None:
        """print info on the stderr"""
        print(f'{bcolors.OKGREEN}\tData Summary after silanization:\n'
              f'\t\t# Atoms: {self.NAtoms}, # Atom`s types: {self.NAtomTyp}\n'
              f'\t\t# Bonds: {self.NBonds}, # Bond`s types: {self.NBondTyp}\n'
              f'\t\t# Angles: {self.NAngles}, '
              f'# Angle`s types: {self.NAngleTyp}\n'
              f'\t\t# Dihedrals: {self.NDihedrals}, '
              f'# Dihedral`s types: {self.NDihedralTyp}\n'
              f'\t\tTotal charge: {self.Atoms_df["charge"].sum():.4f}\n'
              f'\t\tMin charge: {self.Atoms_df["charge"].min()}\n'
              f'\t\tMax charge: {self.Atoms_df["charge"].max()}'
              f'{bcolors.ENDC}'
              )


if __name__ == '__main__':
    F_NAME = sys.argv[1]
    update = upcord.UpdateCoords(F_NAME)  # Updated data for silica
    amino = mkchin.GetAmino()
    up_aminos = mkchin.PrepareAmino(update, amino)
    silanized_data = ConcatAll(update, up_aminos)
    CHECK_BOND: bool = False  # To check the bonds
    if CHECK_BOND:
        bchek.CheckBond(silanized_data)
    np_size: int = int(re.findall(r'\d+', F_NAME)[0])
    fout: str = f'silanized_{np_size}nm_g{stinfo.Constants.Coverage}.data'
    wrt = wrlmp.WriteLmp(obj=silanized_data, output=fout)
    wrt.write_lmp()
