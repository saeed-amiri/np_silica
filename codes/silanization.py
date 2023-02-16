import re
import sys
import numpy as np
import pandas as pd
import write_lmp as wrlmp
import bond_check as bchek
import make_chains as mkchin
import update_coords as upcord
from colors_text import TextColor as bcolors


class ConcatAll:
    """append all the aminopropyle to the silicon data file"""
    def __init__(self,
                 silica: upcord.UpdateCoords,  # Silica updated
                 aminos: mkchin.PrepareAmino  # Update aminos
                 ) -> None:
        self.__concate_all(silica, aminos)
        self.__write_infos()

    def __concate_all(self,
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
        print(f'{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\tConcatenating all the atoms{bcolors.ENDC}')
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
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

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
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

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
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

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
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

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

        df = pd.DataFrame(columns=columns)
        df = pd.concat([df_silica, df_amino])
        return df

    def __set_attrs(self) -> None:
        """set attributes to object(self)"""
        self.Nmols: int = np.max(self.Atoms_df['mol'])
        self.NAtoms: int = len(self.Atoms_df)
        self.NBonds: int = len(self.Bonds_df)
        self.NAngles: int = len(self.Angles_df)
        self.NDihedrals: int = len(self.Dihedrals_df)
        self.NAtomTyp: int = np.max(self.Masses_df['typ'])
        self.NBondTyp: int = np.max(self.Bonds_df['typ'])
        self.NAngleTyp: int = np.max(self.Angles_df['typ'])
        self.NDihedralTyp: int = np.max(self.Dihedrals_df['typ'])

    def __write_infos(self) -> None:
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
    fname = sys.argv[1]
    update = upcord.UpdateCoords(fname)  # Updated data for silica
    amino = mkchin.GetAmino()
    up_aminos = mkchin.PrepareAmino(update, amino)
    silanized_data = ConcatAll(update, up_aminos)
    # if need to check the bonds:
    # bc = bchek.CheckBond(silanized_data)
    np_size: int = int(re.findall(r'\d+', fname)[0])
    fout: str = f'silanized_{np_size}nm.data'
    wrt = wrlmp.WriteLmp(obj=silanized_data, output=fout)
    wrt.write_lmp()
