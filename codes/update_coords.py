import re
import sys
import math
import typing
import numpy as np
import pandas as pd
import read_lmp_data as rdlmp
from colors_text import TextColor as bcolors


class Doc:
    """read data files and update positions and index of them"""


class GetData(rdlmp.ReadData):
    """read data file in LAMMPS format"""
    def __init__(self,
                 fname: str  # Name of the data file
                 ) -> None:
        super().__init__(fname)


class GetOGroups:
    """get Oxygen and/or Hydrogen groups to delete them and update
    data file.
    Oxygen is bonded to the silica and Hydrogen if asked for is bonded
    to the Oxygen.
    """
    def __init__(self,
                 silica: rdlmp.ReadData,  # Atoms df in form of lammps fullatom
                 df_si: pd.DataFrame,  # df with selected group[Si]
                 Ogroup:  list[typing.Any],  # Name | index groups[O] to delete
                 Hgroup:  list[typing.Any],  # Name | index groups[H] to delete
                 fraction: float = 1  # Fraction of to select from, 0<fr<=1
                 ) -> None:
        self.df_O = self.__get_oxgygen(silica, df_si, Ogroup)

    def __get_oxgygen(self,
                      silica: rdlmp.ReadData,  # Atoms df in lammps full atom
                      df_si: pd.DataFrame,  # df with selected group[Si]
                      Ogroup:  list[typing.Any],  # Name|index groups[O]
                      ) -> None:
        """Find the hydrogen which have bonds with the selected Silicas"""
        O_list: list[pd.DataFrame] = []  # df of all Oxygen groups
        Atoms = silica.Atoms_df
        for item in Ogroup:
            O_list.append(Atoms[Atoms['name'] == item])
        df: pd.DataFrame = pd.concat(O_list)
        O_delete = self.__get_o_delete(silica.Bonds_df, df_si, df)
        self.__delete_atoms(silica, O_delete)

    def __get_o_delete(self,
                       bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                       df_si: pd.DataFrame,  # df with selected group[Si]
                       df: pd.DataFrame  # DF with selected Oxygen
                       ) -> list[int]:  # index of the O to delete
        # get bonds
        all_l: list[int]  # All atoms in O and Si list
        all_Si = [item for item in df_si['atom_id']]
        all_O = [item for item in df['atom_id']]
        all_l = all_Si
        all_l.extend(all_O)
        delte_list: list[typing.Any] = []  # Rows of bonds DF
        for _, row in bonds_df.iterrows():
            if row['ai'] in all_l and row['aj'] in all_l:
                if row['ai'] in all_O:
                    if row['ai'] not in delte_list:
                        delte_list.append(row['ai'])
                if row['aj'] in all_O:
                    if row['ai'] not in delte_list:
                        delte_list.append(row['ai'])
        return delte_list

    def __delete_atoms(self,
                       silica: rdlmp.ReadData,  # Atoms df in lammps full atom
                       delete_group: list[int]  # Index of the atom to delete
                       ) -> pd.DataFrame:
        """delete atoms and return update version in LAMMPS format"""
        self.__update_atoms(silica.Atoms_df, delete_group)

    def __update_atoms(self,
                       Atoms_df: pd.DataFrame,  # Atoms in LAMMPS format
                       delete_group: list[int]  # Index of atom to delete
                       ) -> pd.DataFrame:
        """delete atoms and return updated one with a new column"""
        Atoms_df['old_atom_id'] = Atoms_df['atom_id']
        df = Atoms_df.copy()
        for item, row in df.iterrows():
            if row['atom_id'] in delete_group:
                Atoms_df.drop(index=[item+1], axis=0, inplace=True)
        print(f'{bcolors.OKBLUE}{self.__class__.__name__}\n'
              f'\t {len(delete_group)} is deleted from data file\n'
              f'{bcolors.ENDC}')
        del df
        Atoms_df.reset_index(inplace=True)
        Atoms_df.index += 1
        Atoms_df['atom_id'] = Atoms_df.index
        Atoms_df.drop(columns=['index'], inplace=True)
        return Atoms_df


class GetSiGroups:
    """get silinol groups to add cahin to them"""
    def __init__(self,
                 Atoms: pd.DataFrame,  # Atoms df in form of lammps fullatom
                 Sigroup: list[typing.Any],  # Name | index to select group[Si]
                 fraction: float = 1  # Fraction of Silica to remove
                 ) -> None:
        self.df_Si = self.__get_silica(Atoms, Sigroup)

    def __get_silica(self,
                     Atoms: pd.DataFrame,  # Atoms df in the lammps fullatom
                     Sigroup: list[typing.Any],  # Name | index of groups[Si]
                     ) -> pd.DataFrame:
        """Get index or name of atoms"""
        Si_list: list[pd.DataFrame] = []  # df with asked Si groups
        for item in Sigroup:
            Si_list.append(Atoms[Atoms['name'] == item])
        df: pd.DataFrame = pd.concat(Si_list)  # Df of all the Si to replace
        # Drop some columns
        df.drop(axis=1, columns=['nx', 'ny', 'nz', 'cmt', 'b_name'],
                inplace=True)
        print(f'{bcolors.OKGREEN}\tThere are: {len(df)} atoms with selected '
              f'atoms name in the file\n{bcolors.ENDC}')
        max_radius: float = self.__get_radius(Atoms)
        df = self.__get_angles(df)
        df = df[(df['rho'] >= max_radius - 5)]
        print(f'{bcolors.OKGREEN}\tThere are: {len(df)} selected atoms '
              f'in the choosen area of the system, Max_radius = {max_radius}'
              f'\n{bcolors.ENDC}')
        # Get Azimuth and Polar angle of each atom in df
        return df

    def __get_radius(self,
                     Atoms: pd.DataFrame,  # Atoms in lammps full atom
                     ) -> float:
        """return the radius of the nano-particles"""
        x_max: float = Atoms['x'].abs().max()
        y_max: float = Atoms['y'].abs().max()
        z_max: float = Atoms['z'].abs().max()
        return np.max([x_max, y_max, z_max])

    def __get_angles(self,
                     df: pd.DataFrame  # The selected Si group
                     ) -> pd.DataFrame:
        """find and set the azimuth and polar angle of the atoms"""
        rho: list[float] = []  # To add to dataframe
        azimuth: list[float] = []  # To add to dataframe
        polar: list[float] = []  # To add to dataframe
        for _, row in df.iterrows():
            x = row['x']
            y = row['y']
            z = row['z']
            i_rho: float = np.sqrt(x*x + y*y + z*z)
            i_azimuth: float = np.arctan2(x, y)
            i_polar: float = np.arccos(z/i_rho)
            rho.append(i_rho)
            azimuth.append(i_azimuth)
            polar.append(i_polar)

        df['rho'] = rho
        df['azimuth'] = azimuth
        df['polar'] = polar
        return df


class GetAmino(rdlmp.ReadData):
    """read the main Aminopropyle coordinates and put the Si position
    to zero"""
    def __init__(self) -> None:
        fname: str = '/scratch/saeed/MyScripts/np_silica/data/aminopropyl.data'
        super().__init__(fname)
        self.Si = 'Si'
        self.__to_origin()

    def __to_origin(self) -> pd.DataFrame:
        """put the coordinate of Si in Aminopropyle to zero"""
        df: pd.DataFrame = self.Atoms_df.copy()
        print(f'{bcolors.OKCYAN}\tMove Aminopropyle [Si] to origin\n'
              f'{bcolors.ENDC}')
        x_si: float = self.Atoms_df[self.Atoms_df['name'] == self.Si]['x'][1]
        y_si: float = self.Atoms_df[self.Atoms_df['name'] == self.Si]['y'][1]
        z_si: float = self.Atoms_df[self.Atoms_df['name'] == self.Si]['z'][1]
        df: pd.DataFrame = self.Atoms_df.copy()
        df['x'] -= x_si
        df['y'] -= y_si
        df['z'] -= z_si
        return df


if __name__ == '__main__':
    fname = sys.argv[1]
    silica = GetData(fname)
    si_groups = GetSiGroups(silica.Atoms_df,
                            Sigroup=['SD'],
                            fraction=1)
    oxides = GetOGroups(silica,
                        si_groups.df_Si,
                        Ogroup=['OD'],
                        Hgroup=['HO'],
                        fraction=1)
    amino = GetAmino()
