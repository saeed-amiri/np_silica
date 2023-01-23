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


class GetGroups:
    """get silinol groups to add cahin to them"""
    def __init__(self,
                 Atoms: pd.DataFrame,  # Atoms df in form of lammps fullatom
                 Sigroup: list[typing.Any],  # Name | index to select group[Si]
                 Ogroup:  list[typing.Any],  # Name | index groups[O] to delete
                 fraction: float = 1  # Fraction of to select from, 0<fr<=1
                 ) -> None:
        self.__get_silica(Atoms, Sigroup)

    def __get_silica(self,
                     Atoms: pd.DataFrame,  # Atoms df in the lammps fullatom
                     Sigroup: list[typing.Any],  # Name | index of groups[Si]
                     ) -> None:
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
        print(df)

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


if __name__ == '__main__':
    fname = sys.argv[1]
    data = GetData(fname)
    groups = GetGroups(data.Atoms_df,
                       Sigroup=['SD'],
                       Ogroup=['OD'],
                       fraction=1)
