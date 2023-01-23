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
                Atoms: pd.DataFrame,  # Atoms df in the form of lammps fullatom
                Sigroup: list[typing.Any],  # Name | index to select groups[Si]
                Ogroup:  list[typing.Any],  # Name | index groups[O] to delete
                fraction: float = 1  # Fraction of to select from, 0<fr<=1
                ) -> None:
        self.__get_silica(Atoms, Sigroup)

    def __get_silica(self,
                    Atoms: pd.DataFrame,  # Atoms df in the form of lammps fullatom
                    Sigroup: list[typing.Any],  # Name | index to select groups[Si]
                    ) -> None:
        """Get index or name of atoms"""
        Si_list: list[pd.DataFrame] = []  # df with asked Si groups
        for item in Sigroup:
            Si_list.append(Atoms[Atoms['name'] == item])
        df: pd.DataFrame = pd.concat(Si_list)  # Df of all the Si to replace
        max_radius: float = self.__get_radius(Atoms)
        self.__set_radii(Atoms, df)
    
    def __set_radii(self,
                    Atoms: pd.DataFrame,  # Atoms in lammps full atom
                    df: pd.DataFrame  # The Atoms_df for asked silicons
                    ) -> pd.DataFrame:
        """add radii for the atoms to drop inner ones later"""
        # Check the center of mass
        x_cm: float = np.average(Atoms['x'])
        y_cm: float = np.average(Atoms['y'])
        z_cm: float = np.average(Atoms['z'])
        EPSILON: float = 1e-5  # to consider as zero
        if x_cm > EPSILON or y_cm > EPSILON or z_cm > EPSILON:
            exit(f'{bcolors.FAIL}{self.__class__.__name__}\n'
                 f'\tThe Center of mass is not in zero, Sth is wrong '
                 f'in reading the main data file\n{bcolors.ENDC}')
        for index, row in df.iterrows():
            vec: np.array = ([row['x'], row['y'], row['z']])
            radii: float = np.linalg.norm(vec)
            print(radii)
    
    def __get_radius(self,
                    Atoms: pd.DataFrame,  # Atoms in lammps full atom
                    ) -> float:
        """return the radius of the nano-particles"""


if __name__ == '__main__':
    fname = sys.argv[1]
    data = GetData(fname)
    groups = GetGroups(data.Atoms_df,
                       Sigroup=['SD'],
                       Ogroup=['OD'],
                       fraction=1)