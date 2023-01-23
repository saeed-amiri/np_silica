import re
import sys
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
        self.__get_atoms(Atoms, Sigroup, Ogroup, fraction)

    def __get_atoms(self,
                    Atoms: pd.DataFrame,  # Atoms df in the form of lammps fullatom
                    Sigroup: list[typing.Any],  # Name | index to select groups[Si]
                    Ogroup:  list[typing.Any],  # Name | index groups[O] to delete
                    fraction: float = 1  # Fraction of to select from, 0<fr<=1
                    ) -> None:
        """Get index or name of atoms"""
        for item in Sigroup:
            df: pd.DataFrame = Atoms[Atoms['name'] == item]
        
        


if __name__ == '__main__':
    fname = sys.argv[1]
    data = GetData(fname)
    groups = GetGroups(data.Atoms_df,
                       Sigroup=['SD'],
                       Ogroup=['OD'],
                       fraction=1)