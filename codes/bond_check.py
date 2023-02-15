import typing
import numpy as np
import pandas as pd
import concurrent.futures
from colors_text import TextColor as bcolors


class Doc:
    """sanity check!
    Check the length of all the bonds.
    There were some strange long bonds in the data file sometimes,
    This could be because of the random selection of dropping Si from
    silanization or maybe because of the python dataframe selection
    algorithm"""


class CheckBond:
    """calculate all length of all the bonds"""
    def __init__(self,
                 data  # In the form of the LAMMPS data, Atoms_df, Bonds_df
                 ) -> None:
        self.__check_bonds(data)

    def __check_bonds(self,
                      data  # In the form of the LAMMPS data
                      ) -> None:
        """get the all the bonds"""
        print(data.Atoms_df)
