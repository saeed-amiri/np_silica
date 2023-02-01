import numpy as np
import pandas as pd
import update_coords as upcord
from colors_text import TextColor as bcolors


class DropOM:
    """drop OM in aminos which have less than three OM atoms
    The extra atoms will drop, and bonds, angles, dihedrals
    must be updated"""
    def __init__(self,
                 amino,  # From Get amino, no annotaion for loop import!
                 si_row: pd.DataFrame,  # One row of si df
                 del_OM: int  # Number of OM to drop
                 ) -> None:
        self.__drop_OM(amino, si_row, del_OM)

    def __drop_OM(self,
                  amino,  # From Get amino, no annotaion for loop import!
                  si_row: pd.DataFrame,  # One row of si df
                  del_OM: int  # Number of OM to drop
                  ) -> None:
        """drop OM atoms and update all the rest"""
        amino_OM: int = 3  # Number of OM in amino data file
        print(si_row['atom_id'], del_OM)
        OM_df: pd.DataFrame = amino.Atoms_df[amino.Atoms_df['name'] == 'OM']
        print(OM_df)
