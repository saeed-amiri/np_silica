import sys
import numpy as np
import pandas as pd
from colors_text import TextColor as bcolors


class Doc:
    """update charges in the data file in which all the selected atoms
    are removed and prepared for the silanizations"""


class UpdateCharge:
    """set charges for all the types in the nano particles!"""
    def __init__(self,
                 atoms_df: pd.DataFrame,  # Updated atoms from update_coords
                 Si_df: pd.DataFrame  # Si df, selected for adding chains
                 ) -> None:
        self.__update_charges(atoms_df, Si_df)

    def __update_charges(self,
                         atoms_df: pd.DataFrame,  # Updated silica atoms coords
                         Si_df: pd.DataFrame  # Si df, for adding amino
                         ) -> pd.DataFrame:
        Si_df = self.__clean_si_df(Si_df)

    def __clean_si_df(self,
                      Si_df: pd.DataFrame
                      ) -> pd.DataFrame:
        """drop extra columns"""
        df = Si_df.copy()
        columns: list[str]  # Columns to drop
        columns = ['rho', 'azimuth', 'polar', 'ddd']
        for col in columns:
            try:
                df.drop(columns=col, axis=1, inplace=True)
            except KeyError:
                pass
        return df
