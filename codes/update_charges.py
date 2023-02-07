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
                 Si_df: pd.DataFrame,  # Si df, selected for adding chains
                 old_new_dict: dict[int, int]  # Old: new atoms id
                 ) -> None:
        self.__update_charges(atoms_df, Si_df, old_new_dict)

    def __update_charges(self,
                         atoms_df: pd.DataFrame,  # Updated silica atoms coords
                         Si_df: pd.DataFrame,  # Si df, for adding amino
                         old_new_dict: dict[int, int]  # Old: new atoms id
                         ) -> pd.DataFrame:
        Si_df = self.__clean_si_df(Si_df)
        Si_df = self.__update_id_si(Si_df, old_new_dict)
        Si_df = self.__update_id_om(Si_df, old_new_dict)

    def __update_id_si(self,
                       Si_df: pd.DataFrame,  # Si df, for adding amino
                       old_new_dict: dict[int, int]  # Old: new atoms id
                       ) -> pd.DataFrame:
        """update the atom_id with the updated ones in the atoms_df"""
        df: pd.DataFrame = Si_df.copy()
        df['old_atom_id']: list[int]  # Old atom_id
        df['old_atom_id'] = df['atom_id']
        for item, _ in Si_df.iterrows():
            df.at[item, 'atom_id'] = old_new_dict[item]
        return df

    def __update_id_om(self,
                       Si_df: pd.DataFrame,  # Si df, with updated Si id
                       old_new_dict: dict[int, int]  # Old: new atoms id
                       ) -> pd.DataFrame:
        """update the atom_id of the OM atoms in the list with new ones
        after the omission of Ox and H atoms"""
        df: pd.DataFrame = Si_df.copy()
        df['OM_list_old'] = df['OM_list0']
        OM_list: list[int]  # To append new atom_id
        for item, row in Si_df.iterrows():
            OM_list = []
            for ind in row['OM_list0']:
                OM_list.append(old_new_dict[ind])
            df.at[item, 'OM_list0'] = OM_list
            del OM_list
        return df

    def __clean_si_df(self,
                      Si_df: pd.DataFrame
                      ) -> pd.DataFrame:
        """drop extra columns"""
        df = Si_df.copy()
        columns: list[str]  # Columns to drop
        columns = ['rho', 'azimuth', 'polar', 'x', 'y', 'z', 'Ox_list']
        for col in columns:
            try:
                df.drop(columns=col, axis=1, inplace=True)
            except KeyError:
                pass
        return df
