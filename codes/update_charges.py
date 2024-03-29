"""update charges in the data file in which all the selected atoms
    are removed and prepared for the silanizations"""


import pandas as pd
import static_info as stinfo
from colors_text import TextColor as bcolors


class UpdateCharge:
    """set charges for all the types in the nano particles!"""
    def __init__(self,
                 atoms_df: pd.DataFrame,  # Updated atoms from update_coords
                 si_df: pd.DataFrame,  # Si df, selected for adding chains
                 old_new_dict: dict[int, int]  # Old: new atoms id
                 ) -> None:
        self.Atoms_df: pd.DataFrame  # Atoms with updated charges
        si_df = self.update_si_df(si_df)
        self.Atoms_df = self.update_charges(atoms_df, si_df, old_new_dict)

    def update_charges(self,
                       atoms_df: pd.DataFrame,  # Updated silica atoms coords
                       si_df: pd.DataFrame,  # Si df, for adding amino
                       old_new_dict: dict[int, int]  # Old: new atoms id
                       ) -> pd.DataFrame:
        si_df = self.__clean_si_df(si_df)
        si_df = self.__update_id_si(si_df, old_new_dict)
        si_df = self.__update_id_om(si_df, old_new_dict)
        atoms_df = self.__update_si_charge(si_df, atoms_df)
        atoms_df = self.__update_om_charge(si_df, atoms_df)
        return atoms_df

    def update_si_df(self,
                     si_df: pd.DataFrame  # Si_df, selected for adding chains
                     ) -> pd.DataFrame:
        """get all the O bonded to the selected Si to update thier
        charges"""
        df_c: pd.DataFrame = si_df.copy()
        df_c['OM_q_list'] = [None for _ in si_df.index]
        for item, row in si_df.iterrows():
            if len(row['OM_replace']) == 3:
                df_c.at[item, 'OM_q_list'] = row['OM_replace']
            elif len(row['OM_name']) < 3:
                df_c.at[item, 'OM_q_list'] = self.__get_si_bonded_o(row)
        return df_c

    def __get_si_bonded_o(self,
                          row: pd.DataFrame  # A row of si_df to get the Ox
                          ) -> list[int]:
        """get the undrop ox from ox_list and add it(them) with
        OM_replace to a new list and return for OM_q_list"""
        undrop_ox: list[int]  # List of undroped Ox
        undrop_ox = [item for item in row['Ox_list'] if item != row['Ox_drop']]
        om_q_list: list[int] = []
        om_q_list.extend(row['OM_replace'])
        om_q_list.extend(undrop_ox)
        return om_q_list

    def __update_om_charge(self,
                           si_df: pd.DataFrame,  # With updated atoms id
                           atoms_df: pd.DataFrame,  # Updated silica atoms
                           ) -> pd.DataFrame:
        """update all the OM atoms in the O or O&H groups in the main
        atoms dataframe"""
        count_om: int = 0  # Get the number of OM with charge-change
        print(f'\t{bcolors.HEADER}{self.__class__.__name__}:')
        if stinfo.UpdateCharge.OM is None:
            print('\t\tOM charges are remain unchanged!')
        else:
            print('\t\tCharges of OM bonded to the Amine groups are set '
                  f'to {stinfo.UpdateCharge.OM}')
            for _, row in si_df.iterrows():
                for ind in row['OM_q_list']:
                    if (atoms_df.iloc[ind-1]['charge'] !=
                       stinfo.UpdateCharge.OM):
                        atoms_df.at[ind, 'charge'] = stinfo.UpdateCharge.OM
                        count_om += 1
        print(f'\t\t["{count_om}" O atoms changed during set charges]'
              f'{bcolors.ENDC}')
        return atoms_df

    def __update_si_charge(self,
                           si_df: pd.DataFrame,  # With updated atoms id
                           atoms_df: pd.DataFrame,  # Updated silica atoms
                           ) -> pd.DataFrame:
        """update all the Si atoms charges which lost O or O&H groups
        in the main atoms dataframe"""
        count_si: int = 0  # Get the number of si atoms with charge-change
        print(f'\t{bcolors.HEADER}{self.__class__.__name__}:')
        if stinfo.UpdateCharge.SI is None:
            print('\t\tSi charges that are remain unchanged!')
        else:
            print('\t\tCharges of Si bonded to the Amine groups are set '
                  f'to {stinfo.UpdateCharge.SI}')
            for item, _ in si_df.iterrows():
                if atoms_df.iloc[item-1]['charge'] != stinfo.UpdateCharge.SI:
                    atoms_df.at[item, 'charge'] = stinfo.UpdateCharge.SI
                    count_si += 1
        print(f'\t\t["{count_si}" Si atoms changed during set charges]'
              f'{bcolors.ENDC}')
        return atoms_df

    def __update_id_si(self,
                       si_df: pd.DataFrame,  # Si df, for adding amino
                       old_new_dict: dict[int, int]  # Old: new atoms id
                       ) -> pd.DataFrame:
        """update the atom_id with the updated ones in the atoms_df"""
        df: pd.DataFrame = si_df.copy()
        df['old_atom_id'] = df['atom_id']  # Old atom_id
        for item, _ in si_df.iterrows():
            df.at[item, 'atom_id'] = old_new_dict[item]
        df.index = df['atom_id']
        return df

    def __update_id_om(self,
                       si_df: pd.DataFrame,  # Si df, with updated Si id
                       old_new_dict: dict[int, int]  # Old: new atoms id
                       ) -> pd.DataFrame:
        """update the atom_id of the OM atoms in the list with new ones
        after the omission of Ox and H atoms"""
        df: pd.DataFrame = si_df.copy()
        df['OM_list_old'] = df['OM_q_list']
        om_list: list[int]  # To append new atom_id
        for item, row in si_df.iterrows():
            om_list = []
            for ind in row['OM_q_list']:
                om_list.append(old_new_dict[ind])
            df.at[item, 'OM_q_list'] = om_list
            del om_list
        return df

    def __clean_si_df(self,
                      si_df: pd.DataFrame
                      ) -> pd.DataFrame:
        """drop extra columns"""
        df = si_df.copy()
        columns: list[str]  # Columns to drop
        columns = ['rho', 'azimuth', 'polar', 'x', 'y', 'z', 'Ox_list']
        for col in columns:
            try:
                df.drop(columns=col, axis=1, inplace=True)
            except KeyError:
                pass
        return df
