import numpy as np
import pandas as pd
from colors_text import TextColor as bcolors


class DropOM:
    """drop OM in aminos which have less than three OM atoms
    The extra atoms will drop, and bonds, angles, dihedrals
    must be updated"""
    def __init__(self,
                 amino,  # From Get amino, no annotaion for loop import!
                 si_row: pd.DataFrame,  # One row of si df
                 del_OM: int,  # Number of OM to drop
                 OM_xyz: pd.DataFrame  # XYZ info of OM atoms for amino
                 ) -> None:
        self.__drop_OM(amino, si_row, del_OM, OM_xyz)

    def __drop_OM(self,
                  amino,  # From Get amino, no annotaion for loop import!
                  si_row: pd.DataFrame,  # One row of si df
                  del_OM: int,  # Number of OM to drop
                  OM_xyz: pd.DataFrame  # XYZ info of OM atoms for amino
                  ) -> None:
        """drop OM atoms and update all the rest"""
        OM_index: list[int]  # Index of OM n amino df
        OM_index = [item for item in
                    amino.Atoms_df[amino.Atoms_df['name'] == 'OM'].index]
        amino.Atoms_df = self.__set_OM_info(amino.Atoms_df,
                                            si_row,
                                            OM_xyz,
                                            OM_index)
        df = self.__drop_om_atoms(del_OM,
                                  OM_index,
                                  si_row,
                                  amino.Atoms_df.copy())

        df = self.__update_atom_ind(df, si_row['OM_list'])
        old_new_dict: dict[int, int]  # Index of new and old index of atoms
        old_new_dict = {k: v for k, v in zip(df['undrop_ind'], df['atom_id'])}
        df.drop(axis=1, columns=['undrop_ind'], inplace=True)
        self.Atoms_df = df
        bonds_df: pd.DataFrame = amino.Bonds_df.copy()
        self.Bonds_df = self.__drop_bonds(old_new_dict,
                                          bonds_df,
                                          OM_index)
        del df
        angles_df: pd.DataFrame = amino.Angles_df.copy()
        self.Angles_df = self.__drop_angles(old_new_dict, OM_index, angles_df)
        dihedrals_df: pd.DataFrame = amino.Dihedrals_df.copy()
        self.Dihedrals_df = self.__drop_dihedrals(old_new_dict,
                                                  OM_index,
                                                  dihedrals_df)

    def __drop_om_atoms(self,
                        del_OM: int,  # Number of OM to drop
                        OM_index: list[int],  # Index of OM in amino
                        si_row: pd.DataFrame,  # One row of si df
                        atom_df: pd.DataFrame  # Amino atoms df
                        ) -> pd.DataFrame:
        """drop extra OM from amino dataframe"""
        if del_OM < 3:
            for i in OM_index:
                if atom_df.iloc[i-1]['atom_id'] \
                   not in si_row['OM_list']:
                    atom_df.drop(axis=0, index=i, inplace=True)
        else:
            pass
        atom_df.reset_index(inplace=True)
        atom_df.rename(columns={'index': 'undrop_ind'}, inplace=True)
        atom_df.index += 1
        return atom_df

    def __update_atom_ind(self,
                          df: pd.DataFrame,  # Amino atoms df with removed OM
                          OM_list: list[int]  # Index of the OM in NP
                          ) -> pd.DataFrame:
        """update the index of df after OM was droped"""
        for item, row in df.iterrows():
            if row['name'] != 'Si':
                if row['atom_id'] not in OM_list:
                    df.at[item, 'atom_id'] = item
                    df.at[item, 'old_id'] = item
        return df

    def __drop_bonds(self,
                     old_new_dict: dict[int, int],  # Of old and updated index
                     df: pd.DataFrame,  # Bonds df of amino
                     OM_index: list[int]  # Index of the OM atoms
                     ) -> pd.DataFrame:
        """drop and update bonds df"""
        df_ = df.copy()
        for item, row in df_.iterrows():
            if row['ai'] in OM_index or row['aj'] in OM_index:
                df.drop(axis=0, index=item, inplace=True)
        del df_
        df_ = df.copy()
        for item, row in df_.iterrows():
            if row['ai'] >= np.min(OM_index):
                df.at[item, 'ai'] = old_new_dict[row['ai']]
            if row['aj'] >= np.min(OM_index):
                df.at[item, 'aj'] = old_new_dict[row['aj']]
        return df

    def __drop_angles(self,
                      old_new_dict: dict[int, int],  # Of old and updated index
                      OM_index: list[int],  # Index of the OM atoms
                      df: pd.DataFrame  # Angles df of amino
                      ) -> pd.DataFrame:
        """drop and update angles df"""
        df_ = df.copy()
        for item, row in df_.iterrows():
            if row['ai'] in OM_index or \
               row['aj'] in OM_index or \
               row['ak'] in OM_index:
                df.drop(axis=0, index=item, inplace=True)
        del df_
        df_ = df.copy()
        for item, row in df_.iterrows():
            if row['ai'] >= np.min(OM_index):
                df.at[item, 'ai'] = old_new_dict[row['ai']]
            if row['aj'] >= np.min(OM_index):
                df.at[item, 'aj'] = old_new_dict[row['aj']]
            if row['ak'] >= np.min(OM_index):
                df.at[item, 'ak'] = old_new_dict[row['ak']]
        return df

    def __drop_dihedrals(self,
                         old_new_dict: dict[int, int],  # Of old and updated id
                         OM_index: list[int],  # Index of the OM atoms
                         df: pd.DataFrame  # Dihedrals df of amino
                         ) -> pd.DataFrame:
        """drop and update dihedrals df"""
        df_ = df.copy()
        for item, row in df_.iterrows():
            if row['ai'] in OM_index or \
               row['aj'] in OM_index or \
               row['ak'] in OM_index or \
               row['ah'] in OM_index:
                df.drop(axis=0, index=item, inplace=True)
        del df_
        df_ = df.copy()
        for item, row in df_.iterrows():
            if row['ai'] >= np.min(OM_index):
                df.at[item, 'ai'] = old_new_dict[row['ai']]
            if row['aj'] >= np.min(OM_index):
                df.at[item, 'aj'] = old_new_dict[row['aj']]
            if row['ak'] >= np.min(OM_index):
                df.at[item, 'ak'] = old_new_dict[row['ak']]
            if row['ah'] >= np.min(OM_index):
                df.at[item, 'ah'] = old_new_dict[row['ah']]
        return df

    def __set_OM_info(self,
                      Atoms_df: pd.DataFrame,  # Atoms of the amino
                      si_row: pd.DataFrame,  # One row of si df
                      OM_xyz: pd.DataFrame,  # XYZ info of OM atoms for amino
                      OM_index: list[int]  # Index of the OM atoms in amino
                      ) -> pd.DataFrame:
        """set info for OM in the amino df"""

        df: pd.DataFrame = Atoms_df.copy()
        column: list[str]  # Columns of the df to replace informations
        column = ['atom_id', 'mol', 'typ', 'x', 'y', 'z', 'rho', 'azimuth',
                  'polar']
        for j, k in zip(si_row['OM_list'], OM_index):
            OM_row: pd.DataFrame = OM_xyz[OM_xyz['atom_id'] == j]
            for col in column:
                df.at[k, col] = OM_row[col][j]
            del OM_row
        return df
