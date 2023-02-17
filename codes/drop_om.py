import numpy as np
import pandas as pd
import static_info as stinfo
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
        del df

        bonds_df: pd.DataFrame = amino.Bonds_df.copy()
        self.Bonds_df = self.__drop_bonds(old_new_dict,
                                          bonds_df,
                                          self.Atoms_df)

        angles_df: pd.DataFrame = amino.Angles_df.copy()
        self.Angles_df = self.__drop_angles(old_new_dict,
                                            OM_index,
                                            angles_df,
                                            self.Atoms_df)

        dihedrals_df: pd.DataFrame = amino.Dihedrals_df.copy()
        self.Dihedrals_df = self.__drop_dihedrals(old_new_dict,
                                                  OM_index,
                                                  dihedrals_df)

    def __drop_replicate_boandi(self,
                                df: pd.DataFrame,  # Data to check
                                n_list: list[str]  # Name of the set
                                ) -> pd.DataFrame:
        """check if there is boandi which all in ['SI', 'OM']"""
        df['name'] = n_list
        df_ = df.copy()
        NP_LIST: list[str]  # Name of the amino's root in nanoparticles
        NP_LIST = [stinfo.Constants.SI_amino, stinfo.Constants.OM_amino]
        for item, row in df_.iterrows():
            names = row['name'].split('_')
            if set(names).issubset(NP_LIST):
                df.drop(axis=0, index=item, inplace=True)
        df.drop(axis=1, columns=['name'], inplace=True)
        return df

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
        # Check duplicacy: in case there is duplication after updating
        atom_df = self.__OM_duplicacy_check(atom_df)
        atom_df.reset_index(inplace=True)
        # Keep the index of all atoms before drop
        atom_df.rename(columns={'index': 'undrop_ind'}, inplace=True)
        atom_df.index += 1
        return atom_df

    def __OM_duplicacy_check(self,
                             atom_df: pd.DataFrame  # Atoms_df of OM after ind
                             ) -> pd.DataFrame:
        """check if there are OM atoms with same atom index and drop
        if so"""
        df = atom_df[atom_df['name'] == 'OM']
        if not df['atom_id'].is_unique:
            atom_df.drop_duplicates(subset=['atom_id', 'name'],
                                    keep='first',
                                    inplace=True)
            print(f'\n{bcolors.CAUTION}{self.__class__.__name__}: '
                  f'({self.__module__})\n'
                  f'\tTaking care of atom_id complication in updating'
                  f' atom_id in the aminopropyl chain{bcolors.ENDC}')
        return atom_df

    def __update_atom_ind(self,
                          df: pd.DataFrame,  # Amino atoms df with removed OM
                          OM_list: list[int]  # Index of the OM in NP
                          ) -> pd.DataFrame:
        """update the index of df after OM was droped"""
        for item, row in df.iterrows():
            if row['name'] != stinfo.Constants.SI_amino:
                if row['name'] != 'OM':
                    df.at[item, 'atom_id'] = item
                    df.at[item, 'old_id'] = item
        if len(df[df['name'] == 'OM']) != len(OM_list):
            print(f'{bcolors.WARNING}{self.__class__.__name__}'
                  f'({self.__module__}):\n'
                  f'\tThere are more OM in the amino list then the bonded'
                  f' OM atoms{bcolors.ENDC}')
        return df

    def __drop_bonds(self,
                     old_new_dict: dict[int, int],  # Of old and updated index
                     df: pd.DataFrame,  # Bonds df of amino
                     atoms_df: pd.DataFrame  # Atoms df
                     ) -> pd.DataFrame:
        """drop and update bonds df"""
        df_ = df.copy()
        for item, row in df_.iterrows():
            try:
                df.at[item, 'ai'] = old_new_dict[row['ai']]
                df.at[item, 'aj'] = old_new_dict[row['aj']]
            except KeyError:
                df.drop(axis=0, index=item, inplace=True)
        names = check_boandi_name(atoms_df, df, ['ai', 'aj'])
        df = self.__drop_replicate_boandi(df, names)
        return df

    def __drop_angles(self,
                      old_new_dict: dict[int, int],  # Of old and updated index
                      OM_index: list[int],  # Index of the OM atoms
                      df: pd.DataFrame,  # Angles df of amino
                      atoms_df: pd.DataFrame  # Atoms df to get the names
                      ) -> pd.DataFrame:
        """drop and update angles df"""
        df_ = df.copy()
        for item, row in df_.iterrows():
            try:
                df.at[item, 'ai'] = old_new_dict[row['ai']]
                df.at[item, 'aj'] = old_new_dict[row['aj']]
                df.at[item, 'ak'] = old_new_dict[row['ak']]
            except KeyError:
                df.drop(axis=0, index=item, inplace=True)
        names = check_boandi_name(atoms_df, df, ['ai', 'aj', 'ak'])
        df = self.__drop_replicate_boandi(df, names)
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


def check_boandi_name(Atoms_df: pd.DataFrame,  # Updated atoms df
                      df: pd.DataFrame,  # The df to make name for
                      a_list: list[str]  # list of the atoms columns: ai, aj
                      ) -> None:
    """check the name of the bonds, angles, dihedrals
        make a name column for the bonds"""
    atom_name: dict[int, str]  # id and name of the atoms
    atom_name = {k: v for k, v in zip(Atoms_df['atom_id'],
                                      Atoms_df['name'])}
    name_list: list[str] = []  # Name of the bo/an/di
    for _, row in df.iterrows():
        names = []
        for a in a_list:
            names.append(atom_name[row[a]])
        name_list.append('_'.join(names))
    df['name'] = name_list
    return name_list
