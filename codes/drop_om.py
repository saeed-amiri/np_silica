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
        amino_OM: int = 3  # Number of OM in amino data file
        OM_index: list[int]  # Index of OM n amino df
        OM_index = [item for item in
                    amino.Atoms_df[amino.Atoms_df['name'] == 'OM'].index]
        amino.Atoms_df = self.__set_OM_info(amino.Atoms_df,
                                            si_row,
                                            OM_xyz,
                                            OM_index)
        print(amino.Atoms_df)
        df = self.__drop_om_atoms(del_OM,
                                  OM_index,
                                  si_row,
                                  amino.Atoms_df.copy())
        # self.__update_atom_ind(df)

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
        atom_df.reset_index(inplace=True)
        atom_df.rename(columns={'index': 'update_ind'}, inplace=True)
        atom_df.index += 1
        return atom_df

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
