"""sanity check!
    Check the length of all the bonds.
    There were some strange long bonds in the data file sometimes,
    This could be because of the random selection of dropping Si from
    silanization or maybe because of the python dataframe selection
    algorithm"""


import numpy as np
import pandas as pd
import read_lmp_data as rdlmp


class CheckBond:
    """calculate all length of all the bonds"""
    def __init__(self,
                 data: rdlmp.ReadData  # In the form of the LAMMPS data
                 ) -> None:
        self.__check_bonds(data)

    def __check_bonds(self,
                      data  # In the form of the LAMMPS data
                      ) -> None:
        """get the all the bonds"""
        Atoms: pd.DataFrame = data.Atoms_df.copy()  # Atoms_df of the data
        Bonds: pd.DataFrame = data.Bonds_df.copy()  # Bonds_df of the data
        Bonds['bond_length']: list[float]  # Lenght of the bonds
        for item, row in data.Bonds_df.iterrows():
            a_ii = row['ai']
            a_jj = row['aj']
            a_i = Atoms.iloc[a_ii-1]
            a_j = Atoms.iloc[a_jj-1]
            if a_ii != a_i['atom_id']:
                print('a_i error!', item)
                print(row)
                print(a_i)
            if a_jj != a_j['atom_id']:
                print('a_j error!', item)
                print(row)
                print(a_j)

    def __calc_len(self,
                   a_i: pd.DataFrame,  # atom a_i in the bond
                   a_j: pd.DataFrame  # atom a_j in the bond
                   ) -> float:
        """claculate the length of bond between atom a_i and a_j"""
        p1 = np.array([a_i['x'], a_i['y'], a_i['z']])
        p2 = np.array([a_j['x'], a_j['y'], a_j['z']])
        return np.linalg.norm(p1 - p2)
