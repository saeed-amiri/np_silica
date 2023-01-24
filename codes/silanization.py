import sys
import numpy as np
import pandas as pd
import read_lmp_data as rdlmp
import update_coords as upcord
from colors_text import TextColor as bcolors


class Doc:
    """read data
        select the gropus to replace
        write the output file
    """


class GetAmino(rdlmp.ReadData):
    """read the main Aminopropyle coordinates and put the Si position
    to zero"""
    def __init__(self) -> None:
        fname: str = '/scratch/saeed/MyScripts/np_silica/data/aminopropyl.data'
        super().__init__(fname)
        self.Si = 'Si'
        self.__to_origin()

    def __to_origin(self) -> pd.DataFrame:
        """put the coordinate of Si in Aminopropyle to zero"""
        df: pd.DataFrame = self.Atoms_df.copy()
        print(f'{bcolors.OKCYAN}\tMove Aminopropyle [Si] to origin\n'
              f'{bcolors.ENDC}')
        x_si: float = self.Atoms_df[self.Atoms_df['name'] == self.Si]['x'][1]
        y_si: float = self.Atoms_df[self.Atoms_df['name'] == self.Si]['y'][1]
        z_si: float = self.Atoms_df[self.Atoms_df['name'] == self.Si]['z'][1]
        df['x'] -= x_si
        df['y'] -= y_si
        df['z'] -= z_si
        return df


class PrepareAmino:
    """get the data for selected Silicons:
        first rotate aminopropyle for each of one the Si group
        replace the Si of each aminopropyle with the with the index of
        the selected group and updated Si
        """
    def __init__(self,
                 si_df: pd.DataFrame,  # DF of Si groups with rotation angles
                 amino: GetAmino  # Information about one aminos
                 ) -> None:
        """apply the position update to the aminos"""
        self.__update_aminos(si_df, amino)

    def __update_aminos(self,
                        si_df: pd.DataFrame,  # Si groups with rotation angles
                        amino: GetAmino  # Information about one aminos
                        ) -> None:
        """do"""
        si_df = self.__order_si_df(si_df)
        for item, row in si_df.iterrows():
            self.__rotate_amino(amino, row)

    def __rotate_aminos(self,
                        amino: GetAmino,  # Amino information,
                        row: pd.DataFrame  # A row of the dataframe
                        ) -> pd.DataFrame:
        """rotate each amino based on the azimuuth polar angle in the row"""

    def __order_si_df(self,
                      si_df: pd.DataFrame  # Si groups with rotation angles
                      ) -> pd.DataFrame:
        # reindex the si_df to get orders
        si_df.reset_index(inplace=True)
        si_df.index += 1
        si_df.drop(columns=['index'], axis=1, inplace=True)
        return si_df


if __name__ == '__main__':
    fname = sys.argv[1]
    update = upcord.UpdateCoords(fname)
    amino = GetAmino()
    up_aminos = PrepareAmino(update.Si_df, amino)
