"""combine all the data given into the system
It combines water_box with the nanoparticle, putting their origin at
zero, for both thier origin must already be at zero.
"""

import numpy as np
import pandas as pd


class MergeAll:
    """get silanized np and water_box and merge all the infos together
    and return in LAMMPS writting format"""
    def __init__(self,
                 water_box,  # Water box all data in LAMMPS full atoms
                 nano_p  # Nano particles all data in LAMMPS full atoms
                 ) -> None:
        self.merge_all(water_box, nano_p)

    def merge_all(self,
                  water_box,  # Water box all data in LAMMPS full atoms
                  nano_p  # Nano particles all data in LAMMPS full atoms
                  ) -> None:
        # First update types
        self.__update_types(water_box, nano_p)

    def __update_types(self,
                       water_box,  # Water box all data in LAMMPS full atoms
                       nano_p  # Nano particles all data in LAMMPS full atoms
                       ) -> None:
        """update atoms, bonds, and angles types of the water_box"""
        # Update atoms type
        water_atoms: pd.DataFrame  # Water atoms with updated atom type
        water_atoms = self.__update_all_types(water_box.Atoms_df,
                                              nano_p.Atoms_df)
        water_bonds = self.__update_all_types(water_box.Bonds_df,
                                              nano_p.Bonds_df)
        water_angles = self.__update_all_types(water_box.Angles_df,
                                               nano_p.Angles_df)

    def __update_all_types(self,
                           water_df: pd.DataFrame,  # Atoms of water
                           nano_p_df: pd.DataFrame  # Atoms of nano_p
                           ) -> pd.DataFrame:
        """update the atoms type"""
        level_up: int = np.max(nano_p_df['typ'])
        df_c: pd.DataFrame = water_df.copy()
        df_c['typ'] += level_up
        return df_c
