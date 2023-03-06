"""combine all the data given into the system
It combines water_box with the nanoparticle, putting their origin at
zero, for both thier origin must already be at zero.
"""

import numpy as np
import pandas as pd
from colors_text import TextColor as bcolors


class MergeAll:
    """get silanized np and water_box and merge all the infos together
    and return in LAMMPS writting format"""
    def __init__(self,
                 water_box,  # Water box all data in LAMMPS full atoms
                 nano_p  # Nano particles all data in LAMMPS full atoms
                 ) -> None:
        self.Atoms_df: pd.DataFrame  # Merged df of water and nanoparticles
        self.Bonds_df: pd.DataFrame  # Merged df of water and nanoparticles
        self.Angles_df: pd.DataFrame  # Merged df of water and nanoparticles
        self.Masses_Df: pd.DataFrame  # Merged df of water and nanoparticles
        self.Atoms_df, self.Bonds_df, self.Angles_df, self.Masses_Df = \
            self.merge_all(water_box, nano_p)
        self.print_info()

    def merge_all(self,
                  water_box,  # Water box all data in LAMMPS full atoms
                  nano_p  # Nano particles all data in LAMMPS full atoms
                  ) -> tuple[pd.DataFrame,
                             pd.DataFrame,
                             pd.DataFrame,
                             pd.DataFrame]:
        """update the typs, then merge them together"""
        # First update types
        water_atoms: pd.DataFrame  # Water atoms with updated atom type
        water_bonds: pd.DataFrame  # Water bonds with updated bond type
        water_angles: pd.DataFrame  # Water angles with updated angle type
        water_masses: pd.DataFrame  # Water masses with updated angle type
        water_atoms, water_bonds, water_angles, water_masses = \
            self.__update_types(water_box, nano_p)
        # Update atoms id
        level_up: int = np.max(nano_p.Atoms_df['atom_id'])
        water_atoms = self.__update_atoms_id(level_up, water_atoms, 'atoms')
        water_bonds = self.__update_atoms_id(level_up, water_bonds, 'bonds')
        water_angles = self.__update_atoms_id(level_up, water_angles, 'angles')
        # Combine them
        merge_atoms: pd.DataFrame  # Merged df of water and nanoparticles
        merge_bonds: pd.DataFrame  # Merged df of water and nanoparticles
        merge_angles: pd.DataFrame  # Merged df of water and nanoparticles
        merge_atoms = self.__combine_df(nano_p.Atoms_df, water_atoms, 'atoms')
        merge_bonds = self.__combine_df(nano_p.Bonds_df, water_bonds, 'bonds')
        merge_angles = self.__combine_df(nano_p.Angles_df,
                                         water_angles,
                                         'angles')
        return merge_atoms, merge_bonds, merge_angles, water_masses

    def __combine_df(self,
                     nano_p_df: pd.DataFrame,  # df of the nanoparticle
                     water_df: pd.DataFrame,  # df of the water
                     target: str  # Atoms, Bonds, Angles
                     ) -> pd.DataFrame:
        """combine and return the dataframes"""
        columns: list[str]  # Name of the columns of the df
        if target == 'atoms':
            columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y',
                       'z', 'nx', 'ny', 'nz', 'cmt', 'name']
        elif target == 'bonds':
            columns = ['typ', 'ai', 'aj', 'cmt', 'name']
        elif target == 'angles':
            columns = ['typ', 'ai', 'aj', 'ak']
        w_df: pd.DataFrame = water_df[columns].copy()
        np_df: pd.DataFrame = nano_p_df[columns].copy()
        return pd.concat([np_df, w_df])

    def __update_atoms_id(self,
                          level_up: int,  # Maximum of atom_id in the atoms_df
                          water_df: pd.DataFrame,  # df of the water
                          target: str  # Atoms, Bonds, Angles
                          ) -> pd.DataFrame:
        """updating atom id of water_box by adding maximum atom id of
        nanoparticle to each atom id"""
        df_c: pd.DataFrame = water_df.copy()
        if target == 'atoms':
            df_c['atom_id'] += level_up
        elif target == 'bonds':
            df_c['ai'] += level_up
            df_c['aj'] += level_up
        elif target == 'angles':
            df_c['ai'] += level_up
            df_c['aj'] += level_up
            df_c['ak'] += level_up
        return df_c

    def __update_types(self,
                       water_box,  # Water box all data in LAMMPS full atoms
                       nano_p  # Nano particles all data in LAMMPS full atoms
                       ) -> tuple[pd.DataFrame,
                                  pd.DataFrame,
                                  pd.DataFrame,
                                  pd.DataFrame]:
        """update atoms, bonds, and angles types of the water_box"""
        # Update atoms type
        water_atoms: pd.DataFrame  # Water atoms with updated atom type
        water_bonds: pd.DataFrame  # Water bonds with updated atom type
        water_angles: pd.DataFrame  # Water angles with updated atom type
        water_masses: pd.DataFrame  # Water masses with updated atom type
        water_atoms = self.__update_all_types(water_box.Atoms_df,
                                              nano_p.Atoms_df)
        water_bonds = self.__update_all_types(water_box.Bonds_df,
                                              nano_p.Bonds_df)
        water_angles = self.__update_all_types(water_box.Angles_df,
                                               nano_p.Angles_df)
        water_masses = self.__update_masses_df(water_atoms,
                                               water_box.Masses_df)
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              '\tUpdate all the types in Atoms, Bonds, and Angles of '
              f'"water_box"{bcolors.ENDC}')
        return water_atoms, water_bonds, water_angles, water_masses

    def __update_all_types(self,
                           water_df: pd.DataFrame,  # Atoms of water
                           nano_p_df: pd.DataFrame  # Atoms of nano_p
                           ) -> pd.DataFrame:
        """update the atoms type"""
        level_up: int = np.max(nano_p_df['typ'])
        df_c: pd.DataFrame = water_df.copy()
        df_c['typ'] += level_up
        return df_c

    def __update_masses_df(self,
                           water_atoms: pd.DataFrame,  # Atoms with updated typ
                           water_masses: pd.DataFrame  # Masses df of water
                           ) -> pd.DataFrame:
        """update the masses for water box based on the updated types"""
        h_type: int  # Type of hydrogen in the df
        o_type: int  # Type of oxygen in the df
        h_type = list(
                      set(water_atoms[water_atoms['name'] == 'H']['typ'])
                      )[0]
        o_type = list(
                      set(water_atoms[water_atoms['name'] == 'O']['typ'])
                      )[0]
        df_c: pd.DataFrame = water_masses.copy()
        df_c.at[df_c[df_c['name'] == 'H'].index[0], 'typ'] = h_type
        df_c.at[df_c[df_c['name'] == 'O'].index[0], 'typ'] = o_type
        return df_c

    def print_info(self) -> None:
        """print infos"""
        print(f'{bcolors.OKCYAN}\tNanoparticle and the water_box are '
              f'merged together{bcolors.ENDC}')
