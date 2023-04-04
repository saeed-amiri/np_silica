"""tools for use in different modules:
get_radius: find the the radius of the silanized or pure nano partilces
"""

import os
import sys
import typing
import numpy as np
import pandas as pd
import static_info as stinfo
from colors_text import TextColor as bcolors


def get_radius(atoms_df: pd.DataFrame  # Atoms infos, must contain x,y,z column
               ) -> float:
    """return radius of the nanoparticle, in simple way by finding the
    max in x, y, and z direction"""
    x_max: float = atoms_df['x'].abs().max()
    y_max: float = atoms_df['y'].abs().max()
    z_max: float = atoms_df['z'].abs().max()
    return np.max([x_max, y_max, z_max])


def com_to_zero(atoms_df: pd.DataFrame  # Atoms df
                ) -> pd.DataFrame:
    """set the center of mass to zero"""
    print(f'\n{bcolors.OKCYAN}'
          '\tMove the center of mass to zero'
          f'{bcolors.ENDC}')
    x_cm: float = np.average(atoms_df['x'])
    y_cm: float = np.average(atoms_df['y'])
    z_cm: float = np.average(atoms_df['z'])
    df_c: pd.DataFrame = atoms_df.copy()
    df_c['x'] -= x_cm
    df_c['y'] -= y_cm
    df_c['z'] -= z_cm
    return df_c


def drop_duplicate(l_to_set: list[typing.Any]
                   ) -> list[typing.Any]:
    """drop duplicated item with keeping order"""
    seen: set[str] = set()
    seen_add = seen.add
    return [x for x in l_to_set if not (x in seen or seen_add(x))]


def check_file(c_file: str,  # Name of the out put of PACKMOL
               delete: bool = False  # Keep the file or not
               ) -> int:
    """check if a file exist, if delete"""
    pack_flag: int = -1  # Check if exist
    if os.path.isfile(c_file):
        print(f'{bcolors.CAUTION}\t'
              f'`{sys.modules[__name__]}:`\n'
              f'\t\tAn old "{c_file}" exists'
              f'{bcolors.ENDC}')
        if delete:
            print(f'{bcolors.CAUTION}'
                  f'\t\tThe old "{c_file}" is deleted'
                  f'{bcolors.ENDC}')
            os.remove(c_file)
    else:
        if os.path.isfile(c_file):
            pack_flag = 0
    return pack_flag


def oil_depth(radius: float  # Radius of the NP
              ) -> float:
    """calculate and return the the depth of oil phase `h` in system
    box"""
    angle_rad: float  # Contact angle in radian
    angle_rad = np.radians(stinfo.Hydration.CONATCT_ANGLE)
    return radius * np.tan(angle_rad/2)
