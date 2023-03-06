"""tools for use in different modules:
get_radius: find the the radius of the silanized or pure nano partilces
"""

import pandas as pd
import numpy as np


def get_radius(atoms_df: pd.DataFrame  # Atoms infos, must contain x,y,z column
               ) -> float:
    """return radius of the nanoparticle, in simple way by finding the
    max in x, y, and z direction"""
    x_max: float = atoms_df['x'].abs().max()
    y_max: float = atoms_df['y'].abs().max()
    z_max: float = atoms_df['z'].abs().max()
    return np.max([x_max, y_max, z_max])
