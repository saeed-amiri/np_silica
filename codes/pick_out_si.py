import typing
import numpy as np
import pandas as pd
from colors_text import TextColor as bcolors


class Doc:
    """to pick selected Si (based on radius and non-body) based on the
    angles (theta and phi). To avoid close grafting chains on the
    surface of the nanoparticles."""


class PickSi:
    """label each si based on thier azimuth and polar angles"""
    def __init__(self,
                 Si_df: pd.DataFrame  # All the non-body Si in the radius
                 ) -> None:
        self.__set_lables(Si_df)

    def __set_lables(self,
                     Si_df: pd.DataFrame  # All the non-body Si in the radius
                     ) -> pd.DataFrame:
        """adding label to all the Si in the dataframe"""
