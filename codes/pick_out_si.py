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
                 Si_df: pd.DataFrame,  # All the non-body Si in the radius
                 diameter: float,  # The diameter of the Nanoparticles
                 Atoms: pd.DataFrame  # All the atoms in the nanoparticles
                 ) -> None:
        self.__de_coverage: float = 3.0  # The desire coverage
        self.__method: str = 'random'  # random or area
        self.__set_coverage(Si_df, diameter, Atoms)

    def __set_coverage(self,
                       Si_df: pd.DataFrame,  # All the non-body Si in the radiu
                       diameter: float,  # The diameter of the Nanoparticles
                       Atoms: pd.DataFrame  # All the atoms in the nanoparticle
                       ) -> pd.DataFrame:
        """select the method of randomly setting the coverage by
        eliminating some of the Si from the data frame or based on the
        area method."""
        si_coverage: float = self.__get_coverage(Si_df, diameter)
        si_de_num: float = self.__get_si_num(diameter)
        if si_coverage <= si_de_num:
            print(f'\n{bcolors.WARNING}{self.__class__.__name__}:'
                  f' ({self.__module__})\n'
                  f'\tGrafting all the Si gives ({si_coverage:.4f}) '
                  f'less or equal to the desire coverage '
                  f'({self.__de_coverage:.4f})! Returns.'
                  f'{bcolors.ENDC}')
        else:
            if self.__method == 'random':
                print('randomly  sparse the coverage')
                self.__random_sparse(Si_df, diameter, Atoms)
            else:
                self.__set_lables(Si_df, diameter)

    def __random_sparse(self,
                        Si_df: pd.DataFrame,  # All the non-body Si in radius
                        diameter: float,  # The diameter of the Nanoparticles
                        Atoms: pd.DataFrame  # All the atoms in the nanoparticl
                        ) -> pd.DataFrame:
        """sparse the Si randomly"""
        si_coverage: float = self.__get_coverage(Si_df, diameter)
        si_de_num: float = self.__get_si_num(diameter)
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}:'
              f' ({self.__module__})\n'
              f'\tGrafting all the Si group gives coverage of '
              f'"{si_coverage:.5f}"\n'
              f'\tDroping "{len(Si_df)-si_de_num}" Si to get the coverage'
              f' of "{self.__de_coverage}"'
              f'{bcolors.ENDC}')
        print(-si_de_num+len(Si_df))

    def __set_lables(self,
                     Si_df: pd.DataFrame,  # All the non-body Si in the radius
                     diameter: float  # The diameter of the Nanoparticles
                     ) -> pd.DataFrame:
        """adding label to all the Si in the dataframe"""
        si_coverage: float = self.__get_coverage(Si_df, diameter)
        N: int = 112
        counter: int = 0
        az_list: list[float] = []
        for item, row in Si_df.iterrows():
            label_flag: bool = False
            theta_i: float = -np.pi
            for i in range(N):
                theta = 2*np.pi/N + theta_i
                if row['azimuth'] < theta and row['azimuth'] >= theta_i:
                    label_flag = True
                    counter += 1
                    az_list.append(i)
                theta_i = theta
                if label_flag:
                    break
        Si_df['az_lable'] = az_list
        for i in set(az_list):
            df_ = Si_df[Si_df['az_lable'] == i]
            del df_

    def __get_coverage(self,
                       Si_df: pd.DataFrame,  # All the non-body Si in the radiu
                       diameter: float  # The diameter of the Nanoparticles
                       ) -> float:
        """return the coverage of the selected Si are on the np"""
        return len(Si_df)/(np.pi*diameter*diameter)

    def __get_si_num(self,
                     diameter: float  # The diameter of the Nanoparticles
                     ) -> float:
        """return number of Si for desire coverage of chains on the np"""
        return int(np.pi*self.__de_coverage*diameter*diameter) + 1