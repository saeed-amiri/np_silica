import random
import numpy as np
import pandas as pd
import static_info as stinfo
from colors_text import TextColor as bcolors


class Doc:
    """Sparse the Grafting on the top of the nanoparticles, randomly
    or with any other method.
    Right now it only uses the random sparsing!"""


class PickSi:
    """sparse the nanoparticle"""
    def __init__(self,
                 Si_df: pd.DataFrame,  # All the non-body Si in the radius
                 diameter: float  # The diameter of the Nanoparticles
                 ) -> None:
        self.__method: str = 'random'  # random or area
        self.Si_df: pd.DataFrame = self.__set_coverage(Si_df, diameter)

    def __set_coverage(self,
                       Si_df: pd.DataFrame,  # All the non-body Si in the radiu
                       diameter: float  # The diameter of the Nanoparticles
                       ) -> pd.DataFrame:
        """select the method of randomly setting the coverage by
        eliminating some of the Si from the data frame or based on the
        area method."""
        si_coverage: float = self.__get_coverage(Si_df, diameter)
        si_de_num: int = self.__get_si_num(diameter)
        if si_coverage <= stinfo.Constants.Coverage:
            print(f'\n{bcolors.WARNING}{self.__class__.__name__}:'
                  f' ({self.__module__})\n'
                  f'\tGrafting all the Si gives "{si_coverage:.3f}" '
                  f'less or equal to the desire coverage '
                  f'"{stinfo.Constants.Coverage:.3f}"! Returns'
                  f'{bcolors.ENDC}')
        else:
            if self.__method == 'random':
                Si_df = self.__random_sparse(Si_df, si_de_num, si_coverage)
            else:
                Si_df = self.__set_lables(Si_df, diameter)
        return Si_df

    def __random_sparse(self,
                        Si_df: pd.DataFrame,  # All the non-body Si in radius
                        si_de_num: int,  # Number of deleted Si
                        si_coverage: float  # Availabel coverage
                        ) -> pd.DataFrame:
        """sparse the Si randomly"""
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}:'
              f' ({self.__module__})\n'
              f'\tGrafting all the Si group gives coverage of '
              f'"{si_coverage:.5f}"\n'
              f'\tDroping "{len(Si_df)-si_de_num}" Si to get the coverage'
              f' of "{stinfo.Constants.Coverage}"'
              f'{bcolors.ENDC}')
        no_od_si: list[int]  # All the Si which are not bonded to an OD
        no_od_si = self.__get_ox_list(Si_df)
        drop_si: list[int]  # Si id to drop from the Si_df
        drop_si = random.sample(no_od_si, len(Si_df)-si_de_num)
        df: pd.DataFrame = Si_df.copy()
        for item in Si_df['atom_id']:
            if item in drop_si:
                df.drop(axis=0, index=item, inplace=True)
        return df

    def __get_ox_list(self,
                      Si_df: pd.DataFrame  # Si df
                      ) -> list[int]:
        """return atom_id of the Si which thier Ox is not OD
        To make sure all the OD will be replaced with Amine chain"""
        return [item for item in Si_df['atom_id'] if
                Si_df['Ox_drop_name'][item] != 'OD']

    def __set_lables(self,
                     Si_df: pd.DataFrame,  # All the non-body Si in the radius
                     diameter: float  # The diameter of the Nanoparticles
                     ) -> pd.DataFrame:
        """adding label to all the Si in the dataframe"""
        N: int = 112
        counter: int = 0
        az_list: list[int] = []
        for _, row in Si_df.iterrows():
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
                     ) -> int:
        """return number of Si for desire coverage of chains on the np"""
        return int(np.floor(
                   np.pi*stinfo.Constants.Coverage*diameter*diameter)) + 1
