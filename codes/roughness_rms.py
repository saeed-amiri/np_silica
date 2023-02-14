import numpy as np
import pandas as pd
import static_info as stinfo
import read_lmp_data as rdlmp
from colors_text import TextColor as bcolors


class Doc:
    """calculate roughness of the nanoparticles
    Input: Read data: the DataFrames with Lammps format of Full atoms
        format,
    Output:
        Root mean square roughness RMS or rq
    """


class Roughness:
    """get RMS"""
    def __init__(self,
                 np_data  # Data of nanoparticles in LAMMPS format
                 ) -> None:
        self.__roughness(np_data)

    def __roughness(self,
                    np_data  # Data of nanoparticles in LAMMPS format
                    ) -> None:
        """Calculate the roughness by finding atoms on the shell"""
        df: pd.DataFrame = self.__mk_rho(np_data.Atoms_df)
        df = self.__select_atoms(df)
        self.rms: float = self.__get_roughness(df)
        self.__print_info()
        del df

    def __print_info(self) -> None:
        """print info on std"""
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\tRoot mean square roughness for shell of size '
              f'"{stinfo.Constants.Shell_radius}" is "{self.rms: .3f}"'
              f'{bcolors.ENDC}')

    def __get_roughness(self,
                        df: pd.DataFrame  # Atoms df with select atoms
                        ) -> float:
        """calculate the roughness of the atoms inside the shell"""
        rms: list[float] = list(df['rho']*df['rho'])
        return np.sqrt(np.mean(rms))

    def __select_atoms(self,
                       df: pd.DataFrame  # Select atoms inside the shell
                       ) -> pd.DataFrame:
        """return atoms which ate inside the shell"""
        return df[df['rho'] >= 0]

    def __mk_rho(self,
                 Atoms_df: pd.DataFrame  # Atoms info in Lammps format
                 ) -> pd.DataFrame:
        """calculate and add rho column (distance from origin) of all
        the atom"""
        df: pd.DataFrame = Atoms_df.copy()
        df['rho']: list[float] = [-1 for _ in df['atom_id']]
        for item, row in Atoms_df.iterrows():
            df.at[item, 'rho'] = self.__get_rho(row['x'], row['y'], row['z'])
        df = self.__apply_shell(df)
        return df

    def __apply_shell(self,
                      Atoms_df: pd.DataFrame  # Atoms df with rho
                      ) -> pd.DataFrame:
        """cut the out the radius of the shell from rho to find atoms
        inside the shell"""
        inner_r: float  # min of the shell (inner radius)
        inner_r = np.max(Atoms_df['rho'] - stinfo.Constants.Shell_radius)
        df: pd.DataFrame = Atoms_df.copy()
        df['rho'] -= inner_r
        return df

    def __get_rho(self,
                  x: float,  # X component of the atoms
                  y: float,  # y component of the atoms
                  z: float,  # z component of the atoms
                  ) -> float:
        """return the rho of the atom"""
        return np.sqrt(x*x + y*y + z*z)
