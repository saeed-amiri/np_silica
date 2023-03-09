"""setting masses section for silanized data file before writing it to
a file.
These data are needed for converting the LAMMPS data file to PDB."""

import typing
import pandas as pd
import static_info as stinfo
from colors_text import TextColor as bcolors


class SetMasses:
    """get info and set the name section"""
    def __init__(self,
                 masses_df: pd.DataFrame  # The masses df after combination
                 ) -> None:
        self.df_masses: pd.DataFrame = self.set_masses(masses_df)
        self.print_info()

    def set_masses(self,
                   masses_df: pd.DataFrame  # The masses df after combination
                   ) -> pd.DataFrame:
        """update the combined masses section"""
        df_c: pd.DataFrame = masses_df.copy()
        for item, row in masses_df.iterrows():
            info: dict[str, typing.Any] = stinfo.PdbMass.ATOMS.get(row['name'])
            info_list: list[typing.Any]
            new_name: str  # To replace the name section with complete info
            info_list = [info.get("Atoms_names"),
                         info.get("Residue"),
                         info.get("Element_symbol"),
                         info.get("RECORD"),
                         info.get("ff_type")]
            new_name = ' '.join(info_list)
            df_c.at[item, 'name'] = new_name.strip('"')
            del new_name
        return df_c

    def print_info(self) -> None:
        """print some info"""
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              '\tUpdateing the Masses section in the silanized output file'
              f'{bcolors.ENDC}\n')
