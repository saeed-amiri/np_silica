import typing
import numpy as np
import pandas as pd
import read_lmp_data as rdlmp
from colors_text import TextColor as bcolors


class GetData(rdlmp.ReadData):
    """read data file in LAMMPS format"""
    def __init__(self,
                 fname: str  # Name of the data file
                 ) -> None:
        super().__init__(fname)


class GetSiGroups:
    """get silinol groups to add cahin to them"""
    def __init__(self,
                 Atoms: pd.DataFrame,  # Atoms df in form of lammps fullatom
                 Sigroup: list[typing.Any],  # Name | index to select group[Si]
                 fraction: float = 1  # Fraction of Silica to remove
                 ) -> None:
        self.df_Si = self.__get_silica(Atoms, Sigroup)
        self.Si_delete: list[int] = [item for item in self.df_Si['atom_id']]

    def __get_silica(self,
                     Atoms: pd.DataFrame,  # Atoms df in the lammps fullatom
                     Sigroup: list[typing.Any],  # Name | index of groups[Si]
                     ) -> pd.DataFrame:
        """Get index or name of atoms"""
        Si_list: list[pd.DataFrame] = []  # df with asked Si groups
        for item in Sigroup:
            Si_list.append(Atoms[Atoms['name'] == item])
        # DF of all the Si to replace:
        df: pd.DataFrame = pd.concat(Si_list)
        # Drop unwanted columns:
        df = self.__drop_cols(df)
        # Check if contain atoms:
        df = self.__check_si_df(df, Sigroup)
        # Get the nano particle radius:
        max_radius: float = self.__get_sys_radius(Atoms)
        # Get spherical coords for Si atoms:
        df = self.__get_sphere_coord(df)
        # Drop atom outside the shell:
        df = self.__apply_radius(df, max_radius)
        return df

    def __drop_cols(self,
                    df: pd.DataFrame,  # Dataframe from selected Si atoms
                    ) -> pd.DataFrame:
        """drop some unwanted columns"""
        df.drop(axis=1, columns=['nx', 'ny', 'nz', 'cmt', 'b_name'],
                inplace=True)
        return df

    def __check_si_df(self,
                      df: pd.DataFrame,  # Dataframe from selected Si atoms
                      Sigroup: list[typing.Any],  # Name | index of groups[Si]
                      ) -> pd.DataFrame:
        """check the data frame"""
        if len(df) > 0:
            print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
                  f'({self.__module__})\n'
                  f'\tThere are: {len(df)} atoms with selected '
                  f'atoms name [Si] in the file{bcolors.ENDC}')
        else:
            exit(f'{bcolors.FAIL}Error! The selcted atom group `{Sigroup}`'
                 f'does not exist!\n{bcolors.ENDC}')
        return df

    def __apply_radius(self,
                       df: pd.DataFrame,  # Checked df of selected Si
                       radius: float  # Max radius to check the shell
                       ) -> pd.DataFrame:
        """keep the Si on the shell"""
        df = df[(df['rho'] >= radius - 6)]
        print(f'{bcolors.OKBLUE}\tThere are: {len(df)} Si atoms in the '
              f'choosen area of the system, Max_radius = {radius:.3f}'
              f'{bcolors.ENDC}')
        return df

    def __get_sys_radius(self,
                         Atoms: pd.DataFrame,  # Atoms in lammps full atom
                         ) -> float:
        """return the radius of the nano-particles"""
        x_max: float = Atoms['x'].abs().max()
        y_max: float = Atoms['y'].abs().max()
        z_max: float = Atoms['z'].abs().max()
        return np.max([x_max, y_max, z_max])

    def __get_sphere_coord(self,
                           df: pd.DataFrame  # The selected Si group
                           ) -> pd.DataFrame:
        """find and set the azimuth and polar angle of the atoms"""
        rho: list[float] = []  # To add to dataframe
        azimuth: list[float] = []  # To add to dataframe
        polar: list[float] = []  # To add to dataframe
        for _, row in df.iterrows():
            x = row['x']
            y = row['y']
            z = row['z']
            i_rho: float = np.sqrt(x*x + y*y + z*z)
            i_azimuth: float = np.arctan2(x, y)
            i_polar: float = np.arccos(z/i_rho)
            rho.append(i_rho)
            azimuth.append(i_azimuth)
            polar.append(i_polar)
        df['rho'] = rho
        df['azimuth'] = azimuth
        df['polar'] = polar
        return df


class GetOmGroups:
    """Get OM groups bonded to the selcted Silicon groups in the radius
     to replace with the ones in the aminopropyles"""
    def __init__(self,
                 silica: rdlmp.ReadData,  # Atoms in form of lammps
                 Si_df: pd.DataFrame,  # Silicon in the shell
                 OMgroup: list[str]  # Name of the OM oxygen to get
                 ) -> None:
        self.replace_oxy: dict[int, list[int]]  # Si with bonded O to replace
        self.replace_oxy = self.__get_OMgroups(silica,
                                               Si_df,
                                               OMgroup)
        self.OM_list: list[int] = self.__get_OM_list()  # All OM atoms
        self.Si_df: pd.DataFrame  # Si df with droped unbonded Si
        self.Si_OM: list[int]  # Si bonded to OM
        self.Si_df = self.__update_si_df(Si_df)
        self.Si_OM = self.__get_Si_OM()

    def __get_OM_list(self) -> list[int]:
        """return list of OM atoms"""
        ll: list[int] = []  # Of all OM
        for _, v in self.replace_oxy.items():
            ll.extend(v)
        return ll

    def __get_OMgroups(self,
                       silica: rdlmp.ReadData,  # All df in form of lammps
                       Si_df: pd.DataFrame,  # DF of Si in the radius
                       OMgroup: list[str]  # Name of the OM oxygen to get
                       ) -> pd.DataFrame:
        df_om: pd.DataFrame = self.__find_OMgroups(silica.Atoms_df, OMgroup)
        replace_o_dict: dict[int, list[int]]  # Get OM of each selected Si
        replace_o_dict = self.__get_OmSi(silica.Bonds_df, df_om, Si_df)
        return replace_o_dict

    def __find_OMgroups(self,
                        Atoms: pd.DataFrame,  # Atoms df in form of lammps
                        OMgroup: list[str]  # Name of the OM oxygen to get
                        ) -> pd.DataFrame:
        """find the OM atoms"""
        OM_list: list[pd.DataFrame] = []  # Of OM groups
        for item in OMgroup:
            OM_list.append(Atoms[Atoms['name'] == item])
        df: pd.DataFrame = pd.concat(OM_list)  # All the OM atoms
        df = self.__drop_cols(df)
        return df

    def __get_OmSi(self,
                   bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                   df_om: pd.DataFrame,  # DF with selected Oxygen
                   Si_df: pd.DataFrame  # DF of Si in the radius
                   ) -> dict[int, list[int]]:  # index of the O to delete
        # get bonds
        all_o = [item for item in df_om['atom_id']]
        replace_o_dict: dict[int, list[int]] = {}  # Get OM of each selected Si
        replace_o_dict = {item: [] for item in Si_df['atom_id']}
        for _, row in bonds_df.iterrows():
            if row['ai'] in Si_df['atom_id']:
                if row['aj'] in all_o:
                    replace_o_dict[row['ai']].append(row['aj'])
            if row['aj'] in Si_df['atom_id']:
                if row['ai'] in all_o:
                    replace_o_dict[row['aj']].append(row['ai'])
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\tThere are {len(replace_o_dict)} `OM` atoms bonded to the '
              f'selected Si{bcolors.ENDC}')
        return replace_o_dict

    def __update_si_df(self,
                       Si_df: pd.DataFrame  # DF of selected Si
                       ) -> pd.DataFrame:
        """drop the silicons which are not bonded from Si_df and/or
        have more then three OM bonds, which means they are body Si"""
        df: pd.DataFrame = Si_df.copy()
        df['OM_list0']: list[typing.Any]  # Index of OM atoms bonded to the Si
        df['OM_list0'] = [None for _ in self.replace_oxy]
        for item, _ in Si_df.iterrows():
            if item not in self.replace_oxy.keys():
                df.drop(index=[item], axis=0, inplace=True)
            elif len(self.replace_oxy[item]) > 3:
                df.drop(index=[item], axis=0, inplace=True)
            else:
                df.at[item, 'OM_list0'] = self.replace_oxy[item]
        print(f'{bcolors.OKBLUE}'
              f'\tThere are {len(df)} `Si` atoms bonded to less then '
              f'three OM atoms{bcolors.ENDC}\n')
        return df

    def __get_Si_OM(self) -> list[int]:
        """return list of the Si bonded to the OM"""
        return [item for item in self.Si_df['atom_id']]

    def __drop_cols(self,
                    df: pd.DataFrame,  # Dataframe from selected Si atoms
                    ) -> pd.DataFrame:
        """drop some unwanted columns"""
        df.drop(axis=1, columns=['nx', 'ny', 'nz', 'cmt', 'b_name'],
                inplace=True)
        return df


class GetOxGroups:
    """get Oxygen and groups to delete them and update
    data file.
    Oxygen is bonded to the silica and Hydrogen, if any, bonded to the
    oxygens.
    """
    def __init__(self,
                 silica: rdlmp.ReadData,  # Atoms df in form of lammps fullatom
                 Si_OM: list[int],  # With selected group[Si] & OM bonded
                 Si_df: pd.DataFrame,  # All selected Si atoms
                 Ogroup: list[str],  # Name groups[O] to delete
                 fraction: float = 1  # Fraction of to select from, 0<fr<=1
                 ) -> None:
        self.O_delete: list[int]  # All the OD atoms to delete
        self.O_delete, self.bonded_si = self.__get_oxgygen(silica,
                                                           Si_OM,
                                                           Ogroup,
                                                           Si_df)

    def __get_oxgygen(self,
                      silica: rdlmp.ReadData,  # Atoms df in lammps full atom
                      Si_OM:  list[int],  # With selected group[Si] & OM bonded
                      Ogroup:  list[typing.Any],  # Name|index groups[O]
                      Si_df: pd.DataFrame  # All selected Si atoms
                      ) -> list[int]:
        """Find the hydrogen which have bonds with the selected Silicons"""
        O_list: list[pd.DataFrame] = []  # df of all Oxygen groups
        Atoms = silica.Atoms_df.copy()
        for item in Ogroup:
            O_list.append(Atoms[Atoms['name'] == item])
        df_o: pd.DataFrame = pd.concat(O_list)  # All the O atoms with names
        O_delete, bonded_si = self.__get_o_delete(silica.Bonds_df,
                                                  Si_OM,
                                                  df_o,
                                                  Si_df)
        return O_delete, bonded_si

    def __get_o_delete(self,
                       bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                       Si_OM: list[int],  # With selected group[Si] & OM bonded
                       df_o: pd.DataFrame,  # DF with selected Oxygen
                       Si_df: pd.DataFrame  # All selected Si atoms
                       ) -> list[int]:  # index of the O to delete
        # get bonded Ox atoms
        check_dict: dict[int, list[int]]  # Si: [ox bondec]
        check_dict = self.__get_check_dict(bonds_df, Si_OM, df_o)
        # Add a column with all Ox bonded to them
        Si_df = self.__get_ox_column(Si_df, check_dict)
        bonded_si: list[int] = []  # Si bonded to the oxygens
        bonded_O: list[int] = []  # O bonded to the Silica
        for k, v in check_dict.items():
            bonded_si.append(k)
            bonded_O.extend(v)
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\tThere are {len(bonded_O)} `O` atoms bonded to the '
              f'slected Si\n'
              f'\t-> There are {len(bonded_si)} `Si` atoms bonded to the '
              f'selected `O` atoms\n{bcolors.ENDC}')
        return bonded_O, bonded_si

    def __get_check_dict(self,
                         bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                         Si_OM: list[int],  # With selected group[Si] & OM bond
                         df_o: pd.DataFrame,  # DF with selected Oxygen
                         ) -> dict[int, list[int]]:
        """a dict with Si and bonded Ox to them """
        all_o = [item for item in df_o['atom_id']]
        check_dict: dict[int, list[int]]  # to check if all si have the bond Ox
        check_dict = {item: [] for item in Si_OM}
        # Make a dictionary of Si: Oxygens
        for _, row in bonds_df.iterrows():
            if row['ai'] in Si_OM:
                if row['aj'] in all_o:
                    check_dict[row['ai']].append(row['aj'])
            if row['aj'] in Si_OM:
                if row['ai'] in all_o:
                    check_dict[row['aj']].append(row['ai'])
        # Check if there is Si without wanted O groups
        dict_cp = check_dict.copy()
        for k, v in dict_cp.items():
            if not v:
                check_dict.pop(k)
        return check_dict

    def __get_ox_column(self,
                        Si_df: pd.DataFrame,  # All selected Si
                        check_dict: dict[int, list[int]]  # Si: Ox bonded
                        ) -> pd.DataFrame:
        """return the Si_df wiht Ox info column"""
        Si_df['Ox_list']: list[typing.Any]  # Ox atoms bonded to the Si
        Si_df['Ox_list'] = [None for _ in range(len(Si_df))]
        ii: int = 0  # Count the number of Ox bonded to the Si
        iii: int = 0  # Count the number of Ox bonded to the Si
        for item, row in Si_df.iterrows():
            Si_df.at[item, 'Ox_list'] = check_dict[item]
            if len(check_dict[item]) > 1:
                if len(check_dict[item]) == 2:
                    ii += 1
                elif len(check_dict[item]) == 3:
                    iii += 1
        print(f'{bcolors.WARNING}{self.__class__.__name__}: \n'
              f'\tThere are "{ii}" `Si` atoms bonded to '
              f'two, and "{iii}" `Si` atoms bonded to three Ox groups'
              f'{bcolors.ENDC}')
        return Si_df


class GetHyGroups:
    """Find Hydrogen groups bonded to the Oxygen atoms"""
    def __init__(self,
                 silica: rdlmp.ReadData,  # All silica atoms
                 O_delete: list[int],  # Index of O atoms to delete
                 Hgroup: list[str]  # name of the H atoms to check for bonds
                 ) -> None:
        self.H_delete: list[int]  # List of all the H atoms to delete
        self.H_delete = self.__do_hydrogens(silica, O_delete, Hgroup)

    def __do_hydrogens(self,
                       silica: rdlmp.ReadData,  # All silica atoms
                       O_delete: list[int],  # Index of O atoms to delete
                       Hgroup: list[str]  # name of the H atoms to check
                       ) -> list[int]:
        df_H: pd.DataFrame   # All the hydrogens with the selcted type
        df_H = self.__get_hydrogens(silica.Atoms_df, O_delete, Hgroup)
        H_delete: list[int] = self.__H_delete(silica.Bonds_df, O_delete, df_H)
        return H_delete

    def __get_hydrogens(self,
                        silica_atoms: pd.DataFrame,  # All silica atoms
                        O_delete: list[int],  # Index of O atoms to delete
                        Hgroup: list[str]  # name of the H atoms to check for
                        ) -> pd.DataFrame:
        # Get Hydrogen atoms
        H_list: list[pd.DataFrame] = []  # df of all Hydrogens
        for item in Hgroup:
            H_list.append(silica_atoms[silica_atoms['name'] == item])
        return pd.concat(H_list)

    def __H_delete(self,
                   bonds_df: pd.DataFrame,  # All the bonds in silica
                   O_delete: list[int],  # Index of the O atoms
                   df_H: pd.DataFrame  # Df of all the Hydrogens
                   ) -> None:
        all_h = [item for item in df_H['atom_id']]
        delete_list: list[int] = []  # index of H atoms to delete
        for _, row in bonds_df.iterrows():
            if row['ai'] in O_delete or row['aj'] in O_delete:
                if row['ai'] in all_h and row['ai'] not in delete_list:
                    delete_list.append(row['ai'])
                if row['aj'] in all_h and row['aj'] not in delete_list:
                    delete_list.append(row['aj'])
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\tThere are {len(delete_list)} `H` atoms bonded to the '
              f'slected O{bcolors.ENDC}')
        return delete_list
