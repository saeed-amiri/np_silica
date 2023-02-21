"""Getting atoms based on their names:
    Si: a group of silicon atoms on the shell to work with them.
    OM: a group of O atoms bonded to Si to put set the angle and
    torsion with the amino group
    Ox: O atoms to drop and replace with amine chains
    H group: the H atoms bonded to some of the Ox atoms to drop."""


import sys
import typing
import numpy as np
import pandas as pd
import pick_out_si as pickSi
import static_info as stinfo
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
                 si_group: list[typing.Any]  # Name | index to select group[Si]
                 ) -> None:
        self.df_Si = self.__get_silica(Atoms, si_group)
        self.Si_delete: list[int] = list(self.df_Si['atom_id'])

    def __get_silica(self,
                     Atoms: pd.DataFrame,  # Atoms df in the lammps fullatom
                     si_group: list[typing.Any],  # Name | index of groups[Si]
                     ) -> pd.DataFrame:
        """Get index or name of atoms"""
        si_list: list[pd.DataFrame] = []  # df with asked Si groups
        for item in si_group:
            si_list.append(Atoms[Atoms['name'] == item])
        # DF of all the Si to replace:
        df: pd.DataFrame = pd.concat(si_list)
        # Drop unwanted columns:
        df = self.__drop_cols(df)
        # Check if contain atoms:
        df = self.__check_si_df(df, si_group)
        # Check if index of Si is less then 17: the number of atoms in
        # amino atoms, as workaround for similar index problem
        df = self.__check_si_id(df)
        # Get the nano particle radius:
        max_radius: float = self.__get_sys_radius(Atoms)
        # Get spherical coords for Si atoms:
        df = self.__get_sphere_coord(df)
        # Drop atom outside the shell:
        df = self.__apply_radius(df, max_radius)
        return df

    def __check_si_id(self,
                      df: pd.DataFrame  # Silcons in the silica
                      ) -> pd.DataFrame:
        """Check if index of Si is less then 17: the number of atoms in
        amino atoms, as workaround for similar index problem"""
        df_: pd.DataFrame = df[df['atom_id'] > stinfo.Constants.Num_amino]
        print(f'{bcolors.CAUTION}{self.__class__.__name__}:'
              f'({self.__module__}):\n'
              f'\t"{len(df)-len(df_)}" of selected Si atoms are dropped'
              f' due to the small atom_id of Si in df (less then '
              f'{stinfo.Constants.Num_amino})'
              f'{bcolors.ENDC}')
        return df_

    def __drop_cols(self,
                    df: pd.DataFrame,  # Dataframe from selected Si atoms
                    ) -> pd.DataFrame:
        """drop some unwanted columns"""
        df.drop(axis=1, columns=['nx', 'ny', 'nz', 'cmt', 'b_name'],
                inplace=True)
        return df

    def __check_si_df(self,
                      df: pd.DataFrame,  # Dataframe from selected Si atoms
                      si_group: list[typing.Any],  # Name | index of groups[Si]
                      ) -> pd.DataFrame:
        """check the data frame"""
        if len(df) > 0:
            print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
                  f'({self.__module__})\n'
                  f'\tThere are: {len(df)} atoms with selected '
                  f'atoms name [Si] in the file{bcolors.ENDC}')
        else:
            sys.exit(f'{bcolors.FAIL}Error! The selcted atom group '
                     f'`{si_group}` does not exist!\n{bcolors.ENDC}')
        return df

    def __apply_radius(self,
                       df: pd.DataFrame,  # Checked df of selected Si
                       radius: float  # Max radius to check the shell
                       ) -> pd.DataFrame:
        """keep the Si on the shell"""
        # Check logic of the Shell size
        if stinfo.Constants.Shell_radius > np.max(df['rho']):
            print(f'\t{bcolors.CAUTION}Warning: shell size: '
                  f'"{stinfo.Constants.Shell_radius:.3f}" is bigger then '
                  f'maximum radius which is: "{np.max(df["rho"]):.3f}"\n'
                  f'{bcolors.ENDC}')
        df = df[(df['rho'] >= radius - stinfo.Constants.Shell_radius)]
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
                 si_df: pd.DataFrame,  # Silicon in the shell
                 om_group: list[str]  # Name of the OM oxygen to get
                 ) -> None:
        self.replace_oxy: dict[int, list[int]]  # Si with bonded O to replace
        self.replace_oxy_name: dict[int, list[str]]  # Si with bonded OM names
        self.replace_oxy, self.replace_oxy_name = self.__get_omgroups(silica,
                                                                      si_df,
                                                                      om_group)
        self.OM_list: list[int] = self.__get_OM_list()  # All OM atoms
        self.si_df: pd.DataFrame  # Si df with droped unbonded Si
        self.Si_OM: list[int]  # Si bonded to OM
        si_df = self.__update_si_df(si_df)
        si_df = self.__drop_small_id(si_df)
        self.si_df = self.__set_om_names(si_df)
        self.Si_OM = self.__get_si_om()

    def __drop_small_id(self,
                        si_df: pd.DataFrame  # si_df with OM_names
                        ) -> pd.DataFrame:
        """Check if index of OM is less then 17: the number of atoms in
        amino atoms, as workaround for similar index problem"""
        df: pd.DataFrame = si_df.copy()
        for item, row in si_df.iterrows():
            for om in row['OM_replace']:
                if om < stinfo.Constants.Num_amino:
                    try:
                        df.drop(axis=0, index=[item], inplace=True)
                    except KeyError:
                        pass

        print(f'{bcolors.CAUTION}{self.__class__.__name__}:'
              f'({self.__module__}):\n'
              f'\t"{len(si_df)-len(df)}" of selected Si atoms are dropped'
              f' due to the small atom_id of OM atoms in df (less then '
              f'{stinfo.Constants.Num_amino})'
              f'{bcolors.ENDC}')
        return df

    def __get_OM_list(self) -> list[int]:
        """return list of OM atoms"""
        ll: list[int] = []  # Of all OM
        for _, v in self.replace_oxy.items():
            ll.extend(v)
        return ll

    def __get_omgroups(self,
                       silica: rdlmp.ReadData,  # All df in form of lammps
                       si_df: pd.DataFrame,  # DF of Si in the radius
                       om_group: list[str]  # Name of the OM oxygen to get
                       ) -> pd.DataFrame:
        df_om: pd.DataFrame = self.__find_om_groups(silica.Atoms_df, om_group)
        replace_o_dict: dict[int, list[int]]  # Get OM of each selected Si
        name_o_dict: dict[int, list[str]]  # Get the name of OM to distinguish
        replace_o_dict, name_o_dict = self.__get_om_si(silica.Bonds_df,
                                                       df_om,
                                                       si_df)
        return replace_o_dict, name_o_dict

    def __find_om_groups(self,
                         Atoms: pd.DataFrame,  # Atoms df in form of lammps
                         om_group: list[str]  # Name of the OM oxygen to get
                         ) -> pd.DataFrame:
        """find the OM atoms"""
        OM_list: list[pd.DataFrame] = []  # Of OM groups
        for item in om_group:
            OM_list.append(Atoms[Atoms['name'] == item])
        df: pd.DataFrame = pd.concat(OM_list)  # All the OM atoms
        df = self.__drop_cols(df)
        return df

    def __get_om_si(self,
                    bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                    df_om: pd.DataFrame,  # DF with selected Oxygen
                    si_df: pd.DataFrame  # DF of Si in the radius
                    ) -> tuple[dict[int, list[int]], dict[int, list[str]]]:
        """return index & name of the OM atoms to repalce with the ones
        in the nanoparticles"""
        # get bonds
        all_o = [item for item in df_om['atom_id']]
        replace_o_dict: dict[int, list[int]] = {}  # Get OM of each selected Si
        replace_o_dict = {item: [] for item in si_df['atom_id']}
        name_o_dict: dict[int, list[str]]  # Get the name of OM to distinguish
        name_o_dict = {item: [] for item in si_df['atom_id']}
        for _, row in bonds_df.iterrows():
            if row['ai'] in si_df['atom_id']:
                if row['aj'] in all_o:
                    replace_o_dict[row['ai']].append(row['aj'])
                    name: str  # Name of the OM atom
                    name = df_om[df_om.index == row['aj']]['name'][row['aj']]
                    name_o_dict[row['ai']].append(name)
            if row['aj'] in si_df['atom_id']:
                if row['ai'] in all_o:
                    replace_o_dict[row['aj']].append(row['ai'])
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\tThere are {len(replace_o_dict)} `OM` atoms bonded to the '
              f'selected Si{bcolors.ENDC}')
        return replace_o_dict, name_o_dict

    def __update_si_df(self,
                       si_df: pd.DataFrame  # DF of selected Si
                       ) -> pd.DataFrame:
        """drop the silicons which are not bonded from si_df and/or
        have more then three OM bonds, which means they are body Si"""
        df: pd.DataFrame = si_df.copy()
        df['OM_replace']: list[typing.Any]  # Index of OM atoms bonded to the Si
        df['OM_replace'] = [None for _ in self.replace_oxy]
        for item, _ in si_df.iterrows():
            if item not in self.replace_oxy.keys():
                df.drop(index=[item], axis=0, inplace=True)
            elif len(self.replace_oxy[item]) > 3:
                df.drop(index=[item], axis=0, inplace=True)
            else:
                df.at[item, 'OM_replace'] = self.replace_oxy[item]
        print(f'{bcolors.OKBLUE}'
              f'\tThere are {len(df)} `Si` atoms bonded to less then '
              f'three OM atoms{bcolors.ENDC}\n')
        return df

    def __set_om_names(self,
                       si_df: pd.DataFrame  # DF of selected Si with OM list
                       ) -> pd.DataFrame:
        """set the list of OM names to the all Si"""
        df: pd.DataFrame = si_df.copy()
        df['OM_name']: list[list[str]] = [None for _ in range(len(si_df))]
        for item, _ in si_df.iterrows():
            df.at[item, 'OM_name'] = self.replace_oxy_name[item]
        return df

    def __get_si_om(self) -> list[int]:
        """return list of the Si bonded to the OM"""
        return list(self.si_df['atom_id'])

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
                 si_df: pd.DataFrame,  # All selected Si atoms
                 o_group: list[str]  # Name groups[O] to delete
                 ) -> None:
        self.o_delete: list[int]  # All the OD atoms to delete
        self.o_delete, self.bonded_si, self.si_df = self.__get_oxgygen(silica,
                                                                       Si_OM,
                                                                       o_group,
                                                                       si_df)

    def __get_oxgygen(self,
                      silica: rdlmp.ReadData,  # Atoms df in lammps full atom
                      Si_OM:  list[int],  # With selected group[Si] & OM bonded
                      o_group:  list[typing.Any],  # Name|index groups[O]
                      si_df: pd.DataFrame  # All selected Si atoms
                      ) -> tuple[list[int], list[int], pd.DataFrame]:
        """Find the hydrogen which have bonds with the selected Silicons"""
        o_list: list[pd.DataFrame] = []  # df of all Oxygen groups
        Atoms = silica.Atoms_df.copy()
        for item in o_group:
            o_list.append(Atoms[Atoms['name'] == item])
        df_o: pd.DataFrame = pd.concat(o_list)  # All the O atoms with names
        o_delete, bonded_si, si_df = self.__get_o_delete(silica.Bonds_df,
                                                         Si_OM,
                                                         df_o,
                                                         si_df)
        si_df = self.__set_ox_name(df_o, si_df, o_delete)
        pick = pickSi.PickSi(si_df, silica.diameter)
        si_df = pick.si_df
        o_delete = list(si_df['Ox_drop'])
        total_charge: float = self.__get_charges(df_o, o_delete)
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\t-> There are "{len(si_df)}" `Si` atoms bonded to the '
              f'`Ox` atoms\n'
              f'\t"{len(o_delete)}" `O` atoms are selcted to delete'
              f' with total charge: "{total_charge: .4f}"\n'
              f'\tThis gives the total coverage of '
              f'"{len(si_df)/(np.pi*silica.diameter*silica.diameter):.4f}"'
              f'{bcolors.ENDC}')
        return o_delete, bonded_si, si_df

    def __set_ox_name(self,
                      df_o: pd.DataFrame,  # Ox atoms df in lammps full atom
                      si_df: pd.DataFrame,  # Updated df Si with all
                      o_delete: list[int]  # Index of selected Ox atoms to drop
                      ) -> pd.DataFrame:
        """set the index and name of the selcted Ox to drop for each
        Si"""
        si_df['Ox_drop']: list[int] = o_delete
        df: pd.DataFrame = si_df.copy()
        df['Ox_drop_name']: list[str] = [None for _ in o_delete]
        for item, row in si_df.iterrows():
            ind: int = row['Ox_drop']  # atom_id of the Ox
            name: str = df_o[df_o['atom_id'] == ind]['name'][ind]
            df.at[item, 'Ox_drop_name'] = name
        return df

    def __get_o_delete(self,
                       bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                       Si_OM: list[int],  # With selected group[Si] & OM bonded
                       df_o: pd.DataFrame,  # DF with selected Oxygen
                       si_df: pd.DataFrame  # All selected Si atoms
                       ) -> tuple[list[int], list[int], pd.DataFrame]:
        # Return the index of the of Ox to delete, Si to attach chain and df
        # get bonded Ox atoms
        check_dict: dict[int, list[int]]  # Si: [ox bondec]
        check_dict = self.__get_check_dict(bonds_df, Si_OM, df_o)
        si_df = self.__drop_si(si_df, check_dict)
        # Drop Si without bonds to any Ox
        # Add a column with all Ox bonded to them
        si_df = self.__mk_ox_column(si_df, check_dict)
        bonded_si: list[int] = []  # Si bonded to the oxygens
        bonded_o: list[int] = []  # All O bonded to the Silica
        bonded_selected_o: list[int] = []  # All O bonded to the Silica, select
        for k, v in check_dict.items():
            bonded_si.append(k)
            bonded_o.extend(v)
            if len(v) == 1:
                bonded_selected_o.append(v[0])
            elif len(v) > 1:
                bonded_selected_o.append(self.__get_o_drop(v, df_o))
        return bonded_selected_o, bonded_si, si_df

    def __drop_si(self,
                  si_df: pd.DataFrame,  # All the selected Si
                  check_dict: dict[int, list[int]]  # Si: Ox
                  ) -> pd.DataFrame:
        """drop Si that does not have any bonds with Ox"""
        df: pd.DataFrame = si_df.copy()
        si_list: list[int]  # atom_id of the selected Si
        si_list = list(check_dict.keys())
        for item, row in si_df.iterrows():
            if row['atom_id'] not in si_list:
                df.drop(index=item, axis=0, inplace=True)
        return df

    def __get_o_drop(self,
                     o_list: list[int],  # OM atoms bond to the Si
                     df_o: pd.DataFrame  # All OM atoms
                     ) -> int:
        """return the index of the Ox to drop for Si atoms with more
        than one Ox bond to them. If OD is among them, it returns that
        one; otherwise, it returns the atom with a smaller index."""
        o_to_drop = None  # The atom_id of the Ox atom to drop
        o_names: dict[int, str]  # atom_id: Name of the ox atoms in the o_list
        o_names = {item: df_o[df_o['atom_id'] == item]['name'][item]
                   for item in o_list}
        for k, v in o_names.items():
            if v == 'OD':
                o_to_drop = k
        if o_to_drop is None:
            o_to_drop = np.min(o_list)
        return o_to_drop

    def __get_check_dict(self,
                         bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                         Si_OM: list[int],  # With selected group[Si] & OM bond
                         df_o: pd.DataFrame,  # DF with selected Oxygen
                         ) -> dict[int, list[int]]:
        """a dict with Si and bonded Ox to them """
        all_o = list(df_o['atom_id'])
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

    def __mk_ox_column(self,
                       si_df: pd.DataFrame,  # All selected Si
                       check_dict: dict[int, list[int]]  # Si: Ox bonded
                       ) -> pd.DataFrame:
        """return the si_df wiht Ox info column"""
        si_df['Ox_list']: list[typing.Any]  # Ox atoms bonded to the Si
        si_df['Ox_list'] = [None for _ in range(len(si_df))]
        i_i: int = 0  # Count the number of Ox bonded to the Si
        i_ii: int = 0  # Count the number of Ox bonded to the Si
        for item, _ in si_df.iterrows():
            si_df.at[item, 'Ox_list'] = check_dict[item]
            if len(check_dict[item]) > 1:
                if len(check_dict[item]) == 2:
                    i_i += 1
                elif len(check_dict[item]) == 3:
                    i_ii += 1
        print(f'{bcolors.WARNING}{self.__class__.__name__}: \n'
              f'\tThere are "{i_i}" `Si` atoms bonded to '
              f'two, and "{i_ii}" `Si` atoms bonded to three Ox groups'
              f'{bcolors.ENDC}')
        return si_df

    def __get_charges(self,
                      df_o: pd.DataFrame,  # All Ox infos
                      delete_list: list[int]  # Ox atom_id to drop
                      ) -> float:
        """return the total charges of the deleted Oxygens"""
        total_charge: float = 0  # Charge of all deleted Oxygens
        for _, row in df_o.iterrows():
            if row['atom_id'] in delete_list:
                total_charge += row['charge']
        return total_charge


class GetHyGroups:
    """Find Hydrogen groups bonded to the Oxygen atoms"""
    def __init__(self,
                 silica: rdlmp.ReadData,  # All silica atoms
                 o_delete: list[int],  # Index of O atoms to delete
                 h_group: list[str]  # name of the H atoms to check for bonds
                 ) -> None:
        self.h_delete: list[int]  # List of all the H atoms to delete
        self.h_delete = self.__do_hydrogens(silica, o_delete, h_group)

    def __do_hydrogens(self,
                       silica: rdlmp.ReadData,  # All silica atoms
                       o_delete: list[int],  # Index of O atoms to delete
                       h_group: list[str]  # name of the H atoms to check
                       ) -> list[int]:
        df_h: pd.DataFrame   # All the hydrogens with the selcted type
        df_h = self.__get_hydrogens(silica.Atoms_df, h_group)
        h_delete: list[int] = self.__h_delete(silica.Bonds_df, o_delete, df_h)
        return h_delete

    def __get_hydrogens(self,
                        silica_atoms: pd.DataFrame,  # All silica atoms
                        h_group: list[str]  # name of the H atoms to check for
                        ) -> pd.DataFrame:
        # Get Hydrogen atoms
        h_list: list[pd.DataFrame] = []  # df of all Hydrogens
        for item in h_group:
            h_list.append(silica_atoms[silica_atoms['name'] == item])
        return pd.concat(h_list)

    def __h_delete(self,
                   bonds_df: pd.DataFrame,  # All the bonds in silica
                   o_delete: list[int],  # Index of the O atoms
                   df_h: pd.DataFrame  # Df of all the Hydrogens
                   ) -> list[int]:
        all_h = list(df_h['atom_id'])
        delete_list: list[int] = []  # index of H atoms to delete
        for _, row in bonds_df.iterrows():
            if row['ai'] in o_delete:
                if row['aj'] in all_h:
                    delete_list.append(row['aj'])
            if row['aj'] in o_delete:
                if row['ai'] in all_h:
                    delete_list.append(row['ai'])
        total_charges: float = self.__get_charges(df_h, delete_list)
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\t"{len(delete_list)}" `H` atoms bonded to the '
              f'slected O to delete, with total charge of: '
              f'"{total_charges:.4f}"'
              f'{bcolors.ENDC}')
        return delete_list

    def __get_charges(self,
                      df_h: pd.DataFrame,  # All H infos
                      delete_list: list[int]  # H atom_id to drop
                      ) -> float:
        """return the total charges of the deleted Hydrogens"""
        total_charge: float = 0  # Charge of all deleted Hydrogens
        for _, row in df_h.iterrows():
            if row['atom_id'] in delete_list:
                total_charge += row['charge']
        return total_charge
