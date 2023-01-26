import sys
import typing
import numpy as np
import pandas as pd
import read_lmp_data as rdlmp
import write_lmp as wrlmp
from colors_text import TextColor as bcolors


class Doc:
    """read data files and update positions and index of them"""


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
        df: pd.DataFrame = pd.concat(Si_list)  # Df of all the Si to replace
        # Drop some columns
        df.drop(axis=1, columns=['nx', 'ny', 'nz', 'cmt', 'b_name'],
                inplace=True)
        print(f'{bcolors.OKBLUE}{self.__class__.__name__}:\n'
              f'\tThere are: {len(df)} atoms with selected '
              f'atoms name [Si] in the file{bcolors.ENDC}')
        max_radius: float = self.__get_max_radius(Atoms)
        df = self.__get_angles(df)
        df = df[(df['rho'] >= max_radius - 5)]
        print(f'{bcolors.OKBLUE}\tThere are: {len(df)} Si atoms in the '
              f'choosen area of the system, Max_radius = {max_radius:.5f}'
              f'{bcolors.ENDC}')
        return df

    def __get_max_radius(self,
                         Atoms: pd.DataFrame,  # Atoms in lammps full atom
                         ) -> float:
        """return the radius of the nano-particles"""
        x_max: float = Atoms['x'].abs().max()
        y_max: float = Atoms['y'].abs().max()
        z_max: float = Atoms['z'].abs().max()
        return np.max([x_max, y_max, z_max])

    def __get_angles(self,
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


class GetOGroups:
    """get Oxygen and/or Hydrogen groups to delete them and update
    data file.
    Oxygen is bonded to the silica and Hydrogen if asked for is bonded
    to the Oxygen.
    """
    def __init__(self,
                 silica: rdlmp.ReadData,  # Atoms df in form of lammps fullatom
                 Si_delete: list[int],  # With selected group[Si]
                 Ogroup:  list[typing.Any],  # Name | index groups[O] to delete
                 fraction: float = 1  # Fraction of to select from, 0<fr<=1
                 ) -> None:
        self.O_delete: list[int]  # All the O atoms to delete
        self.O_delete = self.__get_oxgygen(silica, Si_delete, Ogroup)

    def __get_oxgygen(self,
                      silica: rdlmp.ReadData,  # Atoms df in lammps full atom
                      Si_delete: list[int],  # With selected group[Si]
                      Ogroup:  list[typing.Any],  # Name|index groups[O]
                      ) -> list[int]:
        """Find the hydrogen which have bonds with the selected Silicons"""
        O_list: list[pd.DataFrame] = []  # df of all Oxygen groups
        Atoms = silica.Atoms_df
        for item in Ogroup:
            O_list.append(Atoms[Atoms['name'] == item])
        df: pd.DataFrame = pd.concat(O_list)
        O_delete = self.__get_o_delete(silica.Bonds_df, Si_delete, df)
        return O_delete

    def __get_o_delete(self,
                       bonds_df: pd.DataFrame,  # Bonds in the LAMMPS format
                       Si_delete: list[int],  # With selected group[Si]
                       df_o: pd.DataFrame  # DF with selected Oxygen
                       ) -> list[int]:  # index of the O to delete
        # get bonds
        all_o = [item for item in df_o['atom_id']]
        delete_list: list[int] = []  # index of the O atoms to delete
        for _, row in bonds_df.iterrows():
            if row['ai'] in Si_delete or row['aj'] in Si_delete:
                if row['ai'] in all_o and row['ai'] not in delete_list:
                    delete_list.append(row['ai'])
                if row['aj'] in all_o and row['aj'] not in delete_list:
                    delete_list.append(row['aj'])
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}:\n'
              f'\tThere are {len(delete_list)} O atoms bonded to the '
              f'slected Si{bcolors.ENDC}')
        return delete_list


class Delete:
    """delete selcted atoms"""
    def __init__(self,
                 silica: rdlmp.ReadData  # Data from LAMMPS
                 ) -> None:
        silicons = GetSiGroups(silica.Atoms_df,
                               Sigroup=['SD'],
                               fraction=1)
        oxygens = GetOGroups(silica,
                             silicons.Si_delete,
                             Ogroup=['OD'],
                             fraction=1)
        self.Si_df = silicons.df_Si
        self.__delete_all(silica, oxygens.O_delete)

    def __delete_all(self,
                     silica: rdlmp.ReadData,  # Data from LAMMPS
                     O_delete: list[int],  # Index of O atoms to delete
                     ) -> None:
        old_new_dict: dict[int, int]  # new and old index of updated atoms df
        self.UAtoms_df: pd.DataFrame  # Atoms  with updated atoms' index
        self.UVelocities: pd.DataFrame  # Velocities with updated atoms' index
        self.UBonds_df: pd.DataFrame  # Bonds with updated atoms' index
        self.UAngles_df: pd.DataFrame  # Angles with updated atoms' index
        delete_group: list[int] = []  # To extend all selected atoms
        delete_group.extend(O_delete)
        old_new_dict, self.UAtoms_df = self.__update_atoms(silica,
                                                           delete_group)
        self.UVelocities = self.__update_velocities(silica.Velocities_df,
                                                    old_new_dict,
                                                    delete_group)
        self.UBonds_df = self.__update_bonds(silica.Bonds_df,
                                             old_new_dict,
                                             delete_group)
        self.UAngles_df = self.__update_angles(silica.Angles_df,
                                               old_new_dict,
                                               delete_group)
        self.USi_df = self.__update_selected_Si(old_new_dict)

    def __update_selected_Si(self,
                             old_new_dict: dict[int, int]  # old:new atom id
                             ) -> pd.DataFrame:
        """update atom index of the deleted atoms indes"""
        df: pd.DataFrame = self.Si_df.copy()
        new_ai = []  # New index for ai
        for item, row in df.iterrows():
            new_ai.append(old_new_dict[item])
        df['atom_id'] = new_ai
        return df

    def __update_atoms(self,
                       silica: rdlmp.ReadData,  # Atoms df in lammps full atom
                       delete_group: list[int]  # Index of the atom to delete
                       ) -> tuple[dict[int, int], pd.DataFrame]:
        """delete atoms and return update version in LAMMPS format"""
        Atoms_df: pd.DataFrame  # DF with removed atoms
        Atoms_df = self.__delete_atoms(silica.Atoms_df, delete_group)
        old_new_dict: dict[int, int]  # Dict with old and new atom id
        old_new_dict = dict(zip(Atoms_df['old_atom_id'], Atoms_df['atom_id']))
        return old_new_dict, Atoms_df

    def __delete_atoms(self,
                       Atoms_df: pd.DataFrame,  # Atoms in LAMMPS format
                       delete_group: list[int]  # Index of atom to delete
                       ) -> pd.DataFrame:
        """delete atoms and return updated one with a new column"""
        Atoms_df['old_atom_id'] = Atoms_df['atom_id']
        df = Atoms_df.copy()
        for item, row in df.iterrows():
            if row['atom_id'] in delete_group:
                Atoms_df.drop(index=[item], axis=0, inplace=True)
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}:\n'
              f'\t {len(delete_group)} atoms is deleted from data file'
              f'{bcolors.ENDC}')
        del df
        Atoms_df.reset_index(inplace=True)
        Atoms_df.index += 1
        Atoms_df['atom_id'] = Atoms_df.index
        Atoms_df.drop(columns=['index'], inplace=True)
        return Atoms_df

    def __update_velocities(self,
                            Velocities_df: pd.DataFrame,  # In LAMMPS format
                            old_new_dict: dict[int, int],  # old:new atom id
                            delete_group: list[int]  # Index of atom to delete
                            ) -> pd.DataFrame:
        """delete velocities for deleted atoms"""
        df = Velocities_df.copy()
        del_counter: int = 0  # count the numbers of deleted velocities
        for item, row in df.iterrows():
            if item in delete_group:
                Velocities_df.drop(index=[item], axis=0, inplace=True)
                del_counter += 1
        print(f'{bcolors.OKBLUE}\t {del_counter} velocities are deleted '
              f'from the data file{bcolors.ENDC}')
        new_ai = []  # New index for ai
        for item, row in Velocities_df.iterrows():
            new_ai.append(old_new_dict[item])
        del df
        Velocities_df.index = new_ai
        return Velocities_df

    def __update_bonds(self,
                       Bonds_df: pd.DataFrame,  # Atoms in LAMMPS format
                       old_new_dict: dict[int, int],  # Dict old: new atom id
                       delete_group: list[int]  # Index of atom to delete
                       ) -> pd.DataFrame:
        """delete bonds for deleted atoms"""
        df = Bonds_df.copy()
        del_counter: int = 0  # count the numbers of deleted bonds
        for item, row in df.iterrows():
            if row['ai'] in delete_group or row['aj'] in delete_group:
                Bonds_df.drop(index=[item], axis=0, inplace=True)
                del_counter += 1
        print(f'{bcolors.OKBLUE}\t {del_counter} bonds are deleted '
              f'from the data file{bcolors.ENDC}')
        new_ai = []  # New index for ai
        new_aj = []  # New index for aj
        for item, row in Bonds_df.iterrows():
            ai = row['ai']
            aj = row['aj']
            new_ai.append(old_new_dict[ai])
            new_aj.append(old_new_dict[aj])
        del df
        Bonds_df['ai'] = new_ai
        Bonds_df['aj'] = new_aj
        return Bonds_df

    def __update_angles(self,
                        Angels_df: pd.DataFrame,  # Atoms in LAMMPS format
                        old_new_dict: dict[int, int],  # Dict old: new atom id
                        delete_group: list[int]  # Index of atom to delete
                        ) -> pd.DataFrame:
        """delete angles for deleted atoms"""
        df = Angels_df.copy()
        del_counter: int = 0  # count the numbers of deleted angles
        for item, row in df.iterrows():
            if row['ai'] in delete_group\
               or row['aj'] in delete_group\
               or row['ak'] in delete_group:
                Angels_df.drop(index=[item], axis=0, inplace=True)
                del_counter += 1
        print(f'{bcolors.OKBLUE}\t {del_counter} angles are deleted '
              f'from the data file{bcolors.ENDC}')
        new_ai = []  # New index for ai
        new_aj = []  # New index for aj
        new_ak = []  # New index for ak
        for item, row in Angels_df.iterrows():
            ai = row['ai']
            aj = row['aj']
            ak = row['ak']
            new_ai.append(old_new_dict[ai])
            new_aj.append(old_new_dict[aj])
            new_ak.append(old_new_dict[ak])
        del df
        Angels_df['ai'] = new_ai
        Angels_df['aj'] = new_aj
        Angels_df['ak'] = new_ak
        return Angels_df


class UpdateCoords:
    """update all the attributes to dataframe to the updated data"""
    def __init__(self,
                 fname: str,  # Main data
                 ) -> None:
        silica = GetData(fname)
        update = Delete(silica)
        self.__set_attrs(silica, update)

    def __set_attrs(self,
                    silica: GetData,  # Tha main data
                    update: Delete  # All the data read from file
                    ) -> None:
        """update all the attrs"""
        self.Atoms_df = update.UAtoms_df
        self.Velocities_df = update.UVelocities
        self.Bonds_df = update.UBonds_df
        self.Angles_df = update.UAngles_df
        self.Masses_df = silica.Masses_df
        self.Si_df = update.USi_df
        self.NAtoms = len(update.UAtoms_df)
        self.Nmols = np.max(update.UAtoms_df['mol'])
        self.Dihedrals_df = self.__set_dihedrlas()

    def __set_dihedrlas(self) -> pd.DataFrame:
        """make empty df"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'ak', 'ah']
        df = pd.DataFrame(columns=columns)
        return df


if __name__ == '__main__':
    fname = sys.argv[1]
    update = UpdateCoords(fname)
    wrt = wrlmp.WriteLmp(obj=update, output='after_del.data')
    wrt.write_lmp()
