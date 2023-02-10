import sys
import numpy as np
import pandas as pd
import write_lmp as wrlmp
import get_atoms as gtatom
import static_info as stinfo
import read_lmp_data as rdlmp
import update_charges as upcharge
from colors_text import TextColor as bcolors


class Doc:
    """read data files and update positions and index of them"""


class Delete:
    """delete selcted atoms"""
    def __init__(self,
                 silica: rdlmp.ReadData  # Data from LAMMPS
                 ) -> None:
        """find silicon atoms on the shell, then find the Oxygen which
        is attached to the silicons and delete them; also, it should
        check if there is any hydrogen attached to the selected O and
        drop them too"""
        # Find Si on the shell
        silicons = gtatom.GetSiGroups(silica.Atoms_df,
                                      Sigroup=stinfo.AtomGroup.SiGroup,
                                      fraction=1)
        # Get OM atoms bonded to the selected Si, and drop the Si in the Body
        om_groups = gtatom.GetOmGroups(silica,
                                       silicons.df_Si,
                                       OMgroup=stinfo.AtomGroup.OMGroup
                                       )
        # Get Ox atoms, which should drop and replace
        oxygens = gtatom.GetOxGroups(silica,
                                     om_groups.Si_OM,
                                     om_groups.Si_df,
                                     Ogroup=stinfo.AtomGroup.OxGroup,
                                     fraction=1)
        # Get hydrogen bonded to the selected oxygen, to drop
        hydrogens = gtatom.GetHyGroups(silica,
                                       oxygens.O_delete,
                                       Hgroup=stinfo.AtomGroup.HyGroup)
        # Drop selected O atached to the Si and if there is H atom bond to them
        # Get the O which bonded to the selected Si, to make angles and torsion
        self.Si_df: pd.DataFrame = oxygens.Si_df
        self.old_new_dict: dict[int, int] = \
            self.__delete_all(silica, oxygens, hydrogens.H_delete, om_groups)

    def __delete_all(self,
                     silica: rdlmp.ReadData,  # Data from LAMMPS
                     oxygens: gtatom.GetOxGroups,  # Index of O atoms to delete
                     H_delete: list[int],  # Index of H atoms to delete
                     om_groups: gtatom.GetOmGroups  # O atoms to replace
                     ) -> dict[int, int]:
        old_new_dict: dict[int, int]  # new and old index of updated atoms df
        self.UAtoms_df: pd.DataFrame  # Atoms with updated atoms' index
        self.UVelocities: pd.DataFrame  # Velocities with updated atoms' index
        self.UBonds_df: pd.DataFrame  # Bonds with updated atoms' index
        self.UAngles_df: pd.DataFrame  # Angles with updated atoms' index
        delete_group: list[int] = []  # To extend all selected atoms
        delete_group.extend(oxygens.O_delete)
        delete_group.extend(H_delete)
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
        USi_df = self.__update_selected_Si(old_new_dict)
        self.USi_df = self.__append_om(USi_df,
                                       old_new_dict,
                                       om_groups.replace_oxy)
        return old_new_dict

    def __update_selected_Si(self,
                             old_new_dict: dict[int, int]  # old:new atom id
                             ) -> pd.DataFrame:
        """update atom index of the deleted atoms indes"""
        df: pd.DataFrame = self.Si_df.copy()
        new_ai = []  # New index for ai
        for item, _ in df.iterrows():
            new_ai.append(old_new_dict[item])
        df['atom_id'] = new_ai
        return df

    def __update_atoms(self,
                       silica: rdlmp.ReadData,  # Atoms df in lammps full atom
                       delete_group: list[int],  # Index of the atom to delete
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
        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\t{len(delete_group)} atoms is deleted from data file'
              f'{bcolors.ENDC}')
        del df
        Atoms_df.reset_index(inplace=True)
        Atoms_df.index += 1
        Atoms_df['atom_id'] = Atoms_df.index
        Atoms_df.drop(columns=['index'], inplace=True)
        return Atoms_df

    def __append_om(self,
                    Si_df: pd.DataFrame,  # Updated Si df
                    old_new_dict: dict[int, int],  # old: new atom id
                    replace_oxy: dict[int, list[int]]  # Si atoms with O atoms
                    ) -> pd.DataFrame:
        """append indices of OM atoms which are bonded to the Si"""
        om: dict[int, list[int]] = self.__update_om_id(old_new_dict,
                                                       replace_oxy)
        df: pd.DataFrame = Si_df.copy()
        df['OM_list'] = [None for _ in df.index]
        for item, row in Si_df.iterrows():
            df.at[item, 'OM_list'] = om[row['atom_id']]
        return df

    def __update_om_id(self,
                       old_new_dict: dict[int, int],  # old: new atom id
                       replace_oxy: dict[int, list[int]]  # Si atom with O atom
                       ) -> dict[int, list[int]]:
        """update the atom_id of the OM atoms bonded to the Si"""
        om: dict[int, list[int]] = {}
        for k, v in replace_oxy.items():
            new_k: int = old_new_dict[k]
            om[new_k] = []
            for item in v:
                om[new_k].append(old_new_dict[item])
        return om

    def __update_velocities(self,
                            Velocities_df: pd.DataFrame,  # In LAMMPS format
                            old_new_dict: dict[int, int],  # old:new atom id
                            delete_group: list[int]  # Index of atom to delete
                            ) -> pd.DataFrame:
        """delete velocities for deleted atoms"""
        df = Velocities_df.copy()
        del_counter: int = 0  # count the numbers of deleted velocities
        for item, _ in df.iterrows():
            if item in delete_group:
                Velocities_df.drop(index=[item], axis=0, inplace=True)
                del_counter += 1
        print(f'{bcolors.OKBLUE}\t{del_counter} velocities are deleted '
              f'from the data file{bcolors.ENDC}')
        new_ai = []  # New index for ai
        for item, _ in Velocities_df.iterrows():
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
        print(f'{bcolors.OKBLUE}\t{del_counter} bonds are deleted '
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
        print(f'{bcolors.OKBLUE}\t{del_counter} angles are deleted '
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
        silica = gtatom.GetData(fname)
        update = Delete(silica)
        upq = upcharge.UpdateCharge(update.UAtoms_df,
                                    update.Si_df,
                                    update.old_new_dict)
        self.__set_attrs(silica, update, upq)
        self.__write_infos()

    def __set_attrs(self,
                    silica: gtatom.GetData,  # Tha main data
                    update: Delete,  # All the data read from file
                    upq: upcharge.UpdateCharge  # Atoms with updated charges
                    ) -> None:
        """update all the attrs"""
        self.Atoms_df = upq.Atoms_df
        self.Velocities_df = update.UVelocities
        self.Bonds_df = update.UBonds_df
        self.Angles_df = update.UAngles_df
        self.Masses_df = silica.Masses_df
        self.Si_df = update.USi_df
        self.NAtoms = len(update.UAtoms_df)
        self.Nmols = np.max(update.UAtoms_df['mol'])
        self.Dihedrals_df = self.__set_dihedrlas()
        self.NAtomTyp: int = np.max(self.Masses_df['typ'])
        self.NBonds: int = len(self.Bonds_df)
        self.NBondTyp: int = np.max(self.Bonds_df['typ'])
        self.NAngles: int = len(self.Angles_df)
        self.NAngleTyp: int = np.max(self.Angles_df['typ'])
        self.NDihedrals: int = len(self.Dihedrals_df)
        self.NDihedralTyp: int = np.max(self.Dihedrals_df['typ'])

    def __set_dihedrlas(self) -> pd.DataFrame:
        """make empty df"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'ak', 'ah']
        df = pd.DataFrame(columns=columns)
        return df

    def __write_infos(self) -> None:
        print(f'{bcolors.OKGREEN}\tData Summary after deleting atoms'
              f' and updataing charges:\n'
              f'\t\t# Atoms: {self.NAtoms}, # Atom`s types: {self.NAtomTyp}\n'
              f'\t\t# Bonds: {self.NBonds}, # Bond`s types: {self.NBondTyp}\n'
              f'\t\t# Angles: {self.NAngles}, '
              f'# Angle`s types: {self.NAngleTyp}\n'
              f'\t\t# Dihedrals: {self.NDihedrals}, '
              f'# Dihedral`s types: {self.NDihedralTyp}\n'
              f'\t\tTotal charge: {self.Atoms_df["charge"].sum():.4f}\n'
              f'\t\tMin charge: {self.Atoms_df["charge"].min()}\n'
              f'\t\tMax charge: {self.Atoms_df["charge"].max()}'
              f'{bcolors.ENDC}'
              )


if __name__ == '__main__':
    fname = sys.argv[1]
    update = UpdateCoords(fname)
    wrt = wrlmp.WriteLmp(obj=update, output='after_del.data')
    wrt.write_lmp()
