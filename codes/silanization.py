import re
import sys
import numpy as np
import pandas as pd
import read_lmp_data as rdlmp
import update_coords as upcord
from colors_text import TextColor as bcolors
import write_lmp as wrlmp


class Doc:
    """read data
        select the gropus to replace
        write the output file
    """


class GetAmino(rdlmp.ReadData):
    """read the main Aminopropyle coordinates and put the Si position
    to zero"""
    def __init__(self) -> None:
        fname: str = '/scratch/saeed/MyScripts/np_silica/data/aminopropyl.data'
        super().__init__(fname)
        self.__set_attr()

    def __set_attr(self) -> None:
        """set some attributes to data file"""
        pass


class PrepareAmino:
    """get the data for selected Silicons:
        first rotate aminopropyle for each of one the Si group
        replace the Si of each aminopropyle with the with the index of
        the selected group and updated Si
        """
    def __init__(self,
                 update: upcord.UpdateCoords,  # All the information of silica
                 amino: GetAmino  # Information about one aminos
                 ) -> None:
        """apply the position update to the aminos"""
        self.__update_aminos(update, amino)
        self.__write_infos()

    def __update_aminos(self,
                        update: upcord.UpdateCoords,  # Information of silica
                        amino: GetAmino  # Information about one aminos
                        ) -> None:
        """do"""
        si_df: pd.DataFrame  # Si groups with rotation angles
        si_df = update.Si_df
        si_df = self.__order_si_df(si_df)
        self.Si = 'Si'
        Atoms_df = self.__to_origin(amino.Atoms_df)
        Atoms_df = self.__get_azimuths(Atoms_df)
        Atoms_list: list[pd.DataFrame] = []  # Keep all the aminos to concate
        Bonds_list: list[pd.DataFrame] = []  # Keep all the aminos to concate
        Angles_list: list[pd.DataFrame] = []  # Keep all the aminos to concate
        Dihedrals_list: list[pd.DataFrame] = []  # Keep all the aminos> concate
        amino_masses: pd.DataFrame  # To update the atom index in masses
        amino_masses = self.__update_mass_index(update.Masses_df,
                                                amino.Masses_df)
        Atoms_df = self.__update_atom_types(Atoms_df, amino_masses)
        amino.Bonds_df = self.__update_boandi_types(update.Bonds_df,
                                                    amino.Bonds_df)
        amino.Angles_df = self.__update_boandi_types(update.Angles_df,
                                                     amino.Angles_df)
        amino.Dihedrals_df = self.__update_boandi_types(update.Dihedrals_df,
                                                        amino.Dihedrals_df)
        for item, row in si_df.iterrows():
            # Si from amino will be deleted later, so the rest of
            # atoms must start on lower
            atom_level: int = (item - 1) * \
                                (amino.NAtoms - 1) + update.NAtoms - 1
            mol_level: int = item + update.Nmols
            i_amino = self.__upgrade_amino(Atoms_df.copy(),
                                           row,
                                           atom_level,
                                           mol_level)
            amino.Atoms_df = i_amino
            Atoms_list.append(self.__drop_si(i_amino))
            boandi = UpdateBoAnDi(amino)  # Update bonds, angles, dihedrals
            Bonds_list.append(boandi.Bonds_df)
            Angles_list.append(boandi.Angles_df)
            Dihedrals_list.append(boandi.Dihedrals_df)
        self.All_amino_atoms = pd.concat(Atoms_list, ignore_index=True)
        self.All_amino_atoms.index += 1
        self.All_amino_bonds = pd.concat(Bonds_list, ignore_index=True)
        self.All_amino_bonds.index += 1
        self.All_amino_angles = pd.concat(Angles_list, ignore_index=True)
        self.All_amino_angles.index += 1
        self.All_amino_dihedrals = pd.concat(Dihedrals_list, ignore_index=True)
        self.All_amino_dihedrals.index += 1
        self.Masses_df = self.__drop_si_mass(amino_masses)

        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\t{len(self.All_amino_atoms)} atoms'
              f', {len(self.All_amino_bonds)} bonds'
              f', {len(self.All_amino_angles)} angles'
              f', {len(self.All_amino_dihedrals)} dihedrals'
              f' is updated{bcolors.ENDC}')

    def __upgrade_amino(self,
                        amino_atoms: pd.DataFrame,  # Amino Atoms information,
                        si_row: pd.DataFrame,  # A row of the dataframe
                        atom_level: int,  # To increase the atoms index not Si
                        mol_level: int  # To increase the mol index
                        ) -> pd.DataFrame:
        """rotate each amino based on the azimuuth polar angle in the row"""
        amino_atoms['old_id'] = amino_atoms['atom_id']
        for item, _ in amino_atoms.iterrows():
            if item != 1:
                amino_atoms.at[item, 'mol'] += mol_level
                amino_atoms.at[item, 'atom_id'] += atom_level
            else:
                amino_atoms.at[item, 'atom_id'] = si_row['atom_id']
                amino_atoms.at[item, 'mol'] = si_row['mol']
                amino_atoms.at[item, 'typ'] = si_row['typ']
                amino_atoms.at[item, 'x'] = si_row['x']
                amino_atoms.at[item, 'y'] = si_row['y']
                amino_atoms.at[item, 'z'] = si_row['z']
                amino_atoms.at[item, 'rho'] = si_row['rho']
                amino_atoms.at[item, 'azimuth'] = si_row['azimuth']
                amino_atoms.at[item, 'polar'] = si_row['polar']
        return self.__rotate_amino(amino_atoms)

    def __drop_si(self,
                  amino_atoms: pd.DataFrame  # Rotated Atoms_df
                  ) -> pd.DataFrame:
        """drop the silicon, since it is already in the main data"""
        df: pd.DataFrame = amino_atoms.copy()
        df.drop(amino_atoms[amino_atoms['name'] == 'Si'].index, inplace=True)
        df.reset_index(inplace=True)
        df.drop(columns=['index', 'old_id'], inplace=True, axis=1)
        df.index += 1
        return df

    def __drop_si_mass(self,
                       amino_masses: pd.DataFrame  # Amino masses
                       ) -> pd.DataFrame:
        """Drop Si from amino masses, it is only one Si in atoms"""
        df: pd.DataFrame = amino_masses.copy()
        df.drop(amino_masses[amino_masses['name'] == 'Si'].index, inplace=True)
        df.reset_index(inplace=True)
        df.drop(columns=['index'], inplace=True, axis=1)
        df.index += 1
        return df

    def __update_mass_index(self,
                            silica_masses: pd.DataFrame,  # DF of masses
                            amino_masses: pd.DataFrame  # DF of masses
                            ) -> pd.DataFrame:
        """update atoms masses index in Masses for Amino"""
        silica_ind: int = np.max(silica_masses['typ'])
        df: pd.DataFrame = amino_masses.copy()
        df['old_typ'] = df['typ']
        for item, _ in amino_masses.iterrows():
            df.at[item, 'typ'] += silica_ind
        return df

    def __update_atom_types(self,
                            Atoms_df: pd.DataFrame,  # Atoms of amino
                            amino_masses: pd.DataFrame  # Updated masses
                            ) -> pd.DataFrame:
        """update type of atoms before rotation and indexing"""
        old_new_dict: dict[int, int]  # To change the types
        old_new_dict = {k: v for k, v in zip(amino_masses['old_typ'],
                                             amino_masses['typ'])}
        df: pd.DataFrame = Atoms_df.copy()
        for item, row in Atoms_df.iterrows():
            df.at[item, 'typ'] = old_new_dict[row['typ']]
        return df

    def __update_boandi_types(self,
                              silica: pd.DataFrame,  # Silica bonds/angles/dihe
                              amino: pd.DataFrame  # Amino bonds/angles/dihedra
                              ) -> pd.DataFrame:
        """update bonds type in aminos"""
        silica_ind: int = np.max(silica['typ'])
        df: pd.DataFrame = amino.copy()
        if not np.isnan(silica_ind):
            for item, _ in amino.iterrows():
                df.at[item, 'typ'] += silica_ind
        return df

    def __rotate_amino(self,
                       amino_atoms: pd.DataFrame  # Updated amino
                       ) -> pd.DataFrame:
        """rotate amino based on the azimuth and polar of the Si"""
        ref_x: float = amino_atoms['x'][1]
        ref_y: float = amino_atoms['y'][1]
        ref_z: float = amino_atoms['z'][1]
        new_x: list[float] = []  # After rotations
        new_y: list[float] = []  # After rotations
        new_z: list[float] = []  # After rotations
        beta: float = amino_atoms['azimuth'][1]  # Si where chain will attach
        gamma: float = amino_atoms['polar'][1]  # Si where chain will attach to
        new_x.append(amino_atoms['x'][1])
        new_y.append(amino_atoms['y'][1])
        new_z.append(amino_atoms['z'][1])
        for item, row in amino_atoms.iterrows():
            if item > 1:
                x_i, y_i, z_i = self.__rot_matrix(beta,
                                                  gamma,
                                                  row['x'],
                                                  row['y'],
                                                  row['z'])
                new_x.append(x_i + ref_x)
                new_y.append(y_i + ref_y)
                new_z.append(z_i + ref_z)
        amino_atoms['x'] = new_x
        amino_atoms['y'] = new_y
        amino_atoms['z'] = new_z
        return amino_atoms

    def __rot_matrix(self,
                     beta: float,  # Beta, the azimuth (pitch) angle
                     gamma: float,  # Gamma, the polar (roll) angle
                     x: float,  # X position
                     y: float,  # Y position
                     z: float  # Z position
                     ) -> tuple[float, float, float]:
        """calculate the rotated coordes"""
        h_pi: float = 0.5*np.pi
        if beta > 0:
            cos_g: float = np.cos(gamma - h_pi)
            sin_g: float = np.sin(gamma - h_pi)
        else:
            cos_g = np.cos(h_pi - gamma)
            sin_g = np.sin(h_pi - gamma)
        cos_b: float = np.cos(h_pi-beta)
        sin_b: float = np.sin(h_pi-beta)
        x_n = x*cos_b-y*sin_b
        y_n = x*sin_b+y*cos_b
        z_n = z
        x_new = x_n*cos_g+z_n*sin_g
        y_new = y_n
        z_new = -x_n*sin_g+z_n*cos_g
        # x_new = x_n*cos_g-z_n*sin_g
        # y_new = x_n*sin_g*sin_g+y_n*cos_g+z_n*sin_g*cos_g
        # z_new = x_n*cos_g*sin_g-y*sin_g+z_n*cos_g*cos_g
        # x_new: float = x*cos_g*cos_b - y*sin_g + z*sin_b*cos_g
        # y_new: float = x*sin_g*cos_b + y*cos_g + z*sin_b*sin_g
        # z_new: float = -x*sin_b + z*cos_b
        return x_new, y_new, z_new

    def __order_si_df(self,
                      si_df: pd.DataFrame  # Si groups with rotation angles
                      ) -> pd.DataFrame:
        # reindex the si_df to get orders
        si_df.reset_index(inplace=True)
        si_df.index += 1
        si_df.drop(columns=['index'], axis=1, inplace=True)
        return si_df

    def __get_azimuths(self,
                       amino_atoms: pd.DataFrame  # Atoms info of Aminopropyle
                       ) -> pd.DataFrame:
        """calculate the azimuth and polar angle of Amino"""
        df: pd.DataFrame = amino_atoms.copy()
        rho: list[float] = []  # Radius for all atoms
        azimuth: list[float] = []  # Azimuth for all atoms
        polar: list[float] = []  # Polar for all atoms
        for item, row in df.iterrows():
            if item == 1:
                rho.append(0)
                azimuth.append(0)
                polar.append(0)
            else:
                x = row['x']
                y = row['y']
                z = row['z']
                i_rho: float = np.sqrt(x*x + y*y + z*z)
                i_azimuth: float = np.arctan2(x, y)
                i_polar: float = np.arccos(z/i_rho)
                rho.append(i_rho)
                azimuth.append(i_azimuth)
                polar.append(i_polar)
        df = df.assign(rho=rho, azimuth=azimuth, polar=polar)
        return df

    def __to_origin(self,
                    amino_atoms: pd.DataFrame  # Df amino Atoms
                    ) -> pd.DataFrame:
        """put the coordinate of Si in Aminopropyle to zero"""
        df: pd.DataFrame = amino_atoms.copy()
        print(f'{bcolors.OKCYAN}\tMove Aminopropyle [Si] to origin'
              f'{bcolors.ENDC}')
        x_si: float = amino_atoms[amino_atoms['name'] == self.Si]['x'][1]
        y_si: float = amino_atoms[amino_atoms['name'] == self.Si]['y'][1]
        z_si: float = amino_atoms[amino_atoms['name'] == self.Si]['z'][1]
        df['x'] -= x_si
        df['y'] -= y_si
        df['z'] -= z_si
        return df

    def __write_infos(self) -> None:
        print(f'{bcolors.OKGREEN}\tData Summary of all aminopropyl:\n'
              f'\t\t# Atoms in aminos` group: {len(self.All_amino_atoms)}\n'
              f'\t\t# Bonds in aminos` group: {len(self.All_amino_bonds)}\n'
              f'\t\t# Angles in aminos` group: {len(self.All_amino_angles)}\n'
              f'\t\t# Dihedrals in aminos` group: '
              f'{len(self.All_amino_dihedrals)}\n'
              f'\t\tTotal charge: {self.All_amino_atoms["charge"].sum()}\n'
              f'\t\tMin charge: {self.All_amino_atoms["charge"].min()}\n'
              f'\t\tMax charge: {self.All_amino_atoms["charge"].max()}'
              f'{bcolors.ENDC}'
              )



class UpdateBoAnDi:
    """update all the bonds, angles, dihedrals, masses"""
    def __init__(self,
                 rot_amino: GetAmino  # Amino data with rotated positions
                 ) -> None:
        self.__update_all(rot_amino)

    def __update_all(self,
                     rot_amino: GetAmino  # Amino data with rotated pos
                     ) -> None:
        """call all the functions"""
        dict_atom_id: dict[int, int]  # to update the indeces
        dict_atom_id = {k: v for k, v in zip(rot_amino.Atoms_df['old_id'],
                                             rot_amino.Atoms_df['atom_id'])}
        self.Bonds_df = self.__update_bonds(
                        rot_amino.Bonds_df, dict_atom_id)
        self.Angles_df = self.__update_angles(
                         rot_amino.Angles_df, dict_atom_id)
        self.Dihedrals_df = self.__update_dihedrals(
                            rot_amino.Dihedrals_df, dict_atom_id)

    def __update_bonds(self,
                       Bonds_df: pd.DataFrame,  # Bonds data frame
                       dict_id: dict[int, int]  # Atoms ids
                       ) -> pd.DataFrame:
        """update atom indexs in the bonds"""
        df = Bonds_df.copy()
        for item, row in Bonds_df.iterrows():
            df.at[item, 'ai'] = dict_id[row['ai']]
            df.at[item, 'aj'] = dict_id[row['aj']]
        return df

    def __update_angles(self,
                        Angles_df: pd.DataFrame,  # Bonds data frame
                        dict_id: dict[int, int]  # Atoms ids
                        ) -> pd.DataFrame:
        """update atom indexs in the angles"""
        df = Angles_df.copy()
        for item, row in Angles_df.iterrows():
            df.at[item, 'ai'] = dict_id[row['ai']]
            df.at[item, 'aj'] = dict_id[row['aj']]
            df.at[item, 'ak'] = dict_id[row['ak']]
        return df

    def __update_dihedrals(self,
                           Dihedrlas_df: pd.DataFrame,  # Bonds data frame
                           dict_id: dict[int, int]  # Atoms ids
                           ) -> pd.DataFrame:
        """update atom indexs in the dihedrals"""
        df = Dihedrlas_df.copy()
        for item, row in Dihedrlas_df.iterrows():
            df.at[item, 'ai'] = dict_id[row['ai']]
            df.at[item, 'aj'] = dict_id[row['aj']]
            df.at[item, 'ak'] = dict_id[row['ak']]
            df.at[item, 'ah'] = dict_id[row['ah']]
        return df


class ConcatAll:
    """append all the aminopropyle to the silicon data file"""
    def __init__(self,
                 silica: upcord.UpdateCoords,  # Silica updated
                 aminos: PrepareAmino  # Update aminos
                 ) -> None:
        self.__concate_all(silica, aminos)
        self.__write_infos()

    def __concate_all(self,
                      silica: upcord.UpdateCoords,  # Silica updated
                      aminos: PrepareAmino  # Update aminos
                      ) -> None:
        """concate the all atoms, bonds, angles, diedrlas"""
        self.Atoms_df: pd.DataFrame  # DF in write_lmp format
        self.Atoms_df = self.__concate_atoms(silica.Atoms_df,
                                             aminos.All_amino_atoms)
        self.Bonds_df = self.__concate_bonds(silica.Bonds_df,
                                             aminos.All_amino_bonds)
        self.Angles_df = self.__concate_angles(silica.Angles_df,
                                               aminos.All_amino_angles)
        self.Dihedrals_df = self.__concate_dihedrals(silica.Dihedrals_df,
                                                     aminos.All_amino_dihedrals
                                                     )
        self.Masses_df = self.__concate_masses(silica.Masses_df,
                                               aminos.Masses_df)
        print(f'{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'Concanating all the atoms\n{bcolors.ENDC}')
        self.__set_attrs()

    def __concate_atoms(self,
                        silica_atoms: pd.DataFrame,  # Silica atoms
                        aminos_atoms: pd.DataFrame  # Aminos atoms
                        ) -> pd.DataFrame:
        """concate all the atoms, make sure they all have same columns"""
        columns: list[str]  # Name of the wanted columns
        columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z',
                   'nx', 'ny', 'nz', 'cmt', 'name']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_atoms[col]
            amino_df[col] = aminos_atoms[col]
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

    def __concate_bonds(self,
                        silica_bonds: pd.DataFrame,  # Silica bonds
                        aminos_bonds: pd.DataFrame  # Aminos bonds
                        ) -> pd.DataFrame:
        """concate all the bonds in the write_lmp format"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'cmt', 'name']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_bonds[col]
            amino_df[col] = aminos_bonds[col]
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

    def __concate_angles(self,
                         silica_angles: pd.DataFrame,  # Silica bonds
                         aminos_angles: pd.DataFrame  # Aminos bonds
                         ) -> pd.DataFrame:
        """concate all the angles in the write_lmp format"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'ak']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_angles[col]
            amino_df[col] = aminos_angles[col]
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

    def __concate_dihedrals(self,
                            silica_dihedrals: pd.DataFrame,  # Silica bonds
                            aminos_dihedrals: pd.DataFrame  # Aminos bonds
                            ) -> pd.DataFrame:
        """concate all the angles in the write_lmp format"""
        columns: list[str]  # Name of the columns
        columns = ['typ', 'ai', 'aj', 'ak', 'ah']
        si_df = pd.DataFrame(columns=columns)
        amino_df = pd.DataFrame(columns=columns)
        for col in columns:
            si_df[col] = silica_dihedrals[col]
            amino_df[col] = aminos_dihedrals[col]
        df: pd.DataFrame  # All atoms dataframe
        df = pd.concat([si_df, amino_df], ignore_index=True)
        df.index += 1
        del si_df
        del amino_df
        return df

    def __concate_masses(self,
                         silica_masses: pd.DataFrame,  # Silica masses
                         amino_masses: pd.DataFrame  # Aminos masses
                         ) -> pd.DataFrame:
        """update and make total masses df"""
        columns: list[str] = ['typ', 'mass', 'cmt', 'name']
        df_silica = pd.DataFrame(columns=columns)
        df_amino = pd.DataFrame(columns=columns)
        for col in columns:
            df_silica[col] = silica_masses[col]
            df_amino[col] = amino_masses[col]

        df = pd.DataFrame(columns=columns)
        df = pd.concat([df_silica, df_amino])
        return df

    def __set_attrs(self) -> None:
        """set attributes to object(self)"""
        self.Nmols: int = np.max(self.Atoms_df['mol'])
        self.NAtoms: int = len(self.Atoms_df)
        self.NBonds: int = len(self.Bonds_df)
        self.NAngles: int = len(self.Angles_df)
        self.NDihedrals: int = len(self.Dihedrals_df)
        self.NAtomTyp: int = np.max(self.Masses_df['typ'])
        self.NBondTyp: int = np.max(self.Bonds_df['typ'])
        self.NAngleTyp: int = np.max(self.Angles_df['typ'])
        self.NDihedralTyp: int = np.max(self.Dihedrals_df['typ'])

    def __write_infos(self) -> None:
        print(f'{bcolors.OKGREEN}\tData Summary after silanization:\n'
              f'\t\t# Atoms: {self.NAtoms}, # Atom`s types: {self.NAtomTyp}\n'
              f'\t\t# Bonds: {self.NBonds}, # Bond`s types: {self.NBondTyp}\n'
              f'\t\t# Angles: {self.NAngles}, '
              f'# Angle`s types: {self.NAngleTyp}\n'
              f'\t\t# Dihedrals: {self.NDihedrals}, '
              f'# Dihedral`s types: {self.NDihedralTyp}\n'
              f'\t\tTotal charge: {self.Atoms_df["charge"].sum()}\n'
              f'\t\tMin charge: {self.Atoms_df["charge"].min()}\n'
              f'\t\tMax charge: {self.Atoms_df["charge"].max()}'
              f'{bcolors.ENDC}'
              )


if __name__ == '__main__':
    fname = sys.argv[1]
    update = upcord.UpdateCoords(fname)  # Updated data for silica
    amino = GetAmino()
    up_aminos = PrepareAmino(update, amino)
    silanized_data = ConcatAll(update, up_aminos)
    np_size: int = int(re.findall(r'\d+', fname)[0])
    fout: str = f'silanized_{np_size}nm.data'
    wrt = wrlmp.WriteLmp(obj=silanized_data, output=fout)
    wrt.write_lmp()
