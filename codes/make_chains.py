import numpy as np
import pandas as pd
import read_lmp_data as rdlmp
import update_coords as upcord
import drop_om as dropOM
from colors_text import TextColor as bcolors


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
        """apply the position update to the aminos' SI & OM"""
        self.Si = 'Si'
        self.OM = 'OM'
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
        Atoms_df = self.__to_origin(amino.Atoms_df)
        Atoms_df = self.__get_azimuths(Atoms_df)

        # Lists to append updated aminos to be replaced
        Atoms_list: list[pd.DataFrame] = []  # Keep all the aminos to concate
        Bonds_list: list[pd.DataFrame] = []  # Keep all the aminos to concate
        Angles_list: list[pd.DataFrame] = []  # Keep all the aminos to concate
        Dihedrals_list: list[pd.DataFrame] = []  # Keep all the aminos> concate

        # Update the type of the atoms before anything
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
             if item <= 2:
                print('OM_list',row['OM_list'], item)
                # calculate the level ups for Aminopropyl
                OM_n: int = 4  # Number of extra atoms (Si, OM) in aminopropyl
                atom_level: int = (item - 1) * (amino.NAtoms - OM_n) + \
                                  update.NAtoms
                mol_level: int = item - 1 + update.Nmols
                
                # Replace the Si atom in Aminopropyl with the proper one
                amino.Atoms_df = self.__set_si_id(Atoms_df.copy(), row)

                # Get atoms info for OM atoms from all NP date
                OM_xyz = self.__get_OM_xyz(row['OM_list'], update.Atoms_df)

                # Repalce the OM atoms info in aminopropyl with proper one
                amino.Atoms_df, amino.Bonds_df, amino.Angles_df, \
                    amino.Dihedrals_df = \
                    self.__set_om_id(amino,
                                     row,
                                     self.__get_azimuths(OM_xyz))

                # Level up the atoms and boandi in aminopropyl data
                i_amino = self.__level_aminos(amino.Atoms_df,
                                                     atom_level,
                                                     mol_level)

                # Rotate each aminopropyl to the direction of the Si
                i_amino = self.__rotate_amino(i_amino)

                # Replace the amino attrs of Atoms with the updated one
                amino.Atoms_df = i_amino

                # Appending all the atoms data
                Atoms_list.append(self.__drop_si_om(i_amino))

                # Update the Bonds/Angles/Dihedrals
                boandi = UpdateBoAnDi(amino)  # Update bonds, angles, dihedrals

                # Append them for concatation
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
        self.Masses_df = self.__drop_si_om_mass(amino_masses)

        print(f'\n{bcolors.OKBLUE}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              f'\t{len(self.All_amino_atoms)} atoms'
              f', {len(self.All_amino_bonds)} bonds'
              f', {len(self.All_amino_angles)} angles'
              f', {len(self.All_amino_dihedrals)} dihedrals'
              f' is updated{bcolors.ENDC}')

    def __set_si_id(self,
                    amino_atoms: pd.DataFrame,  # Amino Atoms information,
                    si_row: pd.DataFrame  # A row of the dataframe
                    ) -> pd.DataFrame:
        """set Si atom_id in amino based on each si_row"""
        amino_atoms['old_id'] = amino_atoms['atom_id']
        column: list[str]  # Columns of the df to replace informations
        column = ['atom_id', 'mol', 'typ', 'x', 'y', 'z', 'rho', 'azimuth',
                  'polar']
        for col in column:
            amino_atoms.at[1, col] = si_row[col]
        return amino_atoms

    def __get_OM_xyz(self,
                     OM_list: list[int],  # Atom id of OM atoms
                     Atoms_df: pd.DataFrame  # Atoms info of NP
                     ) -> pd.DataFrame:
        """return atoms info for OM atoms of each si from NP data file"""
        OM_xyz: list[pd.DataFrame] = []  # Row of each OM atom in NP Atoms_df
        for item in OM_list:
            OM_xyz.append(Atoms_df[Atoms_df['atom_id'] == item])
        return pd.concat(OM_xyz)

    def __set_om_id(self,
                    amino: GetAmino,  # Amino infos
                    si_row: pd.DataFrame,  # One row of Si dataframe
                    OM_xyz: pd.DataFrame  # XYZ info of OM atoms for amino
                    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame,
                               pd.DataFrame]:
        """set the atom ids based on the number of the OM in system"""
        amino_OM: int = 3  # The default numbers of OM in the amino file
        OM_order: int = len(si_row['OM_list'])
        del_OM: int = amino_OM - OM_order
        Atoms_df = amino.Atoms_df.copy()
        Bonds_df = amino.Bonds_df.copy()
        Angles_df = amino.Angles_df.copy()
        Dihedrals_df = amino.Dihedrals_df.copy()
        if del_OM == amino_OM:
            print(f'{bcolors.WARNING}\tNo OM bonded to the SI?\n'
                  f'{bcolors.ENDC}')
        elif del_OM < 0:
            exit(f'{bcolors.FAIL}Error:\n'
                 f'Wrong number of the OM atoms!\n {bcolors.ENDC}')
        elif del_OM == 0:
            Atoms_df = self.__set_OM_info(Atoms_df, si_row, OM_xyz)
            Atoms_df = self.__update_atom_ind(Atoms_df, si_row['OM_list'])
        else:
            print('< 3 : ',si_row['OM_list'], si_row['atom_id'])
            do_om = dropOM.DropOM(amino, si_row, del_OM, OM_xyz)
            Atoms_df = do_om.Atoms_df
            Bonds_df = do_om.Bonds_df
            Angles_df = do_om.Angles_df
            Dihedrals_df = do_om.Dihedrals_df
        return Atoms_df, Bonds_df, Angles_df, Dihedrals_df

    def __drop_si_om(self,
                  amino_atoms: pd.DataFrame  # Rotated Atoms_df
                  ) -> pd.DataFrame:
        """drop the silicon, since it is already in the main data"""
        df: pd.DataFrame = amino_atoms.copy()
        df.drop(amino_atoms[amino_atoms['name'] == 'Si'].index, inplace=True)
        df.drop(amino_atoms[amino_atoms['name'] == 'OM'].index, inplace=True)
        df.reset_index(inplace=True)
        df.drop(columns=['index', 'old_id'], inplace=True, axis=1)
        df.index += 1
        return df

    def __drop_si_om_mass(self,
                       amino_masses: pd.DataFrame  # Amino masses
                       ) -> pd.DataFrame:
        """Drop Si from amino masses, it is only one Si in atoms"""
        df: pd.DataFrame = amino_masses.copy()
        df.drop(amino_masses[amino_masses['name'] == 'Si'].index, inplace=True)
        df.drop(amino_masses[amino_masses['name'] == 'OM'].index, inplace=True)
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
        """update bonds/angels/dihedral type in aminos"""
        silica_ind: int = np.max(silica['typ'])
        df: pd.DataFrame = amino.copy()
        if not np.isnan(silica_ind):
            for item, _ in amino.iterrows():
                df.at[item, 'typ'] += silica_ind
        return df
    
    def __level_aminos(self,
                       Atoms_df: pd.DataFrame, # Amino atoms df
                       atom_level: int,  # To level up the amino atom id
                       mol_level: int   # To level up the amino mol id
                       ) -> pd.DataFrame:
        """update atom id and mol id of atoms in aminopropyl to append
        to silica NP"""
        Si_OM: list[str] = ['Si', 'OM']
        df: pd.DataFrame = Atoms_df.copy()
        root_n: int  # Number of atoms root of the chain (Si-OM)
        root_n = len(Atoms_df[(Atoms_df['name'] == 'Si') | 
                     (Atoms_df['name'] == 'OM')])
        for item, row in Atoms_df.iterrows():
            if row['name'] not in Si_OM:
                df.at[item, 'atom_id'] += atom_level - root_n
                df.at[item, 'mol'] += mol_level
        return df

    def __set_OM_info(self,
                      Atoms_df: pd.DataFrame,  # Atoms of the amino
                      si_row: pd.DataFrame,  # One row of si df
                      OM_xyz: pd.DataFrame  # XYZ info of OM atoms for amino
                      ) -> pd.DataFrame:
        """set info for OM in the amino df"""
        OM_index = [item for item in
                    Atoms_df[Atoms_df['name'] == 'OM'].index]
        df: pd.DataFrame = Atoms_df.copy()
        column: list[str]  # Columns of the df to replace informations
        column = ['atom_id', 'mol', 'typ', 'x', 'y', 'z', 'rho', 'azimuth',
                  'polar']
        for j, k in zip(si_row['OM_list'], OM_index):
            OM_row: pd.DataFrame = OM_xyz[OM_xyz['atom_id'] == j]
            for col in column:
                df.at[k, col] = OM_row[col][j]
            del OM_row
        return df

    def __update_atom_ind(self,
                          df: pd.DataFrame,  # Amino atoms df with removed OM
                          OM_list: list[int]  # Index of the OM in NP
                          ) -> pd.DataFrame:
        """update the index of df after OM index was set"""
        Si_OM_df: pd.DataFrame = df[(df['name'] == 'Si') | (df['name'] == 'OM')]
        Si_OM_len: int = len(Si_OM_df)
        for item, row in df.iterrows():
            if row['name'] != 'Si':
                if row['atom_id'] not in OM_list:
                    df.at[item, 'atom_id'] = item - Si_OM_len
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
                       atoms: pd.DataFrame  # Atoms info of Aminopropyle
                       ) -> pd.DataFrame:
        """calculate the azimuth and polar angle of df"""
        df: pd.DataFrame = atoms.copy()
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
