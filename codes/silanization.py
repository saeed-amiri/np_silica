import sys
import numpy as np
import pandas as pd
import read_lmp_data as rdlmp
import update_coords as upcord
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
        """apply the position update to the aminos"""
        self.__update_aminos(update, amino)

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
        for item, row in si_df.iterrows():
            if item < 8:
                # Si from amino will be deleted later, so the rest of
                # atoms must start on lower
                atom_level: int = (item - 1) * \
                                  (amino.NAtoms - 1) + update.NAtoms - 1
                mol_level: int = item + update.Nmols
                i_amino = self.__upgrade_amino(Atoms_df.copy(),
                                               row,
                                               atom_level,
                                               mol_level)
                update.Atoms_df = i_amino
                self.__update_BoAnDiMa(update)
            else:
                break

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

    def __rotate_amino(self,
                       amino_atoms: pd.DataFrame  # Updated amino
                       ) -> pd.DataFrame:
        """rotate amino based on the azimuth and polar of the Si"""
        # print(amino_atoms)
        ref_x: float = amino_atoms['x'][1]
        ref_y: float = amino_atoms['y'][1]
        ref_z: float = amino_atoms['z'][1]
        new_x: list[float] = []  # After rotations
        new_y: list[float] = []  # After rotations
        new_z: list[float] = []  # After rotations
        beta: float = amino_atoms['azimuth'][1]
        gamma: float = amino_atoms['polar'][1]
        new_x.append(amino_atoms['x'][1])
        new_y.append(amino_atoms['y'][1])
        new_z.append(amino_atoms['z'][1])
        for item, row in amino_atoms.iterrows():
            if item > 1:
                x, y, z = self.__rot_matrix(beta, gamma,
                                            row['x'], row['y'], row['z'])
                new_x.append(x + ref_x)
                new_y.append(y + ref_y)
                new_z.append(z + ref_z)
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
                     ) -> tuple[float]:
        """calculate the rotated coordes"""
        cos_g: float = np.cos(gamma)
        cos_b: float = np.cos(beta)
        sin_g: float = np.sin(gamma)
        sin_b: float = np.sin(beta)

        x_new: float = x*cos_g*cos_b - y*sin_g + z*sin_b*cos_g
        y_new: float = x*sin_g*cos_b + y*cos_g + z*sin_b*sin_g
        z_new: float = -x**sin_b + z*cos_b
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
        print(f'{bcolors.OKCYAN}\tMove Aminopropyle [Si] to origin\n'
              f'{bcolors.ENDC}')
        x_si: float = amino_atoms[amino_atoms['name'] == self.Si]['x'][1]
        y_si: float = amino_atoms[amino_atoms['name'] == self.Si]['y'][1]
        z_si: float = amino_atoms[amino_atoms['name'] == self.Si]['z'][1]
        df['x'] -= x_si
        df['y'] -= y_si
        df['z'] -= z_si
        return df


class UpdateBoAnDiMa:
    """update all the bonds, angles, dihedrals, masses"""
    def __init__(self,
                 rot_amino: PrepareAmino  # Amino data with rotated positions
                 ) -> None:
        self.__update_all(rot_amino)

    def __update_all(self,
                     rot_amino: PrepareAmino  # Amino data with rotated pos
                     ) -> None:
        """call all the functions"""


if __name__ == '__main__':
    fname = sys.argv[1]
    update = upcord.UpdateCoords(fname)
    amino = GetAmino()
    up_aminos = PrepareAmino(update, amino)
