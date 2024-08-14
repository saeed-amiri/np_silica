"""
rotate and order ODA in desired structure in an arbitrary area
This script wrote in a couple hours! with help of chatGpt. There is
no test for it. Just did it for one request.
"""

import sys
import numpy as np
import pandas as pd
import logger


class GetPdb:
    """read and load main pdb data"""

    info_msg: str = 'Message from GetPdb:\n'

    def __init__(self,
                 fname: str,  # PDB file
                 log: logger.logging.Logger
                 ) -> None:
        self.fname = fname
        self.info_msg += f'\tReading `{self.fname}` ...\n'
        self.pdb_df: pd.DataFrame = self._initiate()
        log.info(self.info_msg)

    def _initiate(self) -> pd.DataFrame:
        data_lines: list[str] = self.get_data()
        parsed_lines = self.parse_lines(data_lines)
        return self.mk_dataframe(parsed_lines)

    def get_data(self) -> list[str]:
        """read pdb file"""
        with open(self.fname, 'r', encoding='utf8') as file:
            lines: list[str] = file.readlines()

        # Extract x and y data
        data_lines: list[str] = \
            [line for line in lines if line.startswith("ATOM") and
             line != '&\n']
        return data_lines

    @staticmethod
    def parse_lines(data_lines: list[str]
                    ) -> list[list[str]]:
        """split the lines"""
        splited_line: list[list[str]] = \
            [line.strip().split(' ') for line in data_lines]
        return [[i for i in item if i] for item in splited_line]

    @staticmethod
    def mk_dataframe(lines: list[list[str]]
                     ) -> pd.DataFrame:
        """make dataframe based on the columns name"""
        columns: list[str] = ['ATOM',
                              'atom_id',
                              'atom_name',
                              'residue_name',
                              'residue_number',
                              'x',
                              'y',
                              'z',
                              'temp',
                              'atom_symbol']
        return pd.DataFrame(lines, columns=columns)


class AlignOda:
    """
    Align Oda to a desired axis
    """

    info_msg: str = "Message from AlignOda:\n"
    fout: str = 'aligned_oda.pdb'
    align_axis: str = 'z'

    def __init__(self,
                 fname: str,
                 log: logger.logging.Logger
                 ) -> None:
        pdb_src = GetPdb(fname, log)
        self._initiate(pdb_src, log)

    def _initiate(self,
                  pdb_src: "GetPdb",
                  log: logger.logging.Logger
                  ) -> None:
        oda_xyz: np.ndarray = self.get_xyz(pdb_src.pdb_df)
        aligned_xyz = self.align_to_z_axis(oda_xyz)
        aligned_pdb: pd.DataFrame = \
            self.update_df(pdb_src.pdb_df, aligned_xyz)
        self.aligbed_pdb = self.check_nitrogen_is_down(aligned_pdb)

        self.write_to_pdb(self.aligbed_pdb, self.fout)
        self.info_msg += f'\tOda aligned along `{self.align_axis}` axis\n'
        self.info_msg += f'\tOutput file is writen in `{self.fout}` \n'
        log.info(self.info_msg)
        self.info_msg = ''

    def get_xyz(self,
                pdb_df: pd.DataFrame
                ) -> np.ndarray:
        """
        get coordinates as floats
        """
        columns: list[str] = ['x', 'y', 'z']
        return pdb_df[columns].astype(float).to_numpy()

    def align_to_z_axis(self,
                        oda_xyz: np.ndarray
                        ) -> np.ndarray:
        """
        aligning along axis
        """
        # Compute the centroid of the molecule
        centroid: np.ndarray = np.mean(oda_xyz, axis=0)

        # Center the molecule at the origin
        centered_xyz: np.ndarray = oda_xyz - centroid

        # Compute the inertia tensor
        inertia_tensor: np.ndarray = np.dot(centered_xyz.T, centered_xyz)

        # Compute the eigenvalues and eigenvectors of the inertia tensor
        eigvals, eigvecs = np.linalg.eig(inertia_tensor)
        # Sort eigenvectors by eigenvalues in ascending order
        sorted_indices: np.ndarray = np.argsort(eigvals)
        # Depending on which axis you want to align to:

        # For x-axis: choose the eigenv correspond to largest eigenvalue
        # For y-axis: choose the eigenv correspond to second largest eigenvalue
        # For z-axis: choose the eigenv correspond to smallest eigenvalue
        if self.align_axis == 'z':
            rotation_matrix = eigvecs[:, sorted_indices]
        if self.align_axis == 'y':
            rotation_matrix = eigvecs[:, sorted_indices].T
            rotation_matrix[[0, 1]] = rotation_matrix[[1, 0]]
        elif self.align_axis == 'x':
            rotation_matrix = eigvecs[:, sorted_indices].T

        # Ensure that the rotation matrix is right-handed
        if np.linalg.det(rotation_matrix) > 0:
            rotation_matrix[:, -1] *= -1

        # Rotate molecule to align with z-axis
        return np.dot(centered_xyz, rotation_matrix)

    def update_df(self,
                  pdb_df: pd.DataFrame,
                  aliged_xyz: np.ndarray
                  ) -> pd.DataFrame:
        """replace source dataframe with updated xyz poitions"""
        columns_to_replace = ['x', 'y', 'z']
        aliged_pdb: pd.DataFrame = pdb_df.copy()
        aliged_pdb[columns_to_replace] = aliged_xyz
        return aliged_pdb

    def check_nitrogen_is_down(self,
                               pdb_df: pd.DataFrame
                               ) -> pd.DataFrame:
        """check the nitrogen is down"""
        # Check if 'NH3' is among the four atoms with the minimum 'z' values
        sorted_df = pdb_df.sort_values(by='z')
        top_four_atoms = sorted_df.head(4)
        if 'NH2' not in top_four_atoms['atom_name'].values:
            print("Nitrogen is already pointing down")
            columns: list[str] = ['x', 'y', 'z']
            oda_xyz: np.ndarray = \
                pdb_df[columns].astype(float).to_numpy()
            # Rotate the system by 180 degrees (pi radians)
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, -1]
            ])
            oda_xyz = np.dot(oda_xyz, rotation_matrix)
            pdb_df[columns] = oda_xyz
        return pdb_df

    @staticmethod
    def write_to_pdb(pdb_df: pd.DataFrame,
                     filename: str) -> None:
        """
        write as pdb format
        """
        with open(filename, 'w', encoding='utf8') as f_out:
            for _, row in pdb_df.iterrows():
                atom_line = \
                    (f"ATOM  {int(row['atom_id']):5d} "
                     f"{row['atom_name']:<4s} {row['residue_name']:>3s}  "
                     f"{int(row['residue_number']):>4d}    "
                     f"{float(row['x']):8.3f}{float(row['y']):8.3f}"
                     f"{float(row['z']):8.3f}{float(row['temp']):6.2f}      "
                     f"{row['atom_symbol']:>2s}\n")
                f_out.write(atom_line)


class OrderOda(AlignOda):
    """
    Order the oda in a hexagonal structure in a squre area
    In order to displace z along z axis, there are 3 options:
        normal: z value is moved by z_offset, for all the ODAs,
        zigzag: odd and even rows gets to different values. odd moves
                by z_offset and even remain at zero,
        random: Every Oda move by a random values (zero, abs(z_offset))
    It is just for case of carshes in the running the flat ODA lattice!
    """

    spacing: float = 9.45
    scale_x: float = 1.4
    radius: float = 0
    z_offset: float = 10
    z_style: str = 'random'
    a_x: float = 429
    a_y: float = 429
    desired_oda_nr: int = 200
    output: str = f'ordered_ODA_{z_style}.pdb'

    def __init__(self,
                 fname: str,
                 log: logger.logging.Logger
                 ) -> None:
        super().__init__(fname, log)
        self.info_msg: str = 'Messages From OrderOda:\n'
        self._loging_attributs()
        self.mk_structure()
        log.info(self.info_msg)

    def _loging_attributs(self) -> None:
        """just to save the constant to the log file"""
        self.info_msg += (
            "\tConstant values used to set the structure:\n"
            f"\t\tSize of box in x is `{self.a_x}`, in y is `{self.a_y}`\n"
            f"\t\tThe radius for the nanoparticle is set to `{self.radius}`\n"
            f"\t\tDesired Numbers of ODA is `{self.desired_oda_nr}`\n"
            f"\t\tz_offset is `{self.z_offset}`\n"
            f"\t\tz_offset style is `{self.z_style}`\n"
            f"\t\tSpacing is `{self.spacing}`\n"
            f"\t\tScale_x is `{self.scale_x}`\n"
            f"\t\tOutput file of the final orderd ODA is `{self.output}`\n\n"
        )

    def mk_structure(self) -> None:
        """
        prepare the area
        """
        columns: list[str] = ['x', 'y', 'z']
        oda_xyz: np.ndarray = \
            self.aligbed_pdb[columns].astype(float).to_numpy()
        n_x, n_y = self.calculate_nx_ny(self.a_x, self.a_y)
        z_offset_arr: np.ndarray = self.mk_z_offset(n_x, n_y)
        lattice_points: list[tuple[float, ...]] = \
            self.generate_hexagonal_lattice(n_x, n_y, z_offset_arr)
        if self.radius >= 0:
            excluded_lattice: list[tuple[float, ...]] = \
                self.remove_points_inside_np(lattice_points, self.radius)
            lattice_with_desired_points: list[tuple[float, ...]] =\
                self.adjust_oda_count(excluded_lattice)
        else:
            excluded_lattice = lattice_with_desired_points = \
                self.remove_points_outside_radius(lattice_points, self.radius)
        msg: str = f'Number of ODA is `{len(lattice_with_desired_points)}`'
        self.info_msg += f'\t{msg}\n'
        print(msg)
        oda_lattice = self.position_molecule_on_lattice(
            oda_xyz, lattice_with_desired_points)
        new_df: pd.DataFrame = self.repeat_dataframe(
            self.aligbed_pdb, len(lattice_with_desired_points))
        new_pdb: pd.DataFrame = self.update_df(new_df, np.vstack(oda_lattice))
        self.write_to_pdb(new_pdb, self.output)

    def mk_z_offset(self,
                    n_x: int,
                    n_y: int
                    ) -> np.ndarray:
        """create a n_x by n_y array to save the z_offset"""
        z_offset_arr: np.ndarray = np.zeros((n_x, n_y))
        if (style := self.z_style) == 'normal':
            return z_offset_arr + self.z_offset
        if style == 'zigzag':
            z_offset_arr[::2, :] += self.z_offset
            return z_offset_arr
        if style == 'random':
            return np.random.uniform(0, abs(self.z_offset), (n_x, n_y))
        return z_offset_arr

    def calculate_nx_ny(self,
                        a_x: float,
                        a_y: float
                        ) -> tuple[int, int]:
        """find the number of points in the square area"""
        # Rough estimate without considering exclusion zone
        n_x_rough = int(a_x / (self.spacing * self.scale_x))
        n_y_rough = int(a_y / (self.spacing * np.sqrt(3)))
        lattice_points = \
            self.generate_hexagonal_lattice(n_x_rough,
                                            n_y_rough,
                                            np.zeros((n_x_rough, n_y_rough)))

        # Calculate final n_x and n_y
        all_x = [point[0] for point in lattice_points]
        all_y = [point[1] for point in lattice_points]

        n_x = int(np.ceil(max(all_x) / (self.spacing * self.scale_x)))
        n_y = int(np.ceil(max(all_y) / (self.spacing * np.sqrt(3))))

        return n_x, n_y

    @staticmethod
    def remove_points_inside_np(lattice_points: list[tuple[float, ...]],
                                radius: float
                                ) -> list:
        """
        exclude the points which are inside the NP
        """
        excluded_lattice: list[tuple[float, ...]] = []
        for point in lattice_points:
            r_point = point[0]**2 + point[1]**2
            if r_point > radius**2:
                excluded_lattice.append(point)
        return excluded_lattice

    def adjust_oda_count(self,
                         lattice_points: list[tuple[float, ...]]
                         ) -> list[tuple[float, ...]]:
        """check the number and drop extra points if needed"""
        if (oda_nr := len(lattice_points)) < self.desired_oda_nr:
            msg = ("Number of the ODA is less than expected value: "
                   f"`{oda_nr} < {self.desired_oda_nr}`")
            self.info_msg += f'\t{msg}\n'
            print(msg)
            return lattice_points
        if oda_nr == self.desired_oda_nr:
            msg = "The expected number of ODA is created!"
            self.info_msg += f'\t{msg}\n'
            print(msg)
            return lattice_points
        return self._drop_lattice_points(
                lattice_points, oda_nr-self.desired_oda_nr)

    def _drop_lattice_points(self,
                             lattice_points: list[tuple[float, ...]],
                             drop_nr: int
                             ) -> list[tuple[float, ...]]:
        """drop extra ODA from list"""
        msg: str = f"`{drop_nr}` ODA dropped to get the expected number!"
        self.info_msg += f'\t{msg}\n'
        print(msg)
        xyz_arr: np.ndarray = np.array(lattice_points)
        # Compute the Euclidean distance for each point to the origin
        distances = np.linalg.norm(xyz_arr, axis=1)

        # Find the indices of the n smallest distances
        indices_to_remove: np.ndarray = \
            np.argpartition(distances, drop_nr)[:drop_nr]

        # Remove these indices from the array
        filtered_arr: np.ndarray = \
            np.delete(xyz_arr, indices_to_remove, axis=0)
        return [tuple(row) for row in filtered_arr]

    def remove_points_outside_radius(self,
                                     lattice_points: list[tuple[float, ...]],
                                     radius: float
                                     ) -> list[tuple[float, ...]]:
        """
        exclude the points which are outside the radius
        """
        excluded_lattice: list[tuple[float, ...]] = []
        for point in lattice_points:
            r_point = point[0]**2 + point[1]**2
            if r_point < radius**2:
                excluded_lattice.append(point)
        return excluded_lattice

    def generate_hexagonal_lattice(self,
                                   n_x: int,  # Number in x dirction
                                   n_y: int,
                                   z_offset_arr: np.ndarray
                                   ) -> list[tuple[float, ...]]:
        """
        Generate a hexagonal lattice of size n x m with lattice constant a.

        Returns:
        - list: List of (x, y, z) coordinates for the lattice sites.
        """
        lattice_points: list[tuple[float, ...]] = []
        for i in range(n_x):
            for j in range(n_y):
                if i % 2 == 0:
                    x_i = self.spacing * j * self.scale_x
                else:
                    x_i = self.spacing * j * self.scale_x + 0.5 * self.spacing
                y_i = self.spacing * np.sqrt(3) * i
                z_i = z_offset_arr[i, j]
                lattice_points.append((x_i, y_i, z_i))
        # Calculate centroid of the lattice points
        centroid = np.mean(lattice_points, axis=0)

        # Center the lattice points to (0, 0, 0)
        centered_lattice_points: list[tuple[float, ...]] = \
            [(x-centroid[0], y-centroid[1], z) for x, y, z in lattice_points]
        return centered_lattice_points

    def position_molecule_on_lattice(self,
                                     molecule: np.ndarray,
                                     lattice_points: list[tuple[float, ...]]
                                     ) -> list[np.ndarray]:
        """
        Position the molecule at each site of the lattice.

        Returns:
        - list: List of molecules positioned at each lattice site.
        """
        # Calculate the centroid of the molecule
        centroid: np.ndarray = np.mean(molecule, axis=0)

        positioned_molecules: list[np.ndarray] = []
        for point in lattice_points:
            translation = np.array(point) - centroid
            translated_molecule = molecule + translation
            positioned_molecules.append(translated_molecule)

        return positioned_molecules

    @staticmethod
    def repeat_dataframe(dframe: pd.DataFrame,
                         n_rep: int
                         ) -> pd.DataFrame:
        """repeat the main dataframe"""
        # Ensure 'atom_id' and 'residue_number' are integers
        dframe = dframe.astype({'atom_id': int, 'residue_number': int})

        row_count_original: int = len(dframe)

        tmp_list: list[pd.DataFrame] = []
        for i in range(n_rep):
            dframe_i = dframe.copy()
            dframe_i['residue_number'] += 1 * i
            dframe_i['atom_id'] += i * row_count_original
            tmp_list.append(dframe_i)

        new_dframe: pd.DataFrame = pd.concat([
            dframe.assign(
                residue_number=dframe['residue_number'] + i,
                atom_id=dframe['atom_id'] + i * row_count_original
            )
            for i in range(n_rep)
        ], ignore_index=True)

        return new_dframe


if __name__ == "__main__":
    OrderOda(sys.argv[1], log=logger.setup_logger(log_name='ordered_oda.log'))
