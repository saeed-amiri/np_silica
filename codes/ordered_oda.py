"""
rotate and order ODA in desired structure in an arbitrary area
"""

import sys
import numpy as np
import pandas as pd


class GetPdb:
    """read and load main pdb data"""
    def __init__(self,
                 fname: str  # PDB file
                 ) -> None:
        self.fname = fname
        self.pdb_df: pd.DataFrame = self._initiate()

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
    def __init__(self,
                 fname: str
                 ) -> None:
        pdb_src = GetPdb(fname)
        self._initiate(pdb_src)

    def _initiate(self,
                  pdb_src: "GetPdb"
                  ) -> None:
        oda_xyz: np.ndarray = self.get_xyz(pdb_src.pdb_df)
        aligned_xyz = self.align_to_z_axis(oda_xyz)
        self.aligbed_pdb: pd.DataFrame = \
            self.update_df(pdb_src.pdb_df, aligned_xyz)
        self.write_to_pdb(self.aligbed_pdb, 'aligned_oda.pdb')

    def get_xyz(self,
                pdb_df: pd.DataFrame
                ) -> np.ndarray:
        """
        get coordinates as floats
        """
        columns: list[str] = ['x', 'y', 'z']
        return pdb_df[columns].astype(float).to_numpy()

    def align_to_z_axis(self, oda_xyz, align_axis='z'):
        """
        aligning along axis
        """
        # Compute the centroid of the molecule
        centroid = np.mean(oda_xyz, axis=0)

        # Center the molecule at the origin
        centered_xyz = oda_xyz - centroid

        # Compute the inertia tensor
        inertia_tensor = np.dot(centered_xyz.T, centered_xyz)

        # Compute the eigenvalues and eigenvectors of the inertia tensor
        eigvals, eigvecs = np.linalg.eig(inertia_tensor)
        # Sort eigenvectors by eigenvalues in ascending order
        sorted_indices = np.argsort(eigvals)
        # Depending on which axis you want to align to:

        # For x-axis: choose the eigenv correspond to largest eigenvalue
        # For y-axis: choose the eigenv correspond to second largest eigenvalue
        # For z-axis: choose the eigenv correspond to smallest eigenvalue
        if align_axis == 'z':
            rotation_matrix = eigvecs[:, sorted_indices]
        if align_axis == 'y':
            rotation_matrix = eigvecs[:, sorted_indices].T
            rotation_matrix[[0, 1]] = rotation_matrix[[1, 0]]
        elif align_axis == 'x':
            rotation_matrix = eigvecs[:, sorted_indices].T

        # Ensure that the rotation matrix is right-handed
        if np.linalg.det(rotation_matrix) < 0:
            rotation_matrix[:, -1] *= -1

        # Rotate molecule to align with z-axis
        return np.dot(centered_xyz, rotation_matrix)

    @staticmethod
    def update_df(pdb_df: pd.DataFrame,
                  aliged_xyz: np.ndarray
                  ) -> pd.DataFrame:
        """replace source dataframe with updated xyz poitions"""
        columns_to_replace = ['x', 'y', 'z']
        aliged_pdb: pd.DataFrame = pdb_df.copy()
        aliged_pdb[columns_to_replace] = aliged_xyz
        return aliged_pdb

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
                     f"{row['atom_name']:<4s}{row['residue_name']:>3s}  "
                     f"{int(row['residue_number']):>4d}    "
                     f"{float(row['x']):8.3f}{float(row['y']):8.3f}"
                     f"{float(row['z']):8.3f}{float(row['temp']):6.2f}      "
                     f"{row['atom_symbol']:>2s}\n")
                f_out.write(atom_line)


class OrderOda(AlignOda):
    """
    Order the oda in a hexagonal structure in a squre area
    """

    spacing: float = 4.5
    scale_x: float = 1.5
    radius: float = 35
    z_offset: float = 0

    def __init__(self,
                 fname: str
                 ) -> None:
        super().__init__(fname)
        self.mk_structure()

    def mk_structure(self) -> None:
        """
        prepare the area
        """
        columns: list[str] = ['x', 'y', 'z']
        oda_xyz: np.ndarray = \
            self.aligbed_pdb[columns].astype(float).to_numpy()
        n_x, n_y = self.calculate_nx_ny(a_x=220, a_y=220)
        lattice_points: list[tuple[float, ...]] = \
            self.generate_hexagonal_lattice(n_x, n_y)
        excluded_lattice: list[tuple[float, ...]] = \
            self.exclude_np_zrea(lattice_points, self.radius)
        print(f'Number of ODA is {len(excluded_lattice)}')
        oda_lattice = \
            self.position_molecule_on_lattice(oda_xyz, excluded_lattice)
        new_df: pd.DataFrame = \
            self.repeat_dataframe(self.aligbed_pdb, len(excluded_lattice))
        new_pdb: pd.DataFrame = self.update_df(new_df, np.vstack(oda_lattice))
        self.write_to_pdb(new_pdb, 'structured.pdb')

    def calculate_nx_ny(self,
                        a_x: float,
                        a_y: float
                        ) -> tuple[int, int]:
        """find the number of points in the square area"""
        # Rough estimate without considering exclusion zone
        n_x_rough = int(np.ceil(a_x / (self.spacing * self.scale_x)))
        n_y_rough = int(np.ceil(a_y / (self.spacing * np.sqrt(3))))
        lattice_points = \
            self.generate_hexagonal_lattice(n_x_rough, n_y_rough)

        # Calculate final n_x and n_y
        all_x = [point[0] for point in lattice_points]
        all_y = [point[1] for point in lattice_points]

        n_x = int(np.ceil(max(all_x) / (self.spacing * self.scale_x)))
        n_y = int(np.ceil(max(all_y) / (self.spacing * np.sqrt(3))))

        return n_x, n_y

    @staticmethod
    def exclude_np_zrea(lattice_points: list[tuple[float, ...]],
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

    def generate_hexagonal_lattice(self,
                                   n_x: int,  # Number in x dirction
                                   n_y: int
                                   ) -> list[tuple[float, ...]]:
        """
        Generate a hexagonal lattice of size n x m with lattice constant a.

        Parameters:
        - n (int): Number of rows.
        - m (int): Number of columns.
        - a (float): Lattice constant.
        - z_offset (float): Vertical offset.

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
                z_i = self.z_offset
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
        dframe['atom_id'] = dframe['atom_id'].astype(int)
        dframe['residue_number'] = dframe['residue_number'].astype(int)

        # Create a new dataframe by repeating the original dataframe `n` times
        new_dframe = pd.concat([dframe] * n_rep, ignore_index=True)

        # Compute offset for 'atom_id' & 'residue_number' from repeat number
        ids = (new_dframe.index // len(dframe) * len(dframe)).values
        new_dframe['atom_id'] += ids
        new_dframe['residue_number'] += ids

        return new_dframe


if __name__ == "__main__":
    OrderOda(sys.argv[1])
