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
        aligbed_pdb: pd.DataFrame = self.update_df(pdb_src.pdb_df, aligned_xyz)
        self.write_to_pdb(aligbed_pdb, 'aligned_oda.pdb')

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



if __name__ == "__main__":
    AlignOda(sys.argv[1])
