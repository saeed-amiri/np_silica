import sys
import pandas as pd
import read_lmp_data as redlmp
import write_lmp as wrlmp


class Doc:
    """read odap data file and return all the atomic information
    which are bigger then 47 to make the chain of CH2-CH2-CH2-NH2
    it is a naive way to do it, but must reliable for now
    """


class GetAmino:
    """make 3Aminopropyl from oda"""
    def __init__(self) -> None:
        self.auditor: int = 47  # The atom which the chain strat with
        self.get_data()

    def get_data(self) -> None:
        fname: str  # Name of the odap.data in LAMMPS full atom data
        fname = '/scratch/saeed/MyScripts/np_silica/data/ODAP.data'
        data = self.get_odap(fname)
        atoms_df = self.__get_atoms(data.Atoms_df)
        bonds_df = self.__get_bonds(data.Bonds_df)
        angles_df = self.__get_angles(data.Angles_df)
        dihedrals_df = self.__get_dihedrals(data.Dihedrals_df)
        self.treat_atoms(data.Masses_df, atoms_df)

    def treat_atoms(self,
                    masses_df: pd.DataFrame,  # Masses atomic infos
                    atoms_df: pd.DataFrame  # Selected data from input
                    ) -> pd.DataFrame:
        """update atoms type and indeces"""
        self.update_atom_id(atoms_df)

    def update_atom_id(self,
                       atoms_df: pd.DataFrame  # Selected data from input
                       ) -> pd.DataFrame:
        """add new index"""
        up_atoms: pd.DataFrame = atoms_df.copy()
        up_atoms['old_id'] = up_atoms['atom_id']
        new_id: list[int]  # New atoms' id
        new_id = list(range(1, len(up_atoms)+1))
        up_atoms['atom_id'] = new_id
        up_atoms.reset_index(inplace=True)
        up_atoms.drop('index', inplace=True, axis=1)
        up_atoms.index += 1
        return up_atoms

    def __get_atoms(self,
                    atoms: pd.DataFrame  # Atoms LAMMPS format
                    ) -> pd.DataFrame:
        return atoms[atoms['atom_id'] >= self.auditor]

    def __get_bonds(self,
                    bonds: pd.DataFrame  # Bonds LAMMPS format
                    ) -> pd.DataFrame:
        return bonds[(bonds['ai'] >= self.auditor) &
                     (bonds['aj'] >= self.auditor)]

    def __get_angles(self,
                     angles: pd.DataFrame  # Angles LAMMPS format
                     ) -> pd.DataFrame:
        return angles[(angles['ai'] >= self.auditor) &
                      (angles['aj'] >= self.auditor) &
                      (angles['ak'] >= self.auditor)]

    def __get_dihedrals(self,
                        dihedrals: pd.DataFrame  # Dihedrals LAMMPS format
                        ) -> pd.DataFrame:
        return dihedrals[(dihedrals['ai'] >= self.auditor) &
                         (dihedrals['aj'] >= self.auditor) &
                         (dihedrals['ak'] >= self.auditor) &
                         (dihedrals['ah'] >= self.auditor)]

    def get_odap(self,
                 fname: str  # Name of the odap file
                 ) -> redlmp:
        """read odap file"""
        data = redlmp.ReadData(fname)
        return data


if __name__ == '__main__':
    odap = GetAmino()
