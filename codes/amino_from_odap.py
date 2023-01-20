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
        print(dihedrals_df)

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
