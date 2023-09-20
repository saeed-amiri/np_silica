"""
rotate and order ODA in desired structure in an arbitrary area
"""

import sys
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


if __name__ == "__main__":
    GetPdb(sys.argv[1])
