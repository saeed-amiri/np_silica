"""read water box data, output from PACKMOL package created by
water_box module.
The data is in PDB format with bonds and angles written particular
format.
The data is set in exact locations.
atoms:
    ATOM      1  H   HOH A   1     -24.410  81.001 -67.852
    ATOM      2  H   HOH A   1     -25.876  80.363 -68.186
    ATOM      3  O   HOH A   1     -25.259  80.652 -67.455
The second column of ATOM record are obj, could contains strings, but
they are same everwher.
bonds and angles:
    CONECT    1    3
    CONECT    2    3
    CONECT    3    1    2
format of the atoms section:
(I changed the HETAM to ATOM in main file of water.pdb)
Protein Data Bank Format for ATOM:
    Coordinate Section
      Record Type	Columns	Data 	Justification	Data Type
      ATOM 	1-4	“ATOM”	                    	character
      7-11#	Atom serial number      	right	integer
      13-16	Atom name	                left*	character
      17	    Alternate location indicator		character
      18-20§	Residue name	            right	character
      22	    Chain identifier		            character
      23-26	Residue sequence number	    right	integer
      27	    Code for insertions of residues		character
      31-38	X orthogonal Å coordinate	right	real (8.3)
      39-46	Y orthogonal Å coordinate	right	real (8.3)
      47-54	Z orthogonal Å coordinate	right	real (8.3)
      55-60	Occupancy	                right	real (6.2)
      61-66	Temperature factor	        right	real (6.2)
      73-76	Segment identifier¶	        left	character
      77-78	Element symbol              right	character
      79-80	Charge		                        character
For CONECT:
    COLUMNS DATA TYPE FIELD DEFINITION
    ------------------------------------------------------------------
      1  -  6 Record name "CONECT"
      7  - 11 Integer serial Atom serial number
      12 - 16 Integer serial Serial number of bonded atom
      17 - 21 Integer serial Serial number of bonded atom
      22 - 26 Integer serial Serial number of bonded atom
      27 - 31 Integer serial Serial number of bonded atom
      32 - 36 Integer serial Serial number of hydrogen bonded atom
      37 - 41 Integer serial Serial number of hydrogen bonded atom
      42 - 46 Integer serial Serial number of salt bridged atom
      47 - 51 Integer serial Serial number of hydrogen bonded atom
      52 - 56 Integer serial Serial number of hydrogen bonded atom
      57 - 61 Integer serial Serial number of salt bridged atom54
(ref: https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/
      PDB_format_Dec_1996.pdf)
"""

import typing
import pandas as pd
import static_info as stinfo
from colors_text import TextColor as bcolors


class ReadWater:
    """read the PDB file of the water box"""
    def __init__(self) -> None:
        water_pdb: str = stinfo.Hydration.OUT_FILE  # The file to read
        self.atoms_raw: pd.DataFrame  # Water infos with original atom id(raw)
        self.bonds_raw: pd.DataFrame  # Water infos with original atom id(raw)
        self.angles_raw: pd.DataFrame  # Water infos with original atom id(raw)
        self.atoms_raw, self.bonds_raw, self.angles_raw = \
            self.read_pdb(water_pdb)
        self.print_info(water_pdb)

    def read_pdb(self,
                 water_pdb: str  # Name of the file to read
                 ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """read the pdb file line by line"""
        atoms: list[list[typing.Any]] = []  # atoms info
        bonds: list[list[str]] = []  # Bonds info
        angles: list[list[str]] = []  # Angles info
        atoms, bonds, angles = self.__get_infos(water_pdb)
        atoms_df: pd.DataFrame = self.__mk_df(atoms, 'atoms')
        bonds_df: pd.DataFrame = self.__mk_df(bonds, 'bonds')
        angles_df: pd.DataFrame = self.__mk_df(angles, 'angles')
        return atoms_df, bonds_df, angles_df

    def __mk_df(self,
                data: list[typing.Any],  # Data to mk df, atoms|bonds|angles
                d_type: str  # Name of the data section,  atoms|bonds|angles
                ) -> pd.DataFrame:
        """convert data to df and return"""
        columns: list[str]  # Names of the columns of the df
        if d_type == 'atoms':
            columns = ['atom_id',
                       'atom_name',
                       'residue_name',
                       'chain_identifier',
                       'residue_number',
                       'x',
                       'y',
                       'z']
        elif d_type == 'bonds':
            columns = ['a_i', 'a_j']
        elif d_type == 'angles':
            columns = ['a_i', 'a_j', 'a_k']
        df_data: pd.DataFrame = pd.DataFrame(data, columns=columns)
        return df_data

    def __get_infos(self,
                    water_pdb: str  # Name of the file to read
                    ) -> tuple[list[list[typing.Any]],
                               list[list[str]],
                               list[list[str]]]:
        atoms: list[list[typing.Any]] = []  # Save the atoms section
        bonds: list[list[str]] = []  # Save the bonds
        angles: list[list[str]] = []  # Save the angles
        with open(water_pdb, 'r', encoding="utf8") as f_w:
            while True:
                line: str = f_w.readline()
                line_str: str = line.strip()
                if line_str.startswith('ATOM'):
                    # Pass to procces ATOM
                    atoms.append(self.__process_atom(line_str))
                elif line_str.startswith('CONECT'):
                    # Pass to procces CONECT
                    conect: list[str] = self.__process_connect(line_str)
                    if len(conect) == 2:
                        bonds.append(conect)
                    else:
                        angles.append(conect)
                else:
                    pass
                if not line_str:
                    break
        return atoms, bonds, angles

    def __process_atom(self,
                       line: str  # lines which starts with ATOMS
                       ) -> list[typing.Any]:
        """process lines which are starts with ATOM record"""
        print(line)
        atom_id = line[6:11].strip()
        atom_name: str = line[13:16].strip()
        residue_name: str = line[17:20].strip()
        chain_identifier: str = line[21:22]
        residue_number: str = line[22:27].strip()
        x_i: float = float(line[30:38].strip())
        y_i: float = float(line[39:46].strip())
        z_i: float = float(line[47:55].strip())
        return [atom_id,
                atom_name,
                residue_name,
                chain_identifier,
                residue_number,
                x_i,
                y_i,
                z_i]

    def __process_connect(self,
                          line: str  # lines which starts with CONECT
                          ) -> list[str]:
        """"process line which starts with CONECT"""
        a_i: str = line[6:11]
        a_j: str = line[11:16]
        if len(line) > 16:
            a_k: str = line[16:26]
            return [a_i, a_j, a_k]
        return [a_i, a_j]

    def print_info(self,
                   water_pdb: str  # Name of the file to read
                   ) -> None:
        """print infos"""
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__}):\n'
              f'\tReading water data file: {water_pdb}"'
              f'{bcolors.ENDC}')


class SetAtomId(ReadWater):
    """Update the atoms id in the all the read df, since some of them
    are str"""
    def __init__(self) -> None:
        super().__init__()
        self.update_ids()

    def update_ids(self) -> None:
        """update all the atoms id to an integer, keep the old id in
        atoms dataframe"""
        self.__update_atoms(self.atoms_raw.copy())

    def __update_atoms(self,
                       atoms: pd.DataFrame  # Atoms raw dataframe
                       ) -> None:
        """update atoms"""
        atoms.index += 1
        atoms['old_atom_id'] = atoms['atom_id']
        atoms['atom_id'] = atoms.index
        print(atoms)
