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

import re
import sys
import typing
import my_tools
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
        atom_id = line[6:11].strip()
        atom_name: str = line[13:16].strip()
        residue_name: str = line[17:20].strip()
        chain_identifier: str = line[21:22]
        residue_number: str = line[22:27].strip()
        x_i: float = float(line[30:38].strip())
        y_i: float = float(line[38:46].strip())
        z_i: float = float(line[46:55].strip())
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
        print(f'\n{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__}):\n'
              f'\tReading water data file: "{water_pdb}"'
              f'{bcolors.ENDC}')


class SetAtomId(ReadWater):
    """Update the atoms id in the all the read df, since some of them
    are str"""
    def __init__(self) -> None:
        super().__init__()
        self.atoms_df: pd.DataFrame = self.update_ids()

    def update_ids(self) -> pd.DataFrame:
        """update all the atoms id to an integer, keep the old id in
        atoms dataframe"""
        return self.__update_atoms(self.atoms_raw.copy())

    def __update_atoms(self,
                       atoms: pd.DataFrame  # Atoms raw dataframe
                       ) -> pd.DataFrame:
        """update atoms"""
        atoms['old_atom_id'] = atoms['chain_identifier'] + atoms['atom_id']
        atoms.index += 1
        atoms['atom_id'] = atoms.index
        return self.__update_resid(atoms)

    def __update_resid(self,
                       atoms: pd.DataFrame  # Atoms df
                       ) -> pd.DataFrame:
        """update the residue id to an integer and keep the old one"""
        residue_name: list[str]  # Name of the residue
        residue_name = [i+j for i, j in zip(atoms['chain_identifier'],
                                            atoms['residue_number'])]
        residue_dict: dict[str, int]  # To make new residue id
        residue_dict = {k: v+1 for v, k in
                        enumerate(self.drop_duplicate(residue_name))}
        residue_id: list[int] = []  # To get residue for all the atoms
        for item in residue_name:
            residue_id.append(residue_dict[item])
        atoms['residue_id'] = residue_id
        return atoms

    def drop_duplicate(self,
                       l_to_set: list[typing.Any]
                       ) -> list[typing.Any]:
        """drop duplicated item with keeping order"""
        seen: set[str] = set()
        seen_add = seen.add
        return [x for x in l_to_set if not (x in seen or seen_add(x))]


class GetWaterDf:
    """make bonds and angles df based on the atoms dataframe
    PDB data file (PDB) has some repeated atom ids and since they are
    not distinguishable, the script makes bonds and angles for all the
    residues in the box."""
    def __init__(self,
                 water_moles: int  # Number of water molecules
                 ) -> None:
        atoms: pd.DataFrame = SetAtomId()
        self.Atoms_df: pd.DataFrame  # In lammps version
        self.Bonds_df: pd.DataFrame  # Updated df
        self.Angles_df: pd.DataFrame  # Updated df
        self.Masses_df: pd.DataFrame  # Masses df
        self.Atoms_df, self.Bonds_df, self.Angles_df = \
            self.make_df(atoms, water_moles)
        self.Masses_df = self.__mk_masses_df()
        self.print_info()

    def make_df(self,
                atoms: SetAtomId,  # updated atoms
                water_resid: int  # Number of water molecules in the box
                ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """make bond and angles df"""
        atoms_lmp: pd.DataFrame = self.__lmp_atoms(atoms.atoms_df)
        bonds: pd.DataFrame = self.__mk_bonds(water_resid)
        bonds.index += 1
        bonds['cmt'] = ['#' for _ in bonds.index]
        bonds['typ'] = [1 for _ in bonds.index]  # Only one bonds' type
        angles: pd.DataFrame = self.__mk_angles(water_resid)
        angles.index += 1
        angles['cmt'] = ['#' for _ in angles.index]
        angles['typ'] = [1 for _ in angles.index]  # Only one angles' type
        return atoms_lmp, bonds, angles

    def __lmp_atoms(self,
                    atoms: pd.DataFrame  # Atoms df
                    ) -> pd.DataFrame:
        """convert data to lammps full atoms format, with extera info"""
        columns: list[str]  # Columns in LAMMPS version
        columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z',
                   'nx', 'ny', 'nz', 'cmt', 'name', 'old_atom_id']
        df_lmp: pd.DataFrame  # ATOMS in lammps format
        df_lmp = pd.DataFrame(columns=columns)
        for col in ['atom_id', 'x', 'y', 'z', 'old_atom_id']:
            df_lmp[col] = atoms[col]
        zero_col: list[int] = [0 for _ in atoms.index]
        for col in ['charge', 'nx', 'ny', 'nz']:
            df_lmp[col] = zero_col
        df_lmp['mol'] = atoms['residue_id']
        df_lmp['typ'] = self.__get_atom_type(atoms)
        df_lmp['cmt'] = ['#' for _ in atoms.index]
        df_lmp['name'] = atoms['atom_name']
        return df_lmp

    def __get_atom_type(self,
                        atoms: pd.DataFrame  # Atoms df
                        ) -> list[int]:
        """make a list of atoms type, since there are only two atoms
        type, just using 1 for H and 2 for O"""
        atom_type: list[int] = []  # Type of the atoms
        atom_names: list[str] = my_tools.drop_duplicate(atoms['atom_name'])
        atom_type_dict = {k: v+1 for v, k in enumerate(atom_names)}
        for item in atoms['atom_name']:
            atom_type.append(atom_type_dict[item])
        return atom_type

    def __mk_bonds(self,
                   resid_max: int  # Number of the residues in the box
                   ) -> pd.DataFrame:
        """make bonds df"""
        df_i: pd.DataFrame = self.__one_bond_df()  # df for one residue
        df_list: list[pd.DataFrame] = []  # To append all the dfs
        for i in range(resid_max):
            id_rise: int = i*3  # 3: Number of atom in molecules
            df_c: pd.DataFrame = df_i.copy()
            df_c['ai'] += id_rise
            df_c['aj'] += id_rise
            df_list.append(df_c)
            del df_c
        return pd.concat(df_list, ignore_index=True)

    def __mk_angles(self,
                    resid_max: int  # Number of the residues in the box
                    ) -> pd.DataFrame:
        """make angles df"""
        df_i: pd.DataFrame = self.__one_angle_df()  # df for one residue
        df_list: list[pd.DataFrame] = []  # To append all the dfs
        for i in range(resid_max):
            id_rise: int = i*3  # 3: Number of atom in molecules
            df_c: pd.DataFrame = df_i.copy()
            df_c['ai'] += id_rise
            df_c['aj'] += id_rise
            df_c['ak'] += id_rise
            df_list.append(df_c)
            del df_c
        return pd.concat(df_list, ignore_index=True)

    def __one_bond_df(self) -> pd.DataFrame:
        """set one bond df"""
        a_i: list[int] = [1, 2]  # 1st atoms in the bonds
        a_j: list[int] = [3, 3]  # 2nd atoms in the bonds
        names: list[str] = ['H-O', 'H-O']
        columns: list[str] = ['ai', 'aj', 'name']  # Name of the columns in df
        df_b_one: pd.DataFrame = pd.DataFrame(columns=columns)
        df_b_one['ai'] = a_i
        df_b_one['aj'] = a_j
        df_b_one['name'] = names
        return df_b_one

    def __one_angle_df(self) -> pd.DataFrame:
        """set one angle df"""
        a_i: list[int] = [1]  # 1st atoms in the angles H
        a_j: list[int] = [3]  # 2nd atoms in the angles O
        a_k: list[int] = [2]  # 3rd atoms in the angles H
        names: list[str] = ['H-O-H']
        columns: list[str] = ['ai', 'aj', 'ak', 'name']  # Columns in df
        df_a_one: pd.DataFrame = pd.DataFrame(columns=columns)
        df_a_one['ai'] = a_i
        df_a_one['aj'] = a_j
        df_a_one['ak'] = a_k
        df_a_one['name'] = names
        return df_a_one

    def __mk_masses_df(self) -> pd.DataFrame:
        """make a df for masses in LAMMPS format"""
        atom_names: list[str] = my_tools.drop_duplicate(self.Atoms_df['name'])
        i_type: int  # Type of atom in the df
        i_row: dict[str, typing.Any]  # Row per atom type
        row_list: list[dict[str, typing.Any]] = []  # All the row to convert
        for item in atom_names:
            name = re.sub('[1-9]', '', item)
            if name not in row_list:
                try:
                    i_type = \
                        list(set(self.Atoms_df[self.Atoms_df['name'] == item]
                                 ['typ']))[0]
                    i_row = {'mass': stinfo.Hydration.MASSES[name],
                             'typ': i_type,
                             'cmt': '#',
                             'name': item,
                             'b_name': item}
                    row_list.append(i_row)
                except KeyError:
                    sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}: '
                             f'({self.__module__})\n'
                             f'\tError! Informaton for `{name}` cannot not be '
                             f'found in `static_info` module\n{bcolors.ENDC}')
        df_m = pd.DataFrame(row_list)
        return df_m

    def print_info(self) -> None:
        """print infos"""
        print(f'\n{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__}):\n'
              f'\tAtoms, Bonds and, Angles are created {bcolors.ENDC}')
