"""convert the Atoms section into PDB file
Input:
    LAMMPS data file from lmp_itp_pdb.py:
        Atoms_df, Masses_df
Output:
    pd.DataFrame for pdb file

Main notes:
    ? The names of atoms should be based on the force field:
    > can be set in the input files in the Masses section after a `#.`
    ? The names of the Residues (molecules) should be set somehow
    since, in the lammps, they just separated by numbers:
    > The second name after the atom name in the Masses section
    should be the name of the residue to which the atom belongs
    and should be three or four characters long.
    - Limitations for atoms' names also should be considered.
    ? The element symbols should also be considered:
    The easiest way is to set it in the Masses section after the
    residues' names and in CAPTIAL letters.
    ? The ATOMS or HATOM:
    > It also set after the symbols in Masses section.
    ? TER: This shows that the last atom in the chain prevents the
    display of a connection to the next chain. It specified by the
    index of the aotm:
    > Bit tricky! Do it in the code by selecting the final atom as
    TER if the it is ATOM in the Masses section. For my input should
    work since the last atom is the final atom.
    ? Duplicate Atom Names: One possible editing mistake is the
    failure to uniquely name all atoms within a given residue.
    > Do it in the code by adjusting duplicate atom names by adding
    an index.

    - The following recordes will be left empty chars:
        17	    Alternate location indicator		character
        22	    Chain identifier		            character
        27	    Code for insertions of residues		character
        55-60	Occupancy	                right	real (6.2)
        61-66	Temperature factor	        right	real (6.2)
        73-76	Segment identifier¶	        left	character
    ? For Occupancy:
    > equal to 1, see:
    `https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/
    dealing-with-coordinates`
######################################################################

PdbStyleInfo:
    The PDB file consider is in standard PDB format.

    convert LAMMPS data file to a standard PDB file format based on:
[https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html]
    A PDB file has a format as:
    ATOM 	atomic coordinate record containing the X,Y,Z
        orthogonal Å coordinates for atoms in standard residues]
        (amino acids and nucleic acids).

    HATATM 	atomic coordinate record containing the X,Y,Z orthogonal Å
    coordinates for atoms in nonstandard residues. Nonstandard residues
    include inhibitors, cofactors, ions, and solvent. The only functional
    difference from ATOM records is that HETATM residues are by default
    not connected to other residues. Note that water residues should be
    in HETATM records.

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

    Protein Data Bank Format for HATATM:
        Coordinate Section
        Record Type	Columns	Data 	Justification	Data Type
        1-6	“HETATM”		character
        7-80	same as ATOM records

    #Chimera allows (nonstandard) use of columns 6-11 for the integer
        atom serial number in ATOM records, and in TER records, only the
        “TER” is required.
    *Atom names start with element symbols right-justified in columns
        13-14 as permitted by the length of the name. For example, the
        symbol FE for iron appears in columns 13-14, whereas the symbol
        C for carbon appears in column 14 (see Misaligned Atom Names).
        If an atom name has four characters, however, it must start in
        column 13 even if the element symbol is a single character
        (for example, see Hydrogen Atoms).
    §Chimera allows (nonstandard) use of four-character residue names
        occupying an additional column to the right.
    ¶Segment identifier is obsolete, but still used by some programs.
        Chimera assigns it as the atom attribute pdbSegment to allow
        command-line specification.
    The format of ecah section is (fortran style):
    Format (A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)
"""

import sys
import typing
from collections import Counter
import pandas as pd
from colors_text import TextColor as bcolors


class Pdb:
    """convert to PDB with considers all the above concerns"""
    def __init__(self,
                 Masses_df: pd.DataFrame,  # df contains info about atoms
                 Atoms_df: pd.DataFrame  # df contains atoms' coordinats
                 ) -> None:
        """call the main functions"""
        self.__show_warnings()
        # Sort Atoms_df based on the atom index
        Atoms_df.sort_values(by=['atom_id'], inplace=True)
        self.mk_pdb(Masses_df, Atoms_df)
        self.print_info()

    def mk_pdb(self,
               Masses_df: pd.DataFrame,  # df contains info about atoms
               Atoms_df: pd.DataFrame  # df contains atoms' coordinats
               ) -> None:
        """make pdb data structure"""
        Masses = self.__get_atoms_info(Masses_df)
        pdb_df: pd.DataFrame = self.__mk_pdb_df()  # Empty df
        self.pdb_df = self.__set_pdb(pdb_df, Masses, Atoms_df)
        print(f'{bcolors.WARNING}\tTotal charge is : '
              f'{self.pdb_df["q"].sum():.4f}{bcolors.ENDC}')

    def __show_warnings(self) -> None:
        """show warnings and infos"""
        print(f'{bcolors.CAUTION}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              '\tMasses section in the input file should be in the '
              'following order:\n'
              '\t\tid mass # Atom_names Residue Element_symbol(CAP) '
              f'RECORD ff_type{bcolors.ENDC}')

    def __get_atoms_info(self,
                         Masses_df: pd.DataFrame,  # df contains atoms' info
                         ) -> pd.DataFrame:
        """get data in the Masses_df and check them"""
        Masses: pd.DataFrame = Masses_df
        return Masses

    def __mk_pdb_df(self) -> pd.DataFrame:
        """make pdb dataframe"""
        pdb_df: pd.DataFrame  # The main df for pdb file
        columns: list[str]  # names of the columns of the pdb
        columns = ['records',
                   'atom_id',  # integer
                   'atom_name',  # left character
                   'l_indicator',  # character
                   'residue_name',  # right character
                   'chain_id',  # character
                   'residue_id',  # right integer
                   'Code_residues',  # character
                   'x',  # orthogonal Å coordinate right real (8.3)
                   'y',  # orthogonal Å coordinate right real (8.3)
                   'z',  # orthogonal Å coordinate right real (8.3)
                   'occupancy',  # right real (6.2)
                   'temperature',  # right real (6.2)
                   'Segment_id',  # left character
                   'element',  # right character
                   'charge',  # character
                   'ff_type',  # Type of the atom in Force field: opls_XX
                   'mass',  # Masses of the atoms
                   'q'  # Value of the charge!
                   ]
        pdb_df = pd.DataFrame(columns=columns)
        return pdb_df

    def __set_pdb(self,
                  pdb_df: pd.DataFrame,  # Empty df with columns name
                  Masses: pd.DataFrame,  # Checked Masses section
                  Atoms_df: pd.DataFrame  # Atoms coordinates
                  ) -> pd.DataFrame:
        names: list[str] = []  # Name of the atoms from Masses
        elements: list[str] = []  # Symbole for each atom
        residues: list[str] = []  # Names of each residues
        records: list[str] = []  # Records of each atom, e.g., ATOM, HATOM etc
        ff_type: list[str] = []  # Type of the atom in the FF, e.g., opls_XXX
        atoms_masses: list[float] = []  # Masses of the atoms
        # set columns of the df
        try:
            for item in Atoms_df['typ']:
                df_row = Masses[Masses['typ'] == item]
                names.append(df_row['names'][item])
                elements.append(df_row['elements'][item])
                residues.append(df_row['residues'][item])
                records.append(df_row['records'][item])
                ff_type.append(df_row['ff_type'][item])
                atoms_masses.append(df_row['mass'][item])
        except KeyError:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}: '
                     f'({self.__module__})\n'
                     '\tERROR!\n'
                     '\tThe atoms in the masses section of the input file must\n'
                     '\thave the follwoing order:\n'
                     '\tid mass # Atom_names Residue Element_symbol(CAP) '
                     'RECORD ff_type\n'
                     f'{bcolors.ENDC}'
                     )
        names = self.__fix_atom_names(names,
                                      Atoms_df['mol'],
                                      Atoms_df['atom_id'])
        pdb_df['atom_name'] = names
        empty_data: list[str] = [' ' for _ in names]  # For empty
        pdb_df['element'] = elements
        pdb_df['residue_name'] = residues
        pdb_df['records'] = records
        pdb_df['chain_id'] = empty_data
        # atom_id = [int(item)+32 for item in Atoms_df['atom_id']]
        pdb_df['atom_id'] = Atoms_df['atom_id']
        pdb_df['x'] = Atoms_df['x']
        pdb_df['y'] = Atoms_df['y']
        pdb_df['z'] = Atoms_df['z']
        pdb_df['l_indicator'] = ['A' for _ in names]
        pdb_df['occupancy'] = [1.0 for _ in names]
        pdb_df['temperature'] = empty_data
        pdb_df['Segment_id'] = empty_data
        pdb_df['residue_id'] = Atoms_df['mol']
        pdb_df['Code_residues'] = empty_data
        pdb_df['charge'] = empty_data
        pdb_df['ff_type'] = ff_type
        pdb_df['mass'] = atoms_masses
        pdb_df['q'] = Atoms_df['charge']
        return pdb_df

    def __fix_atom_names(self,
                         names: list[str],  # Name of the atoms from LAMMPS
                         mol_id: list[int],  # Id of each mol
                         atom_id: list[int]  # Id of the atoms
                         ) -> list[str]:
        """make the names by adding index to each similar name"""
        # First seprate residues = having same mol index
        names_id_df: pd.DataFrame  # df of names and mol_id
        names_id_df = pd.DataFrame({'atom_id': atom_id,
                                    'names': names,
                                    'mol_id': mol_id
                                    })
        watch_id: list[int]  # uniqe mol id in the mol_id
        watch_id = list(set(mol_id))
        mol_list: list[pd.DataFrame] = []  # df with one mol id
        for mol in watch_id:
            df1 = pd.DataFrame(columns=['atom_id'
                                        'names'
                                        'mol_id',
                                        'id_name'])
            df_i: pd.DataFrame = names_id_df[names_id_df['mol_id'] == mol]
            id_name: list[str] = self.__rename_atoms(df_i['names'], mol)
            df1['atom_id'] = df_i['atom_id']
            df1['names'] = df_i['names']
            df1['mol_id'] = df_i['mol_id']
            df1['id_name'] = id_name
            mol_list.append(df1)
            del df_i, df1
        rename_df: pd.DataFrame  # df with updated names to orderd with atom id
        rename_df = pd.concat(mol_list)
        rename_df.drop(columns=['names'], inplace=True)
        rename_df.sort_values(by=['atom_id'], inplace=True)
        return rename_df['id_name']

    def __rename_atoms(self,
                       names: list[str],  # Name of the atoms from LAMMPS data
                       mol: int  # Molecule id
                       ) -> list[str]:
        """rename the atoms based on thier repetetion"""
        # Get the repeated item by counter and rename them
        name_id_dict: dict[str, typing.Any]  # Number of repeatation for each
        name_id_dict = {a: list(range(1, b+1)) if b > 1 else ''
                        for a, b in Counter(names).items()}
        name_id: list[str] = [f'{i}{name_id_dict[i].pop(0)}' if
                              len(name_id_dict[i]) else i for i in names]
        for i, name in enumerate(name_id):
            if len(name) > 4:
                print(f'{bcolors.WARNING}\tWarning:\n'
                      f'\t\tLenght of item {i}: {name} '
                      f'is longer than 4, consider renaming the atoms\n')
        if len(names) != len(set(name_id)):
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n\tERROR! '
                     f'There is similar name in a same molecule! Nr.: {mol}'
                     f'{bcolors.ENDC}\n')
        return name_id

    def print_info(self) -> None:
        """Just to subpress the pylint error"""
