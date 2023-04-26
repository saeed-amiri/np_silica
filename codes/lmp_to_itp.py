"""read the LAMMPS data file and get force field parameters and write
        `itp` formats for GROMACS conversion


    informations in the itp file:
    [ moleculetype ] : defines the name of your molecule in this top
    and nrexcl = 3 stands for excluding non-bonded interactions between
    atoms that are no further than 3 bonds away.

    [ atoms ] : defines the molecule, where nr and type are fixed,
    the rest is user defined. So atom can be named as you like,
    cgnr made larger or smaller (if possible, the total charge of
    a charge group should be zero), and charges can be changed here
    too.

    [ bonds ] : no comment.

    [ pairs ] : LJ and Coulomb 1-4 interactions

    [ angles ] : no comment

    [ dihedrals ] : in this case there are 9 proper dihedrals
    (funct = 1), 3 improper (funct = 4) and no Ryckaert-Bellemans type
    dihedrals.
"""


import typing
import pandas as pd
import read_lmp_data as relmp
import static_info as stinfo
from colors_text import TextColor as bcolors

# Check the residues names to see if there is bo/an/di between
# different residues, which is needed here;324
CHECK_RES: bool = True


class Itp:
    """get data from main"""
    def __init__(self,
                 lmp: relmp.ReadData,  # LAMMPS data file
                 pdb_df: pd.DataFrame  # Final df for pdb file
                 ) -> None:
        self.mk_itp(lmp, pdb_df)
        self.print_info()

    def mk_itp(self,
               lmp: relmp.ReadData,  # LAMMPS data file
               pdb_df: pd.DataFrame  # Final df for pdb file
               ) -> None:
        """call functions"""
        self.atoms = self.__mk_atoms(pdb_df)
        self.bonds = self.__mk_bonds(lmp)
        self.angles = self.__mk_angles(lmp)
        self.dihedrals = self.__mk_dihedrals(lmp)

    def __mk_atoms(self,
                   pdb_df: pd.DataFrame  # Final df for pdb file
                   ) -> pd.DataFrame:
        df_i: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['atomnr',  # Index of the atoms in lmp file
                   'atomtype',  # FF type of the atom, e.g. for c, opls_135
                   'resnr',  # Index of the residue which atom is belonged
                   'resname',  # Name of the residue which atom is belonged
                   'atomname',  # Name of the atom as in PDB file
                   'chargegrp',  # Charge group as in PDB file
                   'charge',  # Charge as in the lmp file
                   'mass',  # Mass odf the atom as in the lmp file
                   ' ',  # Comment column for ;
                   'element'  # Second column for the coments
                   ]
        df_i = pd.DataFrame(columns=columns)
        df_i['atomnr'] = [int(item) for item in pdb_df['atom_id']]
        df_i['atomtype'] = [str(item) for item in pdb_df['ff_type']]
        df_i['resnr'] = [int(item) for item in pdb_df['residue_id']]
        df_i['resname'] = [str(item) for item in pdb_df['residue_name']]
        df_i['atomname'] = [str(item) for item in pdb_df['atom_name']]
        df_i['chargegrp'] = [int(1) for _ in pdb_df['charge']]
        df_i['charge'] = [float(item) for item in pdb_df['q']]
        df_i['mass'] = [float(item) for item in pdb_df['mass']]
        df_i[' '] = [';' for _ in df_i['atomnr']]
        df_i['element'] = [str(item) for item in pdb_df["element"]]
        return df_i

    def __mk_bonds(self,
                   lmp: relmp.ReadData  # LAMMPS data file
                   ) -> pd.DataFrame:
        df_i: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['ai',  # 1st atom in bond
                   'aj',  # 2nd atom in bond
                   'funct',  # not sure what is this, just set to 1, or empty!
                   'r',  # Distance parameter in harmonic bond interactions
                   'k',  # Harmonic constant in the harmonic interactions
                   ' ',  # Comment: ;
                   '  ',  # Comment: name of the bond
                   'resname',  # Name of the residue which atoms belonged to
                   'resnr',  # Nr. of the residue which atoms belonged to
                   ]
        df_i = pd.DataFrame(columns=columns)
        Bonds_df: pd.DataFrame = lmp.Bonds_df.copy().sort_values(by='ai')
        df_i['ai'] = Bonds_df['ai']
        df_i['aj'] = Bonds_df['aj']
        # df_i['funct'] = [1 for _ in df_i['ai']]
        try:
            df_i[' '] = [';' for _ in df_i['ai']]
            df_i['name'] = lmp.Bonds_df['name']
        except KeyError:
            df_i['name'] = self.__mk_boandi_name(df_i,
                                                 ['ai', 'aj'],
                                                 lmp.Atoms_df)
            print(f'{bcolors.WARNING}{self.__class__.__name__}: '
                  f'({self.__module__})\n'
                  f'\tThere is no bonds` names in LAMMPS read data'
                  f'{bcolors.ENDC}')
        if stinfo.BoAnDi.BONDS_FLAG:
            df_i = self.__get_boandi_para(df_i,
                                          stinfo.BoAnDi.BONDS,
                                          ['r', 'k', 'funct'])
        if CHECK_RES:
            df_i['resname'], df_i['resnr'] = self.__get_bonds_res(lmp, df_i)
        return df_i

    def __get_bonds_res(self,
                        lmp: relmp.ReadData,  # LAMMPS data file
                        df_i: pd.DataFrame  # df contain itp info
                        ) -> tuple[list, list]:
        """return residues name and index"""
        resname: list[str] = []  # Name of the residues
        resnr: list[int] = []  # index of the residues
        flag_war: bool = False  # To print the warning
        for a_i, a_j in zip(df_i['ai'], df_i['aj']):
            ai_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_i]['typ'][a_i]
            aj_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_j]['typ'][a_j]
            mol_i: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ai_type]['residues'][ai_type]
            mol_j: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == aj_type]['residues'][aj_type]
            mol_iid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_i]['mol'][a_i]
            mol_jid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_j]['mol'][a_j]
            if mol_i != mol_j or mol_iid != mol_jid:
                if not flag_war:
                    print(f'{bcolors.WARNING}{self.__class__.__name__}: '
                          f'({self.__module__})\n'
                          '\tBond between atoms with different residues '
                          f'types{bcolors.ENDC}')
                    flag_war = True
            resnr.append(mol_iid)
            resname.append(mol_i)
        return resname, resnr

    def __mk_angles(self,
                    lmp: relmp.ReadData  # LAMMPS data file
                    ) -> pd.DataFrame:
        df_i: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['ai',  # 1st atom in angle
                   'aj',  # 2nd atom in angle
                   'ak',  # 3rd atom in angle
                   'funct',  # not sure what is this, just set to 1, or empty!
                   'theta',  # The angle between bonds
                   'cth',  # Strength of the bonds
                   ' ',  # Comment: name of the angle
                   'angle_name',  # Comment: name of the angle
                   'resname',  # Name of the residue which atoms belonged to
                   'resnr',  # Nr. of the residue which atoms belonged to
                   ]
        df_i = pd.DataFrame(columns=columns)
        Angles_df: pd.DataFrame = lmp.Angles_df.sort_values(by='ai')
        df_i['ai'] = Angles_df['ai']
        df_i['aj'] = Angles_df['aj']
        df_i['ak'] = Angles_df['ak']
        df_i['funct'] = [1 for _ in df_i['ai']]
        df_i[' '] = [';' for _ in df_i['ai']]
        try:
            df_i['name'] = lmp.Angles_df['name']
        except KeyError:
            df_i['name'] = self.__mk_boandi_name(df_i,
                                                 ['ai', 'aj', 'ak'],
                                                 lmp.Atoms_df)
            print(f'{bcolors.WARNING}{self.__class__.__name__}: '
                  f'({self.__module__})\n'
                  f'\tGetting names for the angles ...'
                  f'{bcolors.ENDC}')
        if stinfo.BoAnDi.ANGLES_FLAG:
            df_i = self.__get_boandi_para(df_i,
                                          stinfo.BoAnDi.ANGLES,
                                          ['theta', 'cth', 'funct'])
        if CHECK_RES:
            df_i['resname'], df_i['resnr'] = self.__get_angles_res(lmp, df_i)
        return df_i

    def __get_angles_res(self,
                         lmp: relmp.ReadData,  # LAMMPS data file
                         df_i: pd.DataFrame  # df contain itp info
                         ) -> tuple[list, list]:
        """return residues name and index"""
        resname: list[str] = []  # Name of the residues
        resnr: list[int] = []  # index of the residues
        flag_war: bool = False  # To print the warning
        for a_i, a_j, a_k in zip(df_i['ai'], df_i['aj'], df_i['ak']):
            ai_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_i]['typ'][a_i]
            aj_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_j]['typ'][a_j]
            ak_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_k]['typ'][a_k]
            mol_i: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ai_type]['residues'][ai_type]
            mol_j: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == aj_type]['residues'][aj_type]
            mol_k: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ak_type]['residues'][ak_type]
            mol_iid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_i]['mol'][a_i]
            mol_jid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_j]['mol'][a_j]
            mol_kid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_k]['mol'][a_k]
            check_list_name = set([mol_i, mol_j, mol_k])
            check_list_id = set([mol_iid, mol_jid, mol_kid])
            if len(check_list_name) != 1 or len(check_list_id) != 1:
                if not flag_war:
                    print(f'{bcolors.WARNING}{self.__class__.__name__}: '
                          f'({self.__module__})\n'
                          f'\tAngles between atoms with different residues '
                          f'types{bcolors.ENDC}')
                flag_war = True
            resnr.append(mol_iid)
            resname.append(mol_i)
        return resname, resnr

    def __mk_dihedrals(self,
                       lmp: relmp.ReadData  # LAMMPS data file
                       ) -> pd.DataFrame:
        df_i: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['ai',  # 1st atom in dihedrals
                   'aj',  # 2nd atom in dihedrals
                   'ak',  # 3rd atom in dihedrals
                   'ah',  # 4th atom in dihedrals
                   'funct',  # not sure what is this, just set to 3, or empty!
                   'C0',  # Dihedrals parameters
                   'C1',  # Dihedrals parameters
                   'C2',  # Dihedrals parameters
                   'C3',  # Dihedrals parameters
                   'C4',  # Dihedrals parameters
                   'C5',  # Dihedrals parameters
                   ' ',  # Comment: name of the dihedrals
                   'dihedral_name',  # names
                   'resname',  # Name of the residue which atoms belonged to
                   'resnr',  # Nr. of the residue which atoms belonged to
                   ]
        df_i = pd.DataFrame(columns=columns)
        try:
            Dihedrals_df: pd.DataFrame = lmp.Dihedrals_df.sort_values(by='ai')
            df_i['ai'] = Dihedrals_df['ai']
            df_i['aj'] = Dihedrals_df['aj']
            df_i['ak'] = Dihedrals_df['ak']
            df_i['ah'] = Dihedrals_df['ah']
            df_i['funct'] = [3 for _ in df_i['ai']]
            df_i[' '] = [';' for _ in df_i['ai']]
            try:
                df_i['name'] = lmp.Dihedrals_df['name']
            except KeyError:
                df_i['name'] = self.__mk_boandi_name(df_i,
                                                     ['ai', 'aj', 'ak', 'ah'],
                                                     lmp.Atoms_df)
                print(f'{bcolors.WARNING}{self.__class__.__name__}: '
                      f'({self.__module__})\n'
                      '\tGetting names for the dihedrals ...'
                      f'{bcolors.ENDC}')
            if stinfo.BoAnDi.DIHEDRALS_FLAG:
                df_i = self.__get_boandi_para(df_i,
                                              stinfo.BoAnDi.DIHEDRLAS,
                                              ['C0', 'C1', 'C2', 'C3',
                                               'C4', 'C5', 'funct'])
            if CHECK_RES:
                df_i['resname'], df_i['resnr'] =\
                    self.__get_dihedrals_res(lmp, df_i)
        except KeyError:
            pass
        return df_i

    def __get_dihedrals_res(self,
                            lmp: relmp.ReadData,  # LAMMPS data file
                            df_i: pd.DataFrame  # df contain itp info
                            ) -> tuple[list, list]:
        """return residues name and index"""
        resname: list[str] = []  # Name of the residues
        resnr: list[int] = []  # index of the residues
        flag_war: bool = False  # To print the warning
        for a_i, a_j, a_k, a_h in zip(df_i['ai'],
                                      df_i['aj'],
                                      df_i['ak'],
                                      df_i['ah']):
            ai_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_i]['typ'][a_i]
            aj_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_j]['typ'][a_j]
            ak_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_k]['typ'][a_k]
            ah_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == a_h]['typ'][a_h]
            mol_i: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ai_type]['residues'][ai_type]
            mol_j: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == aj_type]['residues'][aj_type]
            mol_k: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ak_type]['residues'][ak_type]
            mol_h: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ah_type]['residues'][ah_type]
            mol_iid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_i]['mol'][a_i]
            mol_jid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_j]['mol'][a_j]
            mol_kid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_k]['mol'][a_k]
            mol_hid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == a_h]['mol'][a_h]
            check_list_name = set([mol_i, mol_j, mol_k, mol_h])
            check_list_id = set([mol_iid, mol_jid, mol_kid, mol_hid])
            if len(check_list_name) != 1 or len(check_list_id) != 1:
                if not flag_war:
                    print(f'{bcolors.WARNING}{self.__class__.__name__}: '
                          f'({self.__module__})\n'
                          f'\tDihedrals between atoms with different residues '
                          f'types{bcolors.ENDC}')
                    flag_war = True

            resnr.append(mol_iid)
            resname.append(mol_i)
        return resname, resnr

    def __get_boandi_para(self,
                          df_i: pd.DataFrame,  # To get parameter for ITP
                          data_dict: dict[str, typing.Any],  # Data from stinfo
                          columns: list[str]  # Name of the columns to replace
                          ) -> pd.DataFrame:
        """try to get parameters for the bonds, angles, and dihedrals
        from static_info module"""
        df_c: pd.DataFrame = df_i.copy()
        if not stinfo.Hydration.FFIELD == 'charmm':
            for item, row in df_i.iterrows():
                try:
                    for col in columns:
                        df_c.at[item, col] = data_dict.get(row['name'])[col]
                except (TypeError, KeyError):
                    print(f'{bcolors.CAUTION}{self.__class__.__name__}: '
                          f'({self.__module__})\n'
                          f'\tNo or bad info for the `{col}` with name: '
                          f'{row["name"]}{bcolors.ENDC}')
        else:
            for item, row in df_i.iterrows():
                for col in columns:
                    if col == 'funct':
                        df_c.at[item, col] = data_dict.get(row['name'])[col]
                    else:
                        df_c.at[item, col] = ' '
        return df_c

    def __mk_boandi_name(self,
                         df_i: pd.DataFrame,  # The dataframe
                         a_list: list[str],  # Atoms involved, e.g. ['ai','aj',
                         atoms_df: pd.DataFrame  # All atom info
                         ) -> list[str]:
        """make a name column for the bonds"""
        atom_name: dict[int, str]  # id and name of the atoms
        atom_name = dict(zip(atoms_df['atom_id'], atoms_df['name']))
        name_list: list[str] = []  # Name of the bo/an/di
        for _, row in df_i.iterrows():
            names = []
            for atom in a_list:
                names.append(atom_name[row[atom]])
            name_list.append('_'.join(names))
        return name_list

    def print_info(self) -> None:
        """Just to subpress the pylint error"""
