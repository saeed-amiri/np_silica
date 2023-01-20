import re
import sys
import pandas as pd

class Doc:
    """read Charm Force fields"""


class CHARMM:
    """
    reading charmm36_silica.itp
    to get interactions parameters
    """

    def __init__(self,
                 fcharm: str  # Charmm file
                 ) -> None:
        self.__read_charmm

    def __read_charmm(self,
                    fcharm: str  # Charmm file
                    ) -> pd.DataFrame:
        # Define flages
        atomtypes, nonbond_params, bondtypes, pairtypes, \
        angletypes, dihedraltypes = False, False, False, False, False, False
        # Initiate lists
        atomList, nonbondList, bondList, pairList,  angleList, \
        dihedralList = [], [], [], [], [], []
        with open(fcharm, 'r') as f:
            while True:
                line = f.readline()
                if line.startswith(';'):
                    continue
                if line.startswith('[ atomtypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes,\
                    angletypes, dihedraltypes = True, False, False, False,\
                    False, False
                if line.startswith('[ nonbond_params ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes,\
                    angletypes, dihedraltypes = False, True, False, False,\
                    False, False
                if line.startswith('[ bondtypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes,\
                    angletypes, dihedraltypes = False, False, True, False,\
                    False, False
                if line.startswith('[ pairtypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes,
                    angletypes, dihedraltypes = False, False, False, True,\
                    False, False
                if line.startswith('[ angletypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes,
                    angletypes, dihedraltypes = False, False, False, False,\
                    True, False
                if line.startswith('[ dihedraltypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes,
                    angletypes, dihedraltypes = False, False, False, False,\
                    False, True
                if line.strip() and atomtypes and not line.startswith('[ '):
                    line = self.__drop_semicolon(line).strip()
                    atomList.append(line)
                if line.strip() and nonbond_params and \
                        not line.startswith('[ '):
                    nonbondList.append(line)
                if line.strip() and bondtypes and not line.startswith('[ '):
                    bondList.append(line)
                if line.strip() and pairtypes and not line.startswith('[ '):
                    pairList.append(line)
                if line.strip() and angletypes and not line.startswith('[ '):
                    angleList.append(line)
                if line.strip() and dihedraltypes and \
                    not line.startswith('[ '):
                    dihedralList.append(line)
                if not line: break
        self.atomtyps_df = self.__read_atomtypes(atomList); del atomList
        self.nonbond_df = self.__read_nonbond(nonbondList); del nonbondList
        self.bond_df = self.__read_bondtypes(bondList); del bondList
        self.pair_df = self.__read_pairtypes(pairList); del pairList
        self.angle_df = self.__read_angletypes(angleList); del angleList


    def __read_atomtypes(self, atomList) -> pd.DataFrame:
        atom_dict = dict(name=[],
                        at_num=[],
                        mass=[],
                        charge=[],
                        ptype=[],
                        sigma=[],
                        epsilon=[])
        for item in atomList:
            i_name, i_at_num, i_mass, i_charge, i_ptype, i_sigma,\
                i_epsilon = self.__procces_lines(item,7)
            i_name = i_name.strip()
            i_at_num = int(i_at_num.strip())
            i_mass = float(i_mass.strip())
            i_charge = float(i_charge.strip())
            i_ptype = i_ptype.strip()
            i_sigma = float(i_sigma.strip())
            i_epsilon = float(i_epsilon.strip())
            atom_dict['name'].append(i_name)
            atom_dict['at_num'].append(i_at_num)
            atom_dict['mass'].append(i_mass)
            atom_dict['charge'].append(i_charge)
            atom_dict['ptype'].append(i_ptype)
            atom_dict['sigma'].append(i_sigma)
            atom_dict['epsilon'].append(i_epsilon)
        del atomList
        return pd.DataFrame.from_dict(atom_dict)    
    
    def __read_nonbond(self, nonbondList) -> pd.DataFrame:
        nonbond_dict = dict(ai=[], aj=[], func=[], sigma=[], epsilon=[])
        for item in nonbondList:
            i_ai, i_aj, i_func, i_sigma, i_epsilon\
                = self.__procces_lines(item, 5)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_func = i_func.strip()
            i_sigma = float(i_sigma.strip())
            i_epsilon = float(i_epsilon.strip())
            nonbond_dict['ai'].append(i_ai)
            nonbond_dict['aj'].append(i_aj)
            nonbond_dict['func'].append(i_func)
            nonbond_dict['sigma'].append(i_sigma)
            nonbond_dict['epsilon'].append(i_epsilon)
        del nonbondList
        return pd.DataFrame.from_dict(nonbond_dict)


        
    def __read_bondtypes(self, bondList) -> pd.DataFrame:
        bond_dict = dict(bond=[], func=[], b0=[], Kb=[])
        for item in bondList:
            i_ai, i_aj, i_func, i_b0, i_Kb = self.__procces_lines(item, 5)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_func = i_func.strip()
            i_b0 = float(i_b0.strip())
            i_Kb = float(i_Kb.strip())
            bond_dict['bond'].append(f'{i_ai}_{i_aj}')
            bond_dict['func'].append(i_func)
            bond_dict['b0'].append(i_b0)
            bond_dict['Kb'].append(i_Kb)
        del bondList
        return pd.DataFrame.from_dict(bond_dict)

    def __read_pairtypes(self, pairList) -> pd.DataFrame:
        pair_dict = dict(pairs=[], func=[], sigma=[], epsilon=[])
        for item in pairList:
            i_ai, i_aj, i_func, i_sigma, i_epsilon \
                = self.__procces_lines(item, 5)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_func = i_func.strip()
            i_sigma = float(i_sigma.strip())
            i_epsilon = float(i_epsilon.strip())
            i_pairs = f'{i_ai}_{i_aj}'
            pair_dict['pairs'].append(i_pairs)
            pair_dict['func'].append(i_func)
            pair_dict['sigma'].append(i_sigma)
            pair_dict['epsilon'].append(i_epsilon)
        del pairList
        return pd.DataFrame.from_dict(pair_dict)

    def __read_angletypes(self, angleList) -> pd.DataFrame:
        angle_dict = dict(angle=[], func=[], th0=[], cth=[], S0=[], Kub=[])
        for item in angleList:
            i_ai, i_aj, i_ak, i_func, i_th0, i_cth, i_S0, i_Kub = \
                self.__procces_lines(item, 8)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_ak = i_ak.strip()
            i_func = i_func.strip()
            i_th0 = float(i_th0.strip())
            i_cth = float(i_cth.strip())
            i_S0 = float(i_S0.strip())
            i_Kub = float(i_Kub.strip())
            i_angle = f'{i_ai}_{i_aj}_{i_ak}'
            angle_dict['angle'].append(i_angle)
            angle_dict['func'].append(i_func)
            angle_dict['th0'].append(i_th0)
            angle_dict['cth'].append(i_cth)
            angle_dict['S0'].append(i_S0)
            angle_dict['Kub'].append(i_Kub)
        del angleList
        return pd.DataFrame.from_dict(angle_dict)

    def __read_dihedraltypes(self):
        pass

    def __drop_semicolon(self,
                       line: str # line to replace chrs
                       ) -> list:
        return re.sub(r'\;.*', "", line)

    def __procces_lines(self,
                      line: str,  # Char to clear 
                      lineLen: int  # To check the length
                      ) -> list:
        line =line.strip()
        line = line.split(' ')
        line = [item for item in line if item]
        if len(line) != lineLen :
            exit(f'WRONG LINE in line: {line}, EXIT!')
        return line