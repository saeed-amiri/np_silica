import sys
import itertools
import pandas as pd
import read_ff_charm as charm
import itp_to_df as itp
import pdb_to_df as pdb
import static_info as constants

class Doc:
    """write data for the spcial case of the nano particles"""

class WRITE_DATA :
    """
    write out the output file for lammps
    """
    def __init__(self,
                 pdb: pdb,  # information from pdb file
                 itp: itp,  # information from itp file
                 fout: str   # name of the output file
                 ) -> None:
        self.pdb = pdb
        self.itp = itp
        self.write_file(fout)
        cnt = constants()
        
        del pdb, itp

    def write_file(self,
                   fout: str  # Name of the output file
                   ) -> None:
        with open(fout, 'w') as f:
            f.write(f'#\twrite data file with "{sys.argv[0]}"')
            f.write(f'\n')
            self.write_numbers_types(f)
            self.write_box(f)
            self.write_masses(f)
            # self.write_charges(f)
            self.write_pair_coef(f)
            self.write_bonds_coef(f)
            self.write_angles_coef(f)
            self.write_coords(f)
            self.write_velocity(f)
            self.write_bonds(f)
            self.write_angles(f)

    
    def write_numbers_types(self, f) -> None:
        f.write(f'{self.pdb.NAtoms}\tatoms\n')
        f.write(f'{len(ATOM_TYPE)}\tatom types\n')
        f.write(f'{self.itp.NmBonds}\tbonds\n')
        f.write(f'{itp.TypBonds}\tbond types\n')
        f.write(f'{self.itp.NmAngles}\tangles\n')        
        f.write(f'{itp.TypAngles}\tangle types\n')
        f.write(f'\n')

    def write_box(self, f):
        f.write(f'{self.pdb.xlo}\t{self.pdb.xhi}\txlo\txhi\n')
        f.write(f'{self.pdb.ylo}\t{self.pdb.yhi}\tylo\tyhi\n')
        f.write(f'{self.pdb.zlo}\t{self.pdb.zhi}\tzlo\tzhi\n')
        f.write(f'\n')

    def write_masses(self, f) -> None:
        f.write(f'Masses\n\n')
        for i, (m, t) in enumerate(zip(ATOM_MASS, ATOM_TYPE)):
            if m != t : exit(f"WRONG mass and type sets: {t} :: {m}")
            f.write(f'{ATOM_TYPE[m]}\t{ATOM_MASS[t]}\t# {m}\n')
        f.write(f'\n')    

    def write_charges(self, f):
        """
        added to ATOM card
        """
        f.write(f'# set charges\n\n')
        for t in ATOM_TYPE:
            f.write(f'set type {ATOM_TYPE[t]} charge {ATOM_CHARGE[t]} # {t}\n')
        f.write(f'\n')

    def write_pair_coef(self, f) -> None:
        pairs = self.make_pairs()
        f.write(f'PairIJ Coeffs\n\n')
        for i, item in enumerate(list(pairs)):
            ai = [name for name, id in ATOM_TYPE.items() if id == item[0]][0]
            aj = [name for name, id in ATOM_TYPE.items() if id == item[1]][0]
            _pairIJ = f'{ai}_{aj}'
            try:
                sigma = self.__return_df_value(charmm.pair_df, 'pairs', _pairIJ, 'sigma' )/100
                epsilon = self.__return_df_value(charmm.pair_df, 'pairs', _pairIJ, 'epsilon' )*10
            except:
                sigma = 0.0
                epsilon = 1.0
            if _pairIJ in self.itp.SetBonds: exit('EXIT!! there is a pair with bonding and non-bonding interactions!!')
            f.write(f'{item[0]}\t{item[1]}  {epsilon} {sigma} {14} # {_pairIJ} \n')
        f.write('\n')

    def write_bonds_coef(self, f) -> None:
        f.write('Bond Coeffs # harmonic\n\n')
        for i, item in enumerate(self.itp.SetBonds):
            Kb = self.__return_df_value(charmm.bond_df, 'bond', item, 'Kb')/100
            b0 = self.__return_df_value(charmm.bond_df, 'bond', item, 'b0')*10
            f.write(f'{i+1} {Kb} {b0} # {item}\n')
        f.write('\n')

    def write_angles_coef(self, f) -> None:
        f.write(f'Angle Coeffs\n\n')
        for i, item in enumerate(self.itp.SetAngles):
            try:
                theta = self.__return_df_value(charmm.angle_df, 'angle', item, 'th0')
                K0 = self.__return_df_value(charmm.angle_df, 'angle', item, 'cth')
            except:
                theta = 0.00
                K0 = 0.00
            f.write(f'{i+1}\t{K0}\t{theta} # {item}\n')
        f.write(f'\n')

    def write_coords(self, f) -> None:
        f.write(f'Atoms\t#\tfull\n\n')
        self.pdb.df.to_csv(f, mode='a',header=None, index=False, sep=' ')
        f.write(f'\n')

    def write_velocity(self, f) -> None:
        f.write(f'Velocities\n\n')
        make_velocity_df(pdb.NAtoms).to_csv(f, mode='a', index=True, header=None, sep=' ')
        f.write(f'\n')

    def write_bonds(self, f) -> None:
        f.write(f'Bonds\n\n')
        self.itp.itpBondsDf.to_csv(f, mode='a', header=None, index=False, sep=' ')
        f.write(f'\n')

    def write_angles(self, f) -> None:
        f.write(f'Angles\n\n')
        self.itp.itpAnglesDf.to_csv(f, mode='a', header=None, index=False, sep=' ')

    def make_pairs(self) -> None:
        _type_list = [i+1 for i, _ in enumerate(ATOM_TYPE) ]
        return itertools.combinations_with_replacement(_type_list, 2)

    def __return_df_value(df, name, name_i,  string_re) -> float:
        # df: dataframe
        # name: name of the column to look for a match
        # name_i: the string to look for in the 'name' column
        # string_re: return the value of the selected coulmn
        return df.loc[df[name]==name_i][string_re].values[0]