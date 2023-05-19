"""All the constant values use in the silanization"""

import typing
import os

SOURCE_DIR: str  # SOURCE OF THE DATA FILES
SOURCE_DIR = '/scratch/saeed/MyScripts/np_silica/data'


class AtomMass:
    """atmoic masses of the particles"""
    HO: float = 1.0080
    OB: float = 15.9994
    OH: float = 15.9994
    OM: float = 15.9994
    OMH: float = 15.9994
    OD: float = 15.9994
    SI: float = 28.0860
    SU: float = 28.0860
    SD: float = 28.0860


class AtomCharge:
    """atomic masses of the particles"""
    # Charges of the silica atoms, before silanization
    N: float = 0.000
    HO: float = 0.400
    OH: float = -0.800
    SI: float = 1.600
    OB: float = -0.800
    SD: float = 1.500
    OD: float = -1.000
    OM: float = -0.900
    OMH: float = -0.900
    SU: float = 1.200


class AtomType:
    """atomic types of the particles"""
    HO: int = 1
    OB: int = 2
    OH: int = 3
    OM: int = 4
    OMH: int = 5
    OD: int = 6
    SI: int = 7
    SU: int = 8
    SD: int = 9


class UpdateCharge:
    """charges of the Si-OM groups which lost O or O&H atoms
    If No need to change set them equal to None, ex.:
        OM: float = None
    """
    OM: float = -0.80
    SI: float = 1.600


class Constants:
    """The constants which are used in the script"""
    # The desire coverage for grafting on NP
    Coverage: float = 3.5
    # The thickness of the shell from surface to look for Si atoms
    Shell_radius: float = 6.0
    # calculate the level ups for Aminopropyl
    OM_n: int = 4  # Number of extra atoms (Si, OM) in aminopropyl
    Amino_OM: int = 3  # The default numbers of OM in the amino file
    # Silicon name in Aminopropyl file:
    SI_amino: str = 'SI'
    # OM name in Aminopropyl file:
    OM_amino: str = 'OM'
    # total number of atoms in aminoproyl file
    Num_amino: int = 17
    # Length of the ODAP is around 10.9372292106925 but I set to 11
    ODA_length: float = 11


class AtomGroup:
    """list of name of the atoms to work with"""
    # Si groups to find them in shell and add APTES to them
    SiGroup: list[str] = ['SD', 'SI']
    # Oxygen groups which are bonded to Si on the shell, SHOULD NOT replace!
    OMGroup: list[str] = ['OM', 'OB', 'OMH']
    # Oxygen groups bonded to Si on the shell to drop
    OxGroup: list[str] = ['OH']
    # Hydrogen groups bonded to the Ox groups to drop
    HyGroup: list[str] = ['HO', 'OMH', 'HR']
    # if wanted the exact number of Ox group to be grafted:
    # Sparse the si selected on the surface
    # The optines are: 
    # 'exact' -> must give a integer, otherwise it well stoped
    # 'random' -> sparse the si randomly
    SPARSE_STY: str = 'random' 
    EXACT_NUM: int = 17  # If SPARSE_STY is `exact`


class DataFile:
    """Get data directory which are used in the script"""
    # Prtonated data of the APTES
    APTES: str = os.path.join(SOURCE_DIR, 'aminopropyl_pro.data')
    # Unprtonated data of the APTES
    APTUN: str = os.path.join(SOURCE_DIR, 'aminopropyl_unpro.data')
    SI_DF: str = 'SI_DF'  # File with selected info of si atom in adding APTES
    SI_XYZ: str = 'SI_XYZ'  # File with  info of si atom in adding APTES


class Hydration:
    """set all the info for water box
    Limitation for the box are added to the maximum radius of the NP"""
    # Forcefield type:
    FFIELD: str = 'charmm'
    TOLERANCE: float = 2.0
    # Check and delete files if they are exsit, True -> check and delete
    CHECK_WATER_PDB: bool = False
    # Contact angle, it defeins how much of the nanoparticle should be
    # in the oil phase, in case there is oil phase the APTES on the oil
    # phase are unprotonated
    CONATCT_ANGLE: float = 90  # In degree; If negetive -> no oil, MAX depends!
    # Box dimensions
    # x
    X_MIN: float = -75.0
    X_MAX: float = 75.0
    # y
    Y_MIN: float = -75.0
    Y_MAX: float = 75.0
    # z
    Z_MIN: float = -75.0
    Z_MAX: float = 75.0
    # Constants
    WATER_DENSITY = 0.9998395  # g/ml
    WATER_MOLAR_MASS: float = 18.01528  # g/mol
    OIL_DENSITY = 0.730  # g/ml for Decane (D10)
    OIL_MOLAR_MASS: float = 142.286  # g/mol for Decane (D10)
    AVOGADRO: float = 6.0221408e+23  # 1/mol
    MASSES: dict[str, float] = {'HW': 1.0080,
                                'OW': 15.9994,
                                'OH': 15.9994,
                                'OR': 15.9994,
                                'CL': 35.453,
                                'NA': 22.989769,
                                'H': 1.00080,
                                'C': 12.011,
                                'N': 14.0067,
                                'POT': 39.0983,
                                'CLA': 35.453}
    # Water charges
    # Type of each water model
    WATER_CHARGE: dict[str, dict[str, float]]
    NA_Q: int = +1  # Charge of the na ion
    CL_Q: int = -1  # Charge of the cl ion
    WATER_MODEL: str = 'charmm'  # The model to use in system
    WATER_CHARGE = {'tip3p': {'OW': -0.834,
                              'HW1': 0.417,
                              'HW2': 0.417,
                              'NA': NA_Q,
                              'CL': CL_Q},
                    'spce': {'OW': -0.8476,
                             'HW1': 0.4238,
                             'HW2': 0.4238,
                             'NA': NA_Q,
                             'CL': CL_Q},
                    'spc': {'OW': -0.82,
                            'HW1': 0.41,
                            'HW2': 0.41,
                            'NA': NA_Q,
                            'CL': CL_Q},
                    'tip4p': {'OW': 0.0,
                              'HW1': 0.52,
                              'HW2': 0.52,
                              'MW': -1.04,
                              'NA': NA_Q,
                              'CL': CL_Q},
                    'charmm': {'OH2': 0.0,
                              'H1': 0.52,
                              'H2': 0.52,
                              'POT': NA_Q,
                              'CLA': CL_Q}
                    }
    # PACKMOL files
    
    WATER_PDB: str = os.path.join(SOURCE_DIR, 'water_charmm.pdb')
    ODAP_PDB: str = os.path.join(SOURCE_DIR, 'ODAp_charmm.pdb')
    ODAN_PDB: str = os.path.join(SOURCE_DIR, 'ODAn_charmm.pdb')
    NA_PDB: str = os.path.join(SOURCE_DIR, 'POT.pdb')
    CL_PDB: str = os.path.join(SOURCE_DIR, 'CLA.pdb')
    OIL_PDB: str = os.path.join(SOURCE_DIR, 'D10.pdb')
    ADD_ION: bool = False  # if True it will add the ion to the itp file
    NA_ITP: str = os.path.join(SOURCE_DIR, 'Na.itp')
    CL_ITP: str = os.path.join(SOURCE_DIR, 'Cl.itp')
    ODAP_ITP: str = os.path.join(SOURCE_DIR, 'ODAp_charmm.itp')
    ODAN_ITP: str = os.path.join(SOURCE_DIR, 'ODAn_charmm.itp')
    OIL_ITP: str = os.path.join(SOURCE_DIR, 'D10_charmm.itp')
    INP_FILE: str = 'water_box.inp'
    OUT_FILE: str = 'water_box.pdb'
    WS_INP: str = 'water_silica.inp'  # Input for final water & silanized file
    GRO_PDB: str = 'silica_water.pdb'  # File for GROMACS input
    # PACKMOL lib
    PACKMOL: str = '/home/saeed/Downloads/packmol/packmol'
    # Number or concentration of ODAP and ODAN (in case later wanted)
    # It is used in the write_water and lmp_itp_pdb
    # protonation of ODA in water:
    # If you want the ODA in water to
    # be protonated, set ODAP_PROTONATION to True, and the script will
    # use the data from ODAP. However, remember to set the number of
    # ODA in water still using N_ODAP. The script will select the file
    # based on the protonation state in ODAP_PROTONATION, and N_ODAN
    # should only be used if you want the ODA in the oil phase.
    ODAP_PROTONATION: bool = True
    N_ODAP: int = 50  # Protonated ODA will add to water section
    N_ODAN: int = 0 # Unprotonated ODA will add to if oil section
    # If the protonated ODA should be at the interface at the beginning,
    # set the following to INTERFACE or if they shoud be in the edge of
    # the oil phase set it to `OILDOWN` and in water edge put it to
    # `WATERTOP`:
    ODAP_INTERFACE: bool = 'WATERTOP'
    # Salt (NaCl) parameters
    # Need a tuple type of concentration or molality
    # For now it only support molality
    # Molal: Containing one mole of solute per kilogram of solvent.
    # Molal should be in `MILIMOLAL`:  millimoles per kg !!!
    N_NACL: dict[str, typing.Any]  # Type of concenteration and amount
    N_NACL = {'sty': 'mmolal', 'sum': 0}


class PosRes:
    """write the psition restrians for atoms in the core of the silica
    nanoparticels"""
    RESTRINS_GROUP: list[str] = ['COR', 'SIL']
    POSRES: bool = True  # if want to write it: True
    RES_FILE: str = 'STRONG_POSRES.itp'
    FUNCTION: int = 1  # Type of the function for the restrains
    FX: int = 1000  # Force along x axis
    FY: int = 1000  # Force along y axis
    FZ: int = 1000  # Force along z axis


class GroInp:
    """info for writing inputs for the gromacs input"""
    FORCEFIELD: str = '../../force_field/charmm36_silica.itp'
    WATERITP: str = '../../force_field/TIP3.itp'
    IONITP: str = '../../force_field/CLA.itp'
    TOPFILE: str = 'topol.top'
    NPPOSRES: bool = True  # True if want to set the restraints on NP
    WATERPOSRES: bool = False  #True of want to set restraints on NP


class PdbMass:
    """information for the masses section in the output file of the
    silanization"""
    silica_residue: str = 'SIL'  # Name of the residue for the silica NP
    core_residue: str = 'COR'  # Name of the residue to apply the restraints
    aptes_residue: str = 'APT'  # Name of the residue for the APTES
    odap_residue: str = 'ODA'  # Name of the protenated ODA
    odan_residue: str = 'ODN'  # Name of the unprotenated ODA
    water_residue: str = 'SOL'  # Name of the water residue
    oil_residue: str = 'D10'  # Name of the oil residue
    cl_residue: str = 'CL'  # Name of the cl residue
    na_residue: str = 'NA'  # Name of the na residue
    pot_residue: str = 'POT'  # Name of the pot residue
    cla_residue: str = 'CLA'  # Name of the cla residue
    HO: dict[str, typing.Any] = {'Atoms_names': 'HO',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'HO'
                                 }
    H: dict[str, typing.Any] = {'Atoms_names': 'H',
                                'Residue': odap_residue,
                                'Element_symbol': 'H',
                                'RECORD': 'ATOM',
                                'ff_type': 'opls_140'
                                }
    HR: dict[str, typing.Any] = {'Atoms_names': 'HR',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'HR'
                                 }
    OB: dict[str, typing.Any] = {'Atoms_names': 'OB',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'OB'
                                 }
    OH: dict[str, typing.Any] = {'Atoms_names': 'OH',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'OH'
                                 }
    OR: dict[str, typing.Any] = {'Atoms_names': 'OR',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'OR'
                                 }
    OM: dict[str, typing.Any] = {'Atoms_names': 'OM',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'OM'
                                 }
    OMH: dict[str, typing.Any] = {'Atoms_names': 'OMH',
                                  'Residue': silica_residue,
                                  'Element_symbol': 'O',
                                  'RECORD': 'ATOM',
                                  'ff_type': 'OMH'
                                  }
    OD: dict[str, typing.Any] = {'Atoms_names': 'OD',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'OD'
                                 }
    SI: dict[str, typing.Any] = {'Atoms_names': 'SI',
                                 'Residue':  silica_residue,
                                 'Element_symbol': 'SI',
                                 'RECORD': 'ATOM',
                                 'ff_type':  'SI'
                                 }
    SU: dict[str, typing.Any] = {'Atoms_names': 'SU',
                                 'Residue':  silica_residue,
                                 'Element_symbol': 'SI',
                                 'RECORD': 'ATOM',
                                 'ff_type':  'SU'
                                 }
    SD: dict[str, typing.Any] = {'Atoms_names': 'SD',
                                 'Residue':  silica_residue,
                                 'Element_symbol': 'SI',
                                 'RECORD': 'ATOM',
                                 'ff_type':  'SD'
                                 }
    NH2: dict[str, typing.Any] = {'Atoms_names': 'N',
                                'Residue': aptes_residue,
                                'Element_symbol': 'N',
                                'RECORD': 'ATOM',
                                'ff_type': 'NH2'
                                }
    NH3: dict[str, typing.Any] = {'Atoms_names': 'N',
                                'Residue': aptes_residue,
                                'Element_symbol': 'N',
                                'RECORD': 'ATOM',
                                'ff_type': 'NH3'
                                }
    F: dict[str, typing.Any] = {'Atoms_names': 'F',
                                'Residue': aptes_residue,
                                'Element_symbol': 'F',
                                'RECORD': 'ATOM',
                                'ff_type': 'opls_287'
                                }
    CH: dict[str, typing.Any] = {'Atoms_names': 'CH',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'C',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'CT2'
                                 }
    C: dict[str, typing.Any] = {'Atoms_names': 'C',
                                'Residue': odap_residue,
                                'Element_symbol': 'C',
                                'RECORD': 'ATOM',
                                'ff_type': 'opls_135'
                                }
    HC: dict[str, typing.Any] = {'Atoms_names': 'HC',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'HA2'
                                 }
    CN: dict[str, typing.Any] = {'Atoms_names': 'CT',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'C',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'CT2'
                                 }
    HN2: dict[str, typing.Any] = {'Atoms_names': 'HN',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'HC'
                                 }
    HN3: dict[str, typing.Any] = {'Atoms_names': 'HN',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'HC'
                                 }                                 
    ATOMS: dict[str, dict[str, typing.Any]] = {'HO': HO,
                                               'HR': HR,
                                               'OR': OR,
                                               'OB': OB,
                                               'OH': OH,
                                               'OM': OM,
                                               'OMH': OMH,
                                               'OD': OD,
                                               'SI': SI,
                                               'SU': SU,
                                               'SD': SD,
                                               'NH2': NH2,
                                               'NH3': NH3,
                                               'CH': CH,
                                               'HC': HC,
                                               'CN': CN,
                                               'HN2': HN2,
                                               'HN3': HN3,
                                               'C': C,
                                               'F': F
                                               }
    MOLECULES: dict[str, typing.Any]  # Contains all the molecules' atoms info


class BoAnDi:
    """interaction for bonds, angles, and dihedrals in silanized
    nanoparticles"""
    BONDS_FLAG: bool = True  # if wants to write parameters in the itp file
    ANGLES_FLAG: bool = True  # if wants to write parameters in the itp file
    DIHEDRALS_FLAG: bool = True  # if wants to write parameters in the itp file

    BONDS: dict[str, dict[str, typing.Any]]  # Name and length(nm) and strength
    BONDS = {
             'OMH_HO': {'r': 0.09590,
                        'k': 417560,
                        'funct': 1,
                        'source': 'doi.org/10.1021/la504178g'},

             'SU_OB': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts',
                       'note': '1 times of value of the source'},

             'SD_OD': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts',
                       'note': '1 times of value of the source'},

             'SD_OMH': {'r': 0.151,
                        'k': 417560,
                        'funct': 1,
                        'source': 'charmm from `top` scripts',
                        'note': '1 times of value of the source'},

             'SI_OB': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts',
                       'note': '1 times of value of the source'},

             'SI_OH': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts',
                       'note': '1 times of value of the source'},
            
             'SI_OR': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts',
                       'note': '1 times of value of the source'},

             'OH_HO': {'r': 0.09590,
                       'k': 417560,
                       'funct': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'OR_HR': {'r': 0.09590,
                       'k': 417560,
                       'funct': 1,
                       'source': 'doi.org/10.1021/la504178g'},
             
             'OR_HO': {'r': 0.09500,
                       'k': 0.0,
                       'funct': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'SD_OM': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts'},

             'SD_OR': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts'},

             'SI_OM': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts'},

             'SI_OR': {'r': 0.151,
                       'k': 417560,
                       'funct': 1,
                       'source': 'charmm from `top` scripts'},

             'CH_HC': {'r': 0.10900,
                       'k': 0.0,
                       'funct': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CN_HC': {'r': 0.10900,
                       'k': 0.0,
                       'funct': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'N_HN': {'r': 0.10100,
                      'k': 0,
                      'funct': 1,
                      'source': 'ODA from lammps'},
             

             'N_HN3': {'r': 0.10100,
                      'k': 0,
                      'funct': 1,
                      'source': 'ODA from lammps'},

             'SI_CH': {'r': 0.1800,
                       'k':0,
                       'funct': 1,
                       'source': 'see chatGpt file'},

             'CH_CH': {'r': 0.15290,
                       'k': 0.4,
                       'funct': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CH_CN': {'r': 0.15290,
                       'k': 0.4,
                       'funct': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CN_N': {'r': 0.14480,
                      'k': 0,
                      'funct': 1,
                      'source': 'ODA from lammps'},
            
            'CN_NH2': {'r': 0.14480,
                      'k': 0,
                      'funct': 1,
                      'source': 'ODA from lammps'},
            
            'CN_NH3': {'r': 0.14480,
                      'k': 0,
                      'funct': 1,
                      'source': 'ODA from lammps'},
            
            'NH2_HN2': {'r': 0.14480,
                      'k': 0,
                      'funct': 1,
                      'source': 'ODA from lammps'},

            'NH3_HN3': {'r': 0.14480,
                      'k': 0,
                      'funct': 1,
                      'source': 'ODA from lammps'}
             }
    ANGLES: dict[str, dict[str, typing.Any]]  # Angles names, rad, strength
    ANGLES = {
              'SI_OM_SD': {'theta': 141.00,
                           'cth': 450.50,
                           'funct': 1,
                           'source': 'charmm from `top` scripts'},

              'SI_OB_SI': {'theta': 144.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SD_OMH_HO': {'theta': 119.52,
                            'cth': 0.00,
                            'funct': 1,
                            'source': 'doi.org/10.1021/la504178g'},

              'SU_OB_SI': {'theta': 144.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_OH_HO': {'theta': 119.52,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SD_OM_SI': {'theta': 144.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_OB_SU': {'theta': 144.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_CH_CH': {'theta': 120.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see the chatGpt'},

              'SD_CH_CH': {'theta': 120.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see the chatGpt'},

              'CH_CN_N': {'theta': 109.47,
                          'cth': 0.00,
                          'funct': 5,
                          'source': 'doi.org/10.1021/la504178g'},

              'CH_CH_CN': {'theta': 112.70,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'ODA from lammps'},
                
              'CH_CH_CT': {'theta': 112.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'ODA from lammps'},
                
              'HN_N_HN': {'theta': 106.4,
                          'cth': 0.00,
                          'funct': 5,
                          'source': 'ODA from lammps'},

              'HC_CH_HC': {'theta': 107.80,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CN_HC': {'theta': 107.80,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},

              'SD_CH_HC': {'theta': 120.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see the chatGpt'},

              'SI_CH_HC': {'theta': 109.5,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see the chatGpt'},

              'HC_CN_N': {'theta': 109.5,
                          'cth': 0.00,
                          'funct': 5,
                          'source': 'ODA from lammps'},

              'OB_SI_CH': {'theta': 120.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see chatGpt'},

              'OM_SI_CH': {'theta': 120.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see chatGpt'},

              'OM_SD_CH': {'theta': 120.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see chatGpt'},

              'OMH_SD_CH': {'theta': 120.00,
                           'cth': 0.00,
                           'funct': 1,
                           'source': 'see chatGpt'},

              'CN_N_HN': {'theta': 109.47,
                          'cth': 0.00,
                          'funct': 5,
                          'source': 'ODA from lammps'},

              'CT_N_HN': {'theta': 109.47,
                          'cth': 0.00,
                          'funct': 5,
                          'source': 'ODA from lammps'},

              'CH_CH_HC': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},
                
              'CH_CT_HC': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CH_CH': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},
                
              'HC_CT_HC': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CT_N': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'}, 

              'CH_CN_HC': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CH_CN': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},
                
              'HC_CH_CT': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'},

              'CH_CT_N': {'theta': 110.70,
                           'cth': 0.00,
                           'funct': 5,
                           'source': 'doi.org/10.1021/la504178g'}

              }
    DIHEDRLAS: dict[str, typing.Any]  # Dihedrals names and parameters
    DIHEDRLAS = {
                 'SI_CH_CH_CN': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'SD_CH_CH_CN': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'CH_CH_CN_N': {'funct': 1,
                                'C0': 0.00000,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'C4': 0,
                                'C5': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'SI_CH_CH_HC': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'SD_CH_CH_HC': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'SD_CH_CH_CT': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'SI_CH_CH_CT': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},


                 'HC_CH_CN_N': {'funct': 1,
                                'C0': 0.00000,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'C4': 0,
                                'C5': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'OB_SI_CH_CH': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SD_CH_CH': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SI_CH_CH': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OMH_SD_CH_HC': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OMH_SD_CH_CH': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'CH_CN_N_HN': {'funct': 1,
                                'C0': 0.00000,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'C4': 0,
                                'C5': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'CH_CN_CT_N': {'funct': 1,
                                'C0': 0.00000,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'C4': 0,
                                'C5': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'CH_CN_CT_N': {'funct': 1,
                                'C0': 0.00000,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'C4': 0,
                                'C5': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'CH_CH_CT_N': {'funct': 1,
                                'C0': 0.00000,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'C4': 0,
                                'C5': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'CH_CH_CN_HC': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'CH_CH_CT_HC': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CH_CN': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CH_CT': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'CH_CT_N_HN': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CT_N': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CT_HC': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CT_N_HN': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SI_CH_HC': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OB_SI_CH_HC': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SD_CH_HC': {'funct': 1,
                                 'C0': 0.00000,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CN_N_HN': {'funct': 1,
                                 'C0': 0.00000,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'C4': 0,
                                'C5': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CH_HC': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CN_HC': {'funct': 9,
                                 'C0': 0.62760,
                                 'C1': 1.88280,
                                 'C2': 0,
                                 'C3': -2.51040,
                                 'C4': 0,
                                 'C5': 0,
                                 'source': 'doi.org/10.1021/la504178g'}
                 }
