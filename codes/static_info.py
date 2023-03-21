"""All the constant values use in the silanization"""

import typing


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
    Coverage: float = 3.0
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


class AtomGroup:
    """list of name of the atoms to work with"""
    # Si groups to find them in shell and add APTES to them
    SiGroup: list[str] = ['SD', 'SI']
    # Oxygen groups which are bonded to Si on the shell, SHOULD NOT replace!
    OMGroup: list[str] = ['OM', 'OB']
    # Oxygen groups bonded to Si on the shell to drop
    OxGroup: list[str] = ['OD', 'OH', 'OMH']
    # Hydrogen groups bonded to the Ox groups to drop
    HyGroup: list[str] = ['HO']


class DataFile:
    """Get data directory which are used in the script"""
    APTES: str = '/scratch/saeed/MyScripts/np_silica/data/aminopropyl.data'


class Hydration:
    """set all the info for water box
    Limitation for the box are added to the maximum radius of the NP"""
    TOLERANCE: float = 2.0
    # Box dimensions
    # x
    X_MIN: float = -20.0
    X_MAX: float = 20.0
    # y
    Y_MIN: float = -20.0
    Y_MAX: float = 20.0
    # z
    Z_MIN: float = -20.0
    Z_MAX: float = 20.0
    # Constants
    WATER_DENSITY = 0.9998395  # g/ml
    AVOGADRO: float = 6.0221408e+23  # 1/mol
    WATER_MOLAR_MASS: float = 18.01528  # g/mol
    MASSES: dict[str, float] = {'HW': 1.0080,
                                'OW': 15.9994,
                                'Cl': 35.453,
                                'Na': 22.989769}
    # PACKMOL files
    WATER_PDB: str = '/scratch/saeed/MyScripts/np_silica/data/water.pdb'
    NA_PDB: str = '/scratch/saeed/MyScripts/np_silica/data/Na.pdb'
    CL_PDB: str = '/scratch/saeed/MyScripts/np_silica/data/Cl.pdb'
    INP_FILE: str = 'water_box.inp'
    OUT_FILE: str = 'water_box.pdb'
    WS_INP: str = 'water_silica.inp'  # Input for final water & silanized file
    GRO_PDB: str = 'silica_water.pdb'  # File for GROMACS input
    # PACKMOL lib
    PACKMOL: str = '/home/saeed/Downloads/packmol/packmol'


class PdbMass:
    """information for the masses section in the output file of the
    silanization"""
    silica_residue: str = 'SIL'  # Name of the residue for the silica NP
    aptes_residue: str = 'APT'  # Name of the residue for the APTES
    HO: dict[str, typing.Any] = {'Atoms_names': 'HO',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'oplsaa_xx'
                                 }
    OB: dict[str, typing.Any] = {'Atoms_names': 'OB',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'oplsaa_xx'
                                 }
    OH: dict[str, typing.Any] = {'Atoms_names': 'OH',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'oplsaa_xx'
                                 }
    OM: dict[str, typing.Any] = {'Atoms_names': 'OM',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'oplsaa_xx'
                                 }
    OMH: dict[str, typing.Any] = {'Atoms_names': 'OMH',
                                  'Residue': silica_residue,
                                  'Element_symbol': 'O',
                                  'RECORD': 'ATOM',
                                  'ff_type': 'oplsaa_xx'
                                  }
    OD: dict[str, typing.Any] = {'Atoms_names': 'OD',
                                 'Residue': silica_residue,
                                 'Element_symbol': 'O',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'oplsaa_xx'
                                 }
    SI: dict[str, typing.Any] = {'Atoms_names': 'SI',
                                 'Residue':  silica_residue,
                                 'Element_symbol': 'SI',
                                 'RECORD': 'ATOM',
                                 'ff_type':  'oplsaa_xx'
                                 }
    SU: dict[str, typing.Any] = {'Atoms_names': 'SU',
                                 'Residue':  silica_residue,
                                 'Element_symbol': 'SI',
                                 'RECORD': 'ATOM',
                                 'ff_type':  'oplsaa_xx'
                                 }
    SD: dict[str, typing.Any] = {'Atoms_names': 'SD',
                                 'Residue':  silica_residue,
                                 'Element_symbol': 'SI',
                                 'RECORD': 'ATOM',
                                 'ff_type':  'oplsaa_xx'
                                 }
    N: dict[str, typing.Any] = {'Atoms_names': 'N',
                                'Residue': aptes_residue,
                                'Element_symbol': 'N',
                                'RECORD': 'ATOM',
                                'ff_type': 'opls_287'
                                }
    CH: dict[str, typing.Any] = {'Atoms_names': 'CH',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'C',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'opls_135'
                                 }
    HC: dict[str, typing.Any] = {'Atoms_names': 'HC',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'opls_140'
                                 }
    CN: dict[str, typing.Any] = {'Atoms_names': 'CN',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'C',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'opls_292'
                                 }
    HN: dict[str, typing.Any] = {'Atoms_names': 'HN',
                                 'Residue': aptes_residue,
                                 'Element_symbol': 'H',
                                 'RECORD': 'ATOM',
                                 'ff_type': 'opls_290'
                                 }
    ATOMS: dict[str, dict[str, typing.Any]] = {'HO': HO,
                                               'OB': OB,
                                               'OH': OH,
                                               'OM': OM,
                                               'OMH': OMH,
                                               'OD': OD,
                                               'SI': SI,
                                               'SU': SU,
                                               'SD': SD,
                                               'N': N,
                                               'CH': CH,
                                               'HC': HC,
                                               'CN': CN,
                                               'HN': HN,
                                               }


class BoAnDi:
    """interaction for bonds, angles, and dihedrals in silanized
    nanoparticles"""
    BONDS: dict[str, dict[str, typing.Any]]  # Name and length(nm) and strength
    BONDS = {
             'OMH_HO': {'length': 1,
                        'strength': 1,
                        'source': 'doi.org/10.1021/la504178g'},

             'SU_OB': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'SD_OD': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'SD_OMH': {'length': 1,
                        'strength': 1,
                        'source': 'doi.org/10.1021/la504178g'},

             'SI_OB': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'SI_OH': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'OH_HO': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'SD_OM': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'SI_OM': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CH_HC': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CN_HC': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'N_HN': {'length': 1,
                      'strength': 1,
                      'source': 'doi.org/10.1021/la504178g'},

             'SI_CH': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CH_CH': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CH_CN': {'length': 1,
                       'strength': 1,
                       'source': 'doi.org/10.1021/la504178g'},

             'CN_N': {'length': 1,
                      'strength': 1,
                      'source': 'doi.org/10.1021/la504178g'}

             }
    ANGLES: dict[str, dict[str, typing.Any]]  # Angles names, rad, strength
    ANGLES = {
              'SI_OM_SD': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_OB_SI': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SD_OMH_HO': {'angle': 1,
                            'strength': 1,
                            'source': 'doi.org/10.1021/la504178g'},

              'SU_OB_SI': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_OH_HO': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SD_OM_SI': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_OB_SU': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_CH_CH': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SD_CH_CH': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'CH_CN_N': {'angle': 1,
                          'strength': 1,
                          'source': 'doi.org/10.1021/la504178g'},

              'CH_CH_CN': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'HN_N_HN': {'angle': 1,
                          'strength': 1,
                          'source': 'doi.org/10.1021/la504178g'},

              'HC_CH_HC': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CN_HC': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SD_CH_HC': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'SI_CH_HC': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CN_N': {'angle': 1,
                          'strength': 1,
                          'source': 'doi.org/10.1021/la504178g'},

              'OB_SI_CH': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'OM_SI_CH': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'OM_SD_CH': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'CN_N_HN': {'angle': 1,
                          'strength': 1,
                          'source': 'doi.org/10.1021/la504178g'},

              'CH_CH_HC': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CH_CH': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'CH_CN_HC': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'},

              'HC_CH_CN': {'angle': 1,
                           'strength': 1,
                           'source': 'doi.org/10.1021/la504178g'}

              }
    DIHEDRLAS: dict[str, typing.Any]  # Dihedrals names and parameters
    DIHEDRLAS = {
                 'SI_CH_CH_CN': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'SD_CH_CH_CN': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'CH_CH_CN_N': {'C0': 0,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'SI_CH_CH_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'SD_CH_CH_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CN_N': {'C0': 0,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'OB_SI_CH_CH': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SD_CH_CH': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SI_CH_CH': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'CH_CN_N_HN': {'C0': 0,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'CH_CH_CN_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CH_CN': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SI_CH_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OB_SI_CH_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'OM_SD_CH_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CN_N_HN': {'C0': 0,
                                'C1': 0,
                                'C2': 0,
                                'C3': 0,
                                'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CH_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'},

                 'HC_CH_CN_HC': {'C0': 0,
                                 'C1': 0,
                                 'C2': 0,
                                 'C3': 0,
                                 'source': 'doi.org/10.1021/la504178g'}
                 }
