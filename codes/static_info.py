class Doc:
    """All the constant values use in the silanization"""


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
    OM: float = -0.720
    SI: float = 1.400


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



class AtomGroup:
    """list of name of the atoms to work with"""
    # Si groups to find them in shell and add APTES to them
    SiGroup: list[str] = ['SD', 'SI']
    # Oxygen groups which are bonded to Si on the shell, SHOULD NOT replace!
    OMGroup: list[str] = ['OM', 'OB']
    # Oxygen groups bonded to Si on the shell to drop 
    OxGroup: list[int] = ['OD', 'OH', 'OMH']
    # Hydrogen groups bonded to the Ox groups to drop
    HyGroup: list[str] = ['HO']


class DataFile:
    """Get data directory which are used in the script"""
    APTES: str = '/scratch/saeed/MyScripts/np_silica/data/aminopropyl.data'
