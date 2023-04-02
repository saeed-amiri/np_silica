"""This script reads the contact angle and calculates the part of the
nanoparticle in the oil phase. It will calculate the volume of the
water and oil parts and the number of oil and water molecules and
return them.
If the value for contact angle is a negative number, that means there
is no oil phase in the system.
"""

import sys
import numpy as np
import static_info as stinfo
from colors_text import TextColor as bcolors


def oil_depth(radius: float  # Radius of the nanoparticle
              ) -> float:
    """calculate and return the the depth of oil phase `h` in water"""
    return radius * np.tan(stinfo.Hydration.CONATCT_ANGLE/2)


class BoxEdges:
    """make the calculation for system box based on the contact angle"""
    def __init__(self,
                 radius: float  # Radius of NP after silanization
                 ) -> None:
        self.get_box_edges(radius)

    def get_box_edges(self,
                      radius: float  # Radius of NP after silanization
                      ) -> None:
        """calculate the edges of the water and oil phases"""


class NumMols:
    """calculate the number of molecules for each of the system:
    Numbers of water molecules (In data files named as: SOL)
    Number of decane molecules (In data files named as: D10)
    Numbers of ODAP molecules (In data files named as: ODA)
    Numbers od ODA molecules (In data files named as: ODN)
    Numbers of ION atoms (In data files named as: NA or CL)
    """
    def __init__(self,
                 radius: float,  # Radius of the silanized nanoparticle (NP)
                 net_charge: float  # Charge of the silanized NP with sign!
                 ) -> None:
        self.sol_num: int  # Number of water molecules
        self.d10_num: int  # Number of decane molecules
        self.oda_num: int  # Number of ODAp molecules
        self.odn_num: int  # Number of ODA molecules
        self.ion_num: int  # Number of ion atoms with sign!


if __name__ == '__main__':
    sys_box = BoxEdges(radius=25)
