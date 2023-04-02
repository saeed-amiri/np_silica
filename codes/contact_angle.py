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
    angle_rad: float  # Contact angle in radian
    angle_rad = np.radians(stinfo.Hydration.CONATCT_ANGLE)
    return radius * np.tan(angle_rad/2)


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
        self.cube_edge: float  # Edge of the cube that inscribed the NP
        self.get_numbers(radius, net_charge)

    def get_numbers(self,
                    radius: float,  # Radius of the silanized nanoparticle
                    net_charge: float  # Charge of the silanized NP with sign!
                    ) -> None:
        """clculate the numbers of each moles if asked"""
        box_volume: float  # Volume of the final system's box
        box_volume = self.__box_volumes(radius)

    def __box_volumes(self,
                      radius: float  # Radius of the silanized nanoparticle
                      ) -> float:
        sphere_volume: float  # Volume of the sphere (NP apprx. with sphere)
        box_volume: float  # Volume of the final system's box
        sphere_volume = self.__get_sphere_volume(radius)
        self.cube_edge, box_volume = self.__get_box_volume(sphere_volume)
        return box_volume

    def __get_box_volume(self,
                         sphere_volume: float  # Volume of the sphere
                         ) -> tuple[float, float]:
        """calculate the volume of the box including sphere's area
        For the largest possible sphere is inscribed in cube, the ratio
        of volumes is: V_sphere/V_cube = pi/6"""
        v_inscribed_box: float = 6*sphere_volume/np.pi
        cube_edge: float  # Edge of the cube that inscribed the sphere
        cube_edge = v_inscribed_box**(1/3)
        x_lim: float = (stinfo.Hydration.X_MAX -
                        stinfo.Hydration.X_MIN) + cube_edge
        y_lim: float = (stinfo.Hydration.Y_MAX -
                        stinfo.Hydration.Y_MIN) + cube_edge
        z_lim: float = (stinfo.Hydration.Z_MAX -
                        stinfo.Hydration.Z_MIN) + cube_edge
        box_volume: float = x_lim*y_lim*z_lim
        if box_volume <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                     f'\tZero volume, there in problem in setting box '
                     f'limitaion, box_volume is "{box_volume:.3f}"'
                     f'{bcolors.ENDC}')
        return cube_edge, box_volume

    def __get_sphere_volume(self,
                            radius: float  # Radius of the NP after silanizatio
                            ) -> float:
        """calculate the volume of the sphere for the NP"""
        sphere_volume: float = 4*np.pi*(radius**3)/3
        if sphere_volume <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                     f'\tZero volume, there in problem in setting sphere '
                     f'valume, it is "{sphere_volume:.3f}"'
                     f'{bcolors.ENDC}')
        return sphere_volume


if __name__ == '__main__':
    # sys_box = BoxEdges(radius=25)
    mole_nums = NumMols(radius=25, net_charge=0)
