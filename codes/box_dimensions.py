"""This script reads the contact angle and calculates the part of the
nanoparticle in the oil phase. It will calculate the volume of the
water and oil parts and the number of oil and water molecules and
return them.
If the value for contact angle is a negative number, that means there
is no oil phase in the system.
"""

import sys
import my_tools
import numpy as np
import static_info as stinfo
from colors_text import TextColor as bcolors


class NumMols:
    """calculate the number of molecules for each of the system:
    Numbers of water molecules (In data files named as: SOL)
    Number of decane molecules (In data files named as: oil)
    Numbers of ODAP molecules (In data files named as: ODAP)
    Numbers od ODA molecules (In data files named as: ODAN)
    Numbers of ION atoms (In data files named as: NA or CL)
    """
    def __init__(self,
                 radius: float,  # Radius of the silanized nanoparticle (NP)
                 net_charge: float  # Charge of the silanized NP with sign!
                 ) -> None:
        self.moles_nums: dict[str, int]  # All the needed moles and atoms
        self.moles_nums = {'sol': 0,  # Number of water molecues in the system
                           'oil': 0,  # Number of decane molecues in the system
                           'oda': 0,  # Number of ODAp molecues in the system
                           'odn': 0,  # Number of ODA molecues in the system
                           'ion': 0,  # Number of ION atoms in the system
                           'sal': 0  # Number of NaCl molecule in the system
                           }
        self.box_edges: dict[str, dict[str, float]]  # Edges of system's box
        self.box_edges = {'box': {'x_lim': 0.0,  # System's box
                                  'y_lim': 0.0,
                                  'z_lim': 0.0},
                          'sol': {'x_lim': 0.0,  # Water's section
                                  'y_lim': 0.0,
                                  'z_lim': 0.0},
                          'oil': {'x_lim': 0.0,  # Decane's section
                                  'y_lim': 0.0,
                                  'z_lim': 0.0}}
        self.inscibed_edge: float  # Edge of the cube that inscribed the sphere
        self.get_numbers(radius, net_charge)

    def get_numbers(self,
                    radius: float,  # Radius of the silanized nanoparticle
                    net_charge: float  # Charge of the silanized NP with sign!
                    ) -> None:
        """clculate the numbers of each moles if asked"""
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}:'
              f'({self.__module__}):'
              f'{bcolors.ENDC}')
        box_volume: float  # Volume of the final system's box
        box_volume = self.__box_volume(radius)
        # !!!num_oda must be set before num_ioins!!!
        self.moles_nums['oda'] = self.__get_odap_num()
        self.moles_nums['ion'] = self.__get_ion_num(net_charge)
        if stinfo.Hydration.CONATCT_ANGLE < 0:
            self.__pure_water_system(box_volume)
        else:
            print(f'{bcolors.OKCYAN}'
                  f'\tThe contact angle is: `{stinfo.Hydration.CONATCT_ANGLE}`'
                  ' [degree] which is: '
                  f'"{np.radians(stinfo.Hydration.CONATCT_ANGLE):.4f}" '
                  f'[rad] {bcolors.ENDC}')
            self.__oil_water_system(radius)
            self.moles_nums['odn'] = self.__get_odn_num()
        # Obviously the number of water should be defined before hand
        self.moles_nums['sal'] = int(self.__get_nacl_num())

    def __get_nacl_num(self) -> float:
        """return the number of NaCl molecules based on the concente-
        ration."""
        nacl = SaltSum(self.moles_nums['sol'])
        return nacl.n_nacl

    def __oil_water_system(self,
                           radius: float  # Radius of the silanized NP
                           ) -> None:
        """set the data for the system with oil and water"""
        oil_depth: float  # Depth of the oil phase on the NP
        oil_depth = my_tools.oil_depth(radius)
        self.__set_oil_water_edges(oil_depth)
        self.__get_oil_water_numbers(radius, oil_depth)

    def __get_oil_water_numbers(self,
                                radius: float,  # Radius of the silanized NP
                                oil_depth: float  # Depth of oil around NP
                                ) -> None:
        """set the numbers for the water and oil"""
        water_volume: float = self.__get_ow_volumes(
            self.box_edges['sol'].copy(), 'sol', radius, oil_depth)
        oil_volume: float = self.__get_ow_volumes(
            self.box_edges['oil'].copy(), 'oil', radius, oil_depth)
        self.moles_nums['sol'] = \
            self.__solution_mol_num(water_volume,
                                    stinfo.Hydration.WATER_DENSITY,
                                    stinfo.Hydration.WATER_MOLAR_MASS,
                                    'water')
        self.moles_nums['oil'] = \
            self.__solution_mol_num(oil_volume,
                                    stinfo.Hydration.OIL_DENSITY,
                                    stinfo.Hydration.OIL_MOLAR_MASS,
                                    'oil')

    def __get_ow_volumes(self,
                         edges: dict[str, float],  # Edges of the section
                         section: str,  # Name of section: water or oil
                         radius: float,  # Of silanized NP
                         oil_depth: float  # Depth of oil around NP
                         ) -> float:
        """return the volume of the oil and water phases"""
        box_volume: float  # Volume of the box with the area of NP in it
        box_volume = edges['x_lim']*edges['y_lim']*edges['z_lim']
        net_volume: float  # Net volume of the asked section
        sphere_volume: float = self.__get_sphere_volume(radius)
        sphere_in_water: float  # Volume in water section
        sphere_in_oil: float  # Volume in oil section
        sphere_in_oil = (1/3)*np.pi*(oil_depth**2)*(3*radius-oil_depth)
        sphere_in_water = sphere_volume - sphere_in_oil
        if section == 'sol':
            net_volume = box_volume - sphere_in_water
        elif section == 'oil':
            net_volume = box_volume - sphere_in_oil
        return net_volume

    def __set_oil_water_edges(self,
                              oil_depth: float  # Depth of oil around NP
                              ) -> None:
        """set the edges of the box for water and oil system"""
        self.box_edges['sol'] = self.box_edges['box'].copy()
        self.box_edges['oil'] = self.box_edges['box'].copy()
        self.box_edges['sol']['z_lim'] = \
            self.inscibed_edge - stinfo.Hydration.Z_MIN - oil_depth
        self.box_edges['oil']['z_lim'] = oil_depth+stinfo.Hydration.Z_MAX

    def __pure_water_system(self,
                            box_volume: float  # The volume of the system's box
                            ) -> None:
        """set the data for system with pure water"""
        self.moles_nums['oil'] = 0  # No oil in the system
        self.moles_nums['odn'] = 0  # ODA could be protonated or not
        self.moles_nums['sol'] = \
            self.__solution_mol_num(box_volume,
                                    stinfo.Hydration.WATER_DENSITY,
                                    stinfo.Hydration.WATER_MOLAR_MASS,
                                    'water')
        self.box_edges['sol'] = self.box_edges['box']

    def __get_ion_num(self,
                      net_charge: float,  # Net charge of the NP
                      ) -> int:
        """return the number of ions, with sign"""
        num_odap: int = self.moles_nums['oda']  # Number of ODAp
        charge_floor: float = np.floor(np.abs(net_charge))
        if stinfo.Hydration.ODAP_PROTONATION:
            num_ions: int = int(np.sign(net_charge)*charge_floor) + num_odap
        else:
            num_ions: int = int(np.sign(net_charge)*charge_floor)
        if charge_floor != np.abs(net_charge):
            print(f'{bcolors.CAUTION}{self.__class__.__name__}" '
                  f'({self.__module__}):\n'
                  f'\tNet charge is not a complete number! "{net_charge}"\n'
                  f'{bcolors.ENDC}')
        if num_ions > 0:
            print(f'{bcolors.CAUTION}'
                  f'\tTotal charge of the system is `{charge_floor}`\n'
                  f'\tThe number of ions is set to "{num_ions}"'
                  f'{bcolors.ENDC}')
        return num_ions

    def __solution_mol_num(self,
                           volume: float,  # Volume of the system
                           density: float,  # Density of the soultion
                           molar_mass: float,  # Self explanetory
                           sys_name: str  # Name of the section
                           ) -> int:
        """calculate the number of the water molecules int the volume"""
        lit_m3: float = 1e-24  # convert units
        m_water: float  # Mass of the water in the volume
        m_water = volume * density * lit_m3
        sol_moles: float
        sol_moles = int(m_water * stinfo.Hydration.AVOGADRO /
                        molar_mass) + 1
        print(f'{bcolors.OKCYAN}\tThe number {sys_name} molecules is '
              f'"{sol_moles}"{bcolors.ENDC}')
        return sol_moles

    def __get_odap_num(self) -> int:
        """calculate the number of oda based on the concentration"""
        oda_moles: int  # Number of oda molecules in the system
        # I will add the calculation later, but for now:
        state: str  # If ODA is protonated or not
        charge: int  # Charge of the ODA
        oda_moles = stinfo.Hydration.N_ODAP
        if stinfo.Hydration.ODAP_PROTONATION:
            state = 'protonated'
            charge = oda_moles
        else:
            state = 'unprotonated'
            charge = 0
        print(f'{bcolors.OKCYAN}'
              f'\tNumber of {state} ODA in water is set to '
              f'"{oda_moles}" with total charge of `{charge}`'
              f'{bcolors.ENDC}')
        return oda_moles

    def __get_odn_num(self) -> int:
        """calculate the number of odn based on the concentration in oil"""
        odn_moles: int  # Number of oda molecules in the system
        # I will add the calculation later, but for now:
        odn_moles = stinfo.Hydration.N_ODAN
        print(f'{bcolors.OKCYAN}\tNumber of ODAN (unprotonated ODA) is '
              f'set to "{odn_moles}" with total charge of `{odn_moles}`'
              f'{bcolors.ENDC}')
        return odn_moles

    def __box_volume(self,
                     radius: float  # Radius of the silanized nanoparticle
                     ) -> float:
        """get the box volume with subtracting the volume of NP"""
        sphere_volume: float  # Volume of the sphere (NP apprx. with sphere)
        box_volume: float  # Volume of the final system's box
        sphere_volume = self.__get_sphere_volume(radius)
        self.box_edges['box'], box_volume = \
            self.__get_box_volume(sphere_volume)
        return box_volume-sphere_volume

    def __get_box_volume(self,
                         sphere_volume: float  # Volume of the sphere
                         ) -> tuple[dict[str, float], float]:
        """calculate the volume of the box including sphere's area
        For the largest possible sphere is inscribed in cube, the ratio
        of volumes is: V_sphere/V_cube = pi/6"""
        v_inscribed_box: float = 6*sphere_volume/np.pi
        inscibed_edge: float  # Edge of the cube that inscribed the sphere
        inscibed_edge = v_inscribed_box**(1/3)
        box_edges: dict[str, float]  # Edges of the system box
        x_lim: float = (stinfo.Hydration.X_MAX -
                        stinfo.Hydration.X_MIN) + inscibed_edge
        y_lim: float = (stinfo.Hydration.Y_MAX -
                        stinfo.Hydration.Y_MIN) + inscibed_edge
        z_lim: float = (stinfo.Hydration.Z_MAX -
                        stinfo.Hydration.Z_MIN) + inscibed_edge
        box_volume: float = x_lim*y_lim*z_lim
        box_edges = {'x_lim': x_lim, 'y_lim': y_lim, 'z_lim': z_lim}
        if box_volume <= 0:
            sys.exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                     f'\tZero volume, there in problem in setting box '
                     f'limitaion, box_volume is "{box_volume:.3f}"'
                     f'{bcolors.ENDC}')
        self.inscibed_edge = inscibed_edge
        return box_edges, box_volume

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

    def print_info(self) -> None:
        """pylint"""


class SaltSum:
    """get the number of NaCl molecules based on the concentration
    that asked in the input."""
    def __init__(self,
                 n_water: int  # Number of the water molecules
                 ) -> None:
        self.n_nacl: int  # Number of the NaCl molecules
        self.n_nacl = self.get_salt(n_water)
        self.print_info()

    def get_salt(self,
                 n_water: int  # Number of the water molecules
                 ) -> int:
        """check the input data and select the module to call"""
        style: str  # The type of the concentration
        amount: float  # The amount which is asked for
        style = stinfo.Hydration.N_NACL['sty']
        amount = stinfo.Hydration.N_NACL['sum']
        n_nacl: int  # Number of the NaCl molecules
        if n_water == 0:
            n_nacl = int(amount)
        else:
            if amount == 0:
                n_nacl = int(amount)
            else:
                if style == 'mmolal':
                    n_nacl = self.get_molal(n_water, amount)
        return n_nacl

    def get_molal(self,
                  n_water: int,  # Number of the water molecules
                  amount: float  # Molality of the solution in millimoles/kg
                  ) -> int:
        """cfinding the number of NaCl molecules required to make a
        solution with a given molality in millimoles per kilogram of
        solvent (water),you did
        """
        n_nacl: float  # The claulated number of NaCl molecules
        n_nacl = n_water*stinfo.Hydration.WATER_MOLAR_MASS*amount*1e-6
        return np.ceil(n_nacl).astype(int)

    def print_info(self) -> None:
        """print info"""
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}: '
              f'({self.__module__})\n'
              '\tThe number of NaCl molecules for concentration of '
              f'`{stinfo.Hydration.N_NACL["sum"]} '
              f'{stinfo.Hydration.N_NACL["sty"]}` is set to '
              f'"{self.n_nacl}"\n{bcolors.ENDC}')


class BoxEdges:
    """make the calculation for system box based on the contact angle"""
    def __init__(self,
                 radius: float,  # Radius of NP after silanization
                 net_charge: float  # Total charge of the NP with sign!
                 ) -> None:
        self.radius: float = radius  # To have it as attribute for later use
        self.water_axis: dict[str, float] = {}  # Limitation of the water part
        self.oil_axis: dict[str, float] = {}  # Limitation of the oil part
        num_mols = NumMols(radius, net_charge)
        self.box_edges = num_mols.box_edges  # See NumMols class __init__
        self.num_mols = num_mols.moles_nums  # See NumMols class __init__
        self.get_sections_edge()
        self.print_info()

    def get_sections_edge(self) -> None:
        """set the limits of axis for each sections of the system box"""
        x_lo: float  # Low limit of the system box in x direction
        y_lo: float  # Low limit of the system box in y direction
        x_lo, y_lo = self.__get_xy_lims()
        self.__set_xy_lims(x_lo, y_lo)
        self.__set_z_lims()

    def __set_z_lims(self) -> None:
        """set limtations for the z axis in both water and oil"""
        z_w_lo: float  # Low limit water's section
        z_w_hi: float  # High limit water's section
        z_w_lo = -self.box_edges['box'].copy()['z_lim'] / 2
        z_w_hi = -z_w_lo - self.box_edges['oil'].copy()['z_lim']
        self.water_axis['z_lo'] = z_w_lo
        self.water_axis['z_hi'] = z_w_hi
        if stinfo.Hydration.CONATCT_ANGLE > 0:
            z_o_lo: float  # Low limit of the oil section
            z_o_hi: float  # High limit water's section
            z_o_lo = z_w_hi + stinfo.Hydration.TOLERANCE
            z_o_hi = -z_w_lo + stinfo.Hydration.TOLERANCE
            self.oil_axis['z_lo'] = z_o_lo
            self.oil_axis['z_hi'] = z_o_hi

    def __set_xy_lims(self,
                      x_lo: float,  # Low limit of system box in x direction
                      y_lo: float  # Low limit of system box in y direction
                      ) -> None:
        """set limitations in x and y"""
        self.water_axis['x_lo'] = x_lo
        self.water_axis['x_hi'] = -x_lo
        self.water_axis['y_lo'] = y_lo
        self.water_axis['y_hi'] = -y_lo
        if stinfo.Hydration.CONATCT_ANGLE > 0:
            self.oil_axis['x_lo'] = x_lo
            self.oil_axis['x_hi'] = -x_lo
            self.oil_axis['y_lo'] = y_lo
            self.oil_axis['y_hi'] = -y_lo

    def __get_xy_lims(self) -> tuple[float, float]:
        """find the x and y lowe limitations of the box"""
        x_lo: float  # min of the axis in the x_axis
        x_lo = -self.box_edges['box'].copy()['x_lim'] / 2
        y_lo: float  # min of the axis in the x_axis
        y_lo = -self.box_edges['box'].copy()['y_lim'] / 2
        return x_lo, y_lo

    def print_info(self) -> None:
        """pylit"""


if __name__ == '__main__':
    axis_limits = BoxEdges(radius=25, net_charge=10)
