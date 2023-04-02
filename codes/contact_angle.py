"""This script reads the contact angle and calculates the part of the
nanoparticle in the oil phase. It will calculate the volume of the
water and oil parts and the number of oil and water molecules and
return them.
If the value for contact angle is a negative number, that means there
is no oil phase in the system.
"""


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


if __name__ == '__main__':
    sys_box = BoxEdges(radius=25)
