"""add water box for the NP after silanization,
It gets radius of the silanized NP and prepare a box of water with a
free space around the origin of the box to put the NP there"""

import water_box as wbox


if __name__ == '__main__':
    moles = wbox.NumberMols(radius=50)
    in_file = wbox.InFile(radius=50,
                          num_mols=moles.number_mols,
                          edge=moles.edge_cube)
    water_box = wbox.RunPackMol()
