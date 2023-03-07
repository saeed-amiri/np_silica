"""add water box for the NP after silanization,
It gets radius of the silanized NP and prepare a box of water with a
free space around the origin of the box to put the NP there.
It should be able to run separately and within the silanization script
If executed independently, read the LAMMPS data file, make the box,
put them together, and write both the LAMMPS data file and the PDB
files.
If executed within the silanization script, it needs the data before
writing into the output file."""

import sys
import my_tools
import water_box as wbox
import write_lmp as wrlmp
import combine_all as merge
import read_lmp_data as rdlmp
import read_water_box as rbox


if __name__ == '__main__':
    fname: str = sys.argv[1]
    nano_p = rdlmp.ReadData(fname)
    print(nano_p.Net_charge)
    radius: float = my_tools.get_radius(nano_p.Atoms_df)
    moles = wbox.NumberMols(radius=radius)
    in_file = wbox.InFile(radius=radius,
                          num_mols=moles.number_mols,
                          edge=moles.edge_cube,
                          net_charge=nano_p.Net_charge)
    water_box = wbox.RunPackMol()
    read_box = rbox.GetWaterDf(in_file.num_water)
    combined_box = merge.MergeAll(read_box, nano_p)
    fout: str  # Name of the output file
    fout = f'boxed_{fname}'
    write_lmp = wrlmp.WriteLmp(combined_box, output=fout)
    write_lmp.write_lmp()
