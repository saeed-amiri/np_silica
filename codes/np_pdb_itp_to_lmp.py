import sys
import re
import itertools
import pandas as pd
import numpy as np

class Doc:
    """Reading the silicas' NP `PDB` and `ITP` files generated from
    `XYZ` files by the scripts code of Delle Piane.
    I am still determining the code; it works for the 1 to 9 nm but
    crashes for the 10 nm particles.
    The output of the scripts does not follow standard `PDB` formats;
    we used another one for this script.
    Input:
        Name of the file without an extension; if it had an extension,
        it would be removed; both the `PDB` and `ITP` file is needed.
    Outputs:
        A LAMMPS data file
        A LOG file that reports the numbers of each particle, the total
        charge of each type, and the total charge of the system.
    """