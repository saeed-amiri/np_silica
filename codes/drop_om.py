import numpy as np
import pandas as pd
import update_coords as upcord
from colors_text import TextColor as bcolors


class DropOM:
    """drop OM in aminos which have less than three OM atoms
    The extra atoms will drop, and bonds, angles, dihedrals
    must be updated"""
    def __init__(self) -> None:
        pass
