import pytest

from Environmental_PAH_Mutagenicity.read_ccris_data import CIRconvert

def test__cir_convert():

    smiles = CIRconvert('108-88-3')

    assert smiles == 'Cc1ccccc1'
