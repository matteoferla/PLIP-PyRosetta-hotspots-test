import unittest
from rdkit.Chem import PandasTools
from rdkit import Chem
from plipspots_docking.plipspots.calculator import PlipspotsCalculator


class PlipspotsTest(unittest.TestCase):

    def setUp(self):
        self.pdb_filename: str = 'NUDT7_reference.relax.pdb'
        self.sdf_filename: str = 'NUDT7_hits.sdf'
        self.hits = PandasTools.LoadSDF(self.sdf_filename)

    def test_something(self):
        raise WIPError
        mol: Chem.Mol = self.hits.ROMol[14]
        plipper = PlipspotsCalculator.from_filename(self.pdb_filename)
        self.assertEqual(len(plipper(mol)), 3)

if __name__ == '__main__':
    unittest.main()
