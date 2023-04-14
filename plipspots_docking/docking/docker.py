from plipspots_docking.docking.prep import PlipspotsDockerPrep
from .extras import PlipspotsDockerExtras
# PlipspotsDockerExtras inherits from PlipspotsDockerPrep

import functools
import warnings
import pyrosetta
import pyrosetta_help as ph
from types import ModuleType
from typing import List, Dict, Union

prc: ModuleType = pyrosetta.rosetta.core
prn: ModuleType = pyrosetta.rosetta.numeric
pru: ModuleType = pyrosetta.rosetta.utility  # noqa
prbo: ModuleType = pyrosetta.rosetta.basic.options
pr_chem: ModuleType = pyrosetta.rosetta.core.chemical
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
prp: ModuleType = pyrosetta.rosetta.protocols
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
from rdkit import Chem
from rdkit_to_params import Params

class PlipspotsDocker(PlipspotsDockerPrep, PlipspotsDockerExtras):

    def dock(self):
        # for now
        warnings.warn('This will change in behaviour soon! via Relax for now', DeprecationWarning)
        return self.dock_via_relax()
    def dock_via_relax(self, holo: pyrosetta.Pose) -> float:

        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()  # noqa
        scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, 5)

        resn_sele = prc.select.residue_selector.ResidueNameSelector('LIG')
        neigh_sele = prc.select.residue_selector.NeighborhoodResidueSelector(resn_sele, distance=8,
                                                                             include_focus_in_subset=True)
        n = neigh_sele.apply(holo)
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(allow_chi=n)
        movemap.set_jump(True)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.apply(holo)
        return scorefxn(holo)