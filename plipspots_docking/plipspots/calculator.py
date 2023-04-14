from .serial import SerialPLIPper

import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools

from functools import singledispatchmethod
from typing import Tuple, Dict, List, Union
from collections import Counter, defaultdict
from plip.structure.preparation import PDBComplex, PLInteraction
from openbabel.pybel import Atom, Residue
from openbabel.pybel import ob
from fragmenstein.victor import MinimalPDBParser
import warnings
class PlipspotsCalculator(SerialPLIPper):
    """
    Get the data needed for pyrosetta (different class: ``PlipspotsDocker``). See ``summarize()``.
    Note the instance attribute ``interactions`` containing all the ``summarize`` outputs.
    ``summarize`` summarises one interaction, while the attribute has all of them.
    _i.e._ ``_ = pd.Series.apply(plipspotscalculator); plipspotscalculator.interactions``

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # note: self.interactions contains repeats!
        self.interactions: List[Dict[str, Union[str, int, float]]] = []

    def __call__(self, mol) -> List[Dict[str, Union[str, int, float]]]:
        if mol is None or not isinstance(mol, Chem.Mol) or mol.GetNumAtoms() == 0:
            return []
        holo: str = self.plonk(mol)
        interaction_set: PLInteraction = self.get_interaction_set(holo)
        interxns: List = interaction_set.all_hbonds_ldon + interaction_set.all_hbonds_pdon
        summaries: List[Dict[str, Union[str, int, float]]] = [plipper.summarize(interxn) for interxn in interxns]
        self.interactions.extend(summaries)
        return summaries

    def summarize(self, interaction) -> Dict[str, Union[str, int, float]]:
        """
        Summarize a hydrogen bond interaction for the creation of constraints.

        :param interaction: an iteraction, not a set
        :type interaction: isinstance does not work as module 'plip.structure.detection' has no attribute 'hbond'
                            smells like a named tuple
        :return: something like the code below

        .. code-block::

            {'ligand_donor': False
             'resn': 'ARG',
             'resi': 61,
             'chain': 'A',
             'distance_ad': 3.836168661568465,
             'distance_ca': 4.84658591175273,
             'angle_dha': 122.43975227115855,
             'angle_lcan': 82.24542212804627,
             'angle_lcac': 157.1401711721203,
             'donor_atomname': ' NH1',
             'donor_atomtype': 'Nox',
             'acceptor_atomname': 'N ',
             'acceptor_atomtype': 'Nar',
             'lig_coordinates': [54.58, 67.818, 43.279]}

        What to use?

        * coordinate constraint —annoying to extrapolate closer
        * distance and angle constraints (D-H-A) easier to idealise
        * may require angle constraint of protein neighbor - protein interactor - ligand interactor
        * for coarse grain, N - CA - ligand interactor & C - CA - ligand interactor for a pseudo-dihedral

        But actually the latter is easier to implement in general.
        However, what does the angle mean when I idealise the distance and D-H-A?
        """
        if interaction.__class__.__name__ != 'hbond':
            warnings.warn(f'Interaction {interaction} is not a hydrogen bond')
            return {}
        # ligand donor
        if interaction.d.residue.name in ('LIG', 'UNL'):
            ligand_donor = True
            prot_obresidue = interaction.a.residue.OBResidue
            lig_obatom = interaction.d.OBAtom
        # ligand acceptor
        else:
            ligand_donor = False
            prot_obresidue = interaction.d.residue.OBResidue
            lig_obatom = interaction.a.OBAtom
        ca: ob.OBAtom = self.get_atom_by_atomname(prot_obresidue, 'CA')
        n: ob.OBAtom = self.get_atom_by_atomname(prot_obresidue, 'N')
        c: ob.OBAtom = self.get_atom_by_atomname(prot_obresidue, 'C')
        distance_ca: float = ca.GetDistance(lig_obatom)
        # angle between CA and ligand
        # if the donor is N it would still make sense
        angle_lcan: float = n.GetAngle(ca, lig_obatom)
        angle_lcac: float = c.GetAngle(ca, lig_obatom)
        # all values:
        return dict(ligand_donor=ligand_donor,
                    resn=interaction.restype,
                    resi=interaction.resnr,
                    chain=interaction.reschain,
                    distance_ad=interaction.distance_ad,
                    distance_ca=distance_ca,
                    angle_dha=interaction.angle,
                    angle_lcan=angle_lcan,
                    angle_lcac=angle_lcac,
                    donor_atomname=self.get_atomname(interaction.d),
                    donor_atomtype=interaction.dtype,
                    acceptor_atomname=self.get_atomname(interaction.a),
                    acceptor_atomtype=interaction.atype,
                    lig_coordinates=[lig_obatom.GetX(), lig_obatom.GetY(), lig_obatom.GetZ()]
                    )
