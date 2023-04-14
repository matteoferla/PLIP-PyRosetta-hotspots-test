import functools
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

logger = ph.configure_logger()
logger.setLevel(30)
Params.log.setLevel(30)  # redundant
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                     ex1=None,
                                                                     ex2=None,
                                                                     # mute='all',
                                                                     ignore_unrecognized_res=True,
                                                                     load_PDB_components=True,
                                                                     ignore_waters=False)
                                 )


class PlipspotsDockerPrep:
    """
    This docks using data from ``PlipspotsCalculator.interactions``.
    Which is assumed deduplicated... (TODO implement)
    """

    def __init__(self,
                 pdb: str,
                 interactions: List[Dict[str, Union[str, int, float]]],
                 resn: str = 'LIG',
                 max_con_value: int = 5,
                 idealized: bool = False,
                 preminimize_cycles: int = 0, preminimize_con_weight: int = 5):
        self.pose = pyrosetta.Pose()
        if pdb.find('\n') != -1:  # pdb block
            prc.import_pose.pose_from_pdbstring(self.pose, pdb)
        else:  # pdb filename
            prc.import_pose.pose_from_file(self.pose, pdb)
        assert self.pose.total_residue() > 0, f'Your PDB has issues?'
        self.resn = resn
        self.interactions = interactions
        self.ligand_donor_interactions = [intxn for intxn in self.interactions if intxn['ligand_donor']]
        self.ligand_acceptor_interactions = [intxn for intxn in self.interactions if not intxn['ligand_donor']]
        self.idealized = idealized
        self.max_con_value = int(max_con_value)  # int is confusing
        if preminimize_cycles:
            self.minimize(preminimize_cycles, preminimize_con_weight)

    @functools.cached_property
    def centroid(self) -> prn.xyzVector_double_t:
        # this is an average... not a least squares... bad.
        # todo fix
        # https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.numeric.html
        xyz = prn.xyzVector_double_t(0., 0., 0.)
        for intxn in self.interactions:
            xyz += prn.xyzVector_double_t(*intxn['lig_coordinates'])
        l = len(self.interactions)
        return prn.xyzVector_double_t(xyz.x / l, xyz.y / l, xyz.z / l)

    def minimize(self, cycles: int = 3, con_weight: int = 5) -> float:
        """
        Minimize constraining to starting coordinates
        Note that it adds cons... and then removes them
        """
        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()  # noqa
        scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, con_weight)
        ori_option = prbo.get_boolean_option('relax:constrain_relax_to_start_coords')
        prbo.set_boolean_option('relax:constrain_relax_to_start_coords', True)
        prp.relax.FastRelax.register_options()
        relax = prp.relax.FastRelax(scorefxn, cycles)
        relax.constrain_relax_to_start_coords(True)
        pre: float = scorefxn(self.pose)
        relax.apply(self.pose)
        if not ori_option:
            self.pose.remove_constraints()
            prbo.set_boolean_option('relax:constrain_relax_to_start_coords', False)
            prp.relax.FastRelax.register_options()
        return scorefxn(self.pose) - pre

    def prepare(self, mol: Chem.Mol):
        if mol is None or not isinstance(mol, Chem.Mol) or mol.GetNumAtoms() == 0:
            raise ValueError
        params = Params.from_mol(mol, name='LIG')
        holo: pyrosetta.Pose = self.make_holo(params)
        cons = self.add_constraints(holo)  # returned variable is unused as add adds to pose in place
        return holo

    def make_holo(self, params):
        holo = self.pose.clone()
        # is this a reg ResidueTypeSet or a GlobalResidueTypeSet (its subclass), ehh who cares here
        rts: pr_chem.ResidueTypeSet = params.add_residuetype(holo, True)
        lig: pr_conf.Residue = pr_conf.ResidueFactory.create_residue(rts.name_map(self.resn))
        # fix positions to centroid
        xyz = prn.xyzVector_double_t(0., 0., 0.)
        for i in range(1, lig.natoms() + 1):
            xyz += lig.xyz(i)
        l = lig.natoms()
        offset = prn.xyzVector_double_t(xyz.x / l, xyz.y / l, xyz.z / l)
        for i in range(1, lig.natoms() + 1):
            lig.set_xyz(i, lig.xyz(i) + self.centroid - offset)
        # add mol
        holo.append_residue_by_jump(lig, 1)
        return holo

    def add_constraints(self, holo: pyrosetta.Pose):
        lig_resi = self.get_ligand_index(holo)
        ligand: pr_conf.Residue = self.get_ligand_residue(holo)
        # ## classify the ligand atoms
        ligand_donor_ids: List[pyrosetta.AtomID] = []
        ligand_acceptor_ids: List[pyrosetta.AtomID] = []
        for i in range(1, ligand.natoms() + 1):
            # atom_name: str = ligand.atom_name(i)
            atom_type: pr_chem.AtomType = ligand.atom_type(i)
            if atom_type.is_hydrogen() or atom_type.is_virtual():
                continue
            if atom_type.is_donor():
                ligand_donor_ids.append(pyrosetta.AtomID(atomno_in=i, rsd_in=lig_resi))
            if atom_type.is_acceptor():  # not elif -- e.g. hydroxyl
                ligand_acceptor_ids.append(pyrosetta.AtomID(atomno_in=i, rsd_in=lig_resi))
        # ## make an ambiguous constraint per protein atom
        for c, interaction in enumerate(self.interactions):
            # housekeeping
            if self.idealized:
                d = 3.
            else:
                d = interaction['distance_ad']
            other_resi: int = holo.pdb_info().pdb2pose(res=interaction['resi'], chain=interaction['chain'])
            assert other_resi, f'{interaction["resn"]}{interaction["resi"]} does not exist in the pose!?'
            other_res: pr_conf.Residue = holo.residue(other_resi)
            assert other_res.name3() == interaction[
                'resn'], f'Was expecting {interaction["resn"]} got {holo.residue(other_resi).name3()} for residue PDB:{interaction["resi"]}/pose:{other_resi}'
            if interaction['ligand_donor']:
                other_name = interaction['acceptor_atomname']
            else:
                other_name = interaction['donor_atomname']
            other_id = pyrosetta.AtomID(atomno_in=other_res.atom_index(other_name), rsd_in=other_resi)
            # ### make cons
            # classic madness:
            n_donors = len(ligand_donor_ids)
            cons = pru.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t(n_donors + 1)
            for c, lig_id in enumerate(ligand_donor_ids):
                cons[c + 1] = pr_scoring.constraints.AtomPairConstraint(lig_id, other_id,
                                                                        pr_scoring.func.HarmonicFunc(x0_in=d,
                                                                                                     sd_in=0.2))
            # maxed out case:
            cons[n_donors + 1] = pr_scoring.constraints.ConstantConstraint(
                pr_scoring.func.ConstantFunc(self.max_con_value),
                pr_scoring.ScoreType.atom_pair_constraint
                )
            holo.add_constraint(pr_scoring.constraints.AmbiguousConstraint(cons))

    def get_ligand_index(self, holo: pyrosetta.Pose) -> int:
        # bad... but hey. It is obs the last one added...
        return holo.total_residue()

    def get_ligand_residue(self, holo: pyrosetta.Pose) -> pr_conf.Residue:
        resi: int = self.get_ligand_index(holo)
        res: pr_conf.Residue = holo.residue(resi)
        assert res.name3() == self.resn, f'Impossible: {res.name3()} ≠ {self.resn}'
        return res

    def ligand2pdbblock(self, holo: pyrosetta.Pose) -> str:
        """
        This is the most elegant way to get a pdb block of the ligand I know.
        It feels very hacky, but it works.
        """
        res = self.get_ligand_residue(holo)
        pose = pyrosetta.Pose()
        rts = pose.conformation().modifiable_residue_type_set_for_conf(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)
        rts.add_base_residue_type(pr_chem.MutableResidueType(res.type()))
        pose.append_residue_by_jump(res.clone(), 0)
        return ph.get_pdbstr(pose)
