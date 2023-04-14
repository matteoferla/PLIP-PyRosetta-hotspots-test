import pyrosetta
pr_scoring = pyrosetta.rosetta.core.scoring
from typing import List

class PlipspotsDockerExtras:
    def get_constraint_cost(self, holo: pyrosetta.Pose, weight=1) -> float:
        """
        How many kcal/mol does it cost the holoenyme if the constraints are set to ``weight`` (default = 1).
        """
        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()  # noqa
        without = scorefxn(holo)
        scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, weight)
        avec = scorefxn(holo)
        return avec - without

    def get_constraint_scores(self, holo: pyrosetta.Pose):
        """
        Each row of the output is the score of that constraint split into the individual constraints.
        """
        scores: List[List[float]] = []
        for con in holo.constraint_set().get_all_constraints():
            if not isinstance(con, pr_scoring.constraints.AmbiguousConstraint):
                continue
            scores.append([icon.score(holo) for icon in con.member_constraints()])
