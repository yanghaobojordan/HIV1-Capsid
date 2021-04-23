from pyrosetta import *
from pyrosetta import PyMOLMover
from pyrosetta.toolbox import cleanATOM
from pyrosetta.toolbox import get_secstruct
from pyrosetta.teaching import *
from pyrosetta.toolbox import get_hbonds
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.rigid import *
from pyrosetta.rosetta.protocols.docking import *
import sys
init()

def main():
    translation=1.0
    rotation=1.0
    jobs=500
    job_output = 'Folding_WT_8_dock_output'
    filename=sys.argv[1]
    pose=pose_from_pdb(filename)
    dock_jump=1
    setup_foldtree(pose, 'ABCDEFOPTUV_M', Vector1([dock_jump]))

    test_pose = Pose()
    test_pose.assign(pose)

    scorefxn_high = create_score_function('docking')
    dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)

    movemap = MoveMap()
    movemap.set_jump(dock_jump, True)
    minmover = MinMover()
    minmover.movemap(movemap)
    minmover.score_function(scorefxn_high)

    perturb = SequenceMover()
    perturb.add_mover(dock_pert)
    perturb.add_mover(minmover)

    dock_prot = DockMCMProtocol()
    dock_prot.set_movable_jumps(Vector1([1]))
    dock_prot.set_scorefxn(scorefxn_high)
    dock_prot.set_partners('ABCDEFOPTUV_M')

    jd = PyJobDistributor(job_output, jobs, scorefxn_high)
    temp_pose = Pose()
    temp_pose.assign(pose)
    jd.native_pose = temp_pose

    while not jd.job_complete:
        try:
            test_pose.assign(pose)
            perturb.apply(test_pose)
            dock_prot.apply(test_pose)
            jd.output_decoy(test_pose)
        except RuntimeError:
            print 'ERROR:NAN occured in H-bonding calculations'
main()
