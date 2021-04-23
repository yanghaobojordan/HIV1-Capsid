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
    job_output = 'Folding_Output_P1411Y_dock_output'
    filename=sys.argv[1]
    pose=pose_from_pdb(filename)
    dock_jump=1
    setup_foldtree(pose, 'ABCDEFPMTUV_O', Vector1([dock_jump]))
    to_centroid = SwitchResidueTypeSetMover('centroid')
    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    recover_sidechains = ReturnSidechainMover(pose)

    to_centroid.apply(pose)

    test_pose = Pose()
    test_pose.assign(pose)

    scorefxn_low = create_score_function('interchain_cen')
    scorefxn_high = create_score_function('docking')
    scorefxn_high_min = create_score_function('docking', 'docking_min')


    #use RigidBodyPerturbMover, can opt out RigidBodyTransMover, RigidBodySpinMover
    #use FaDockingSlideIntoContact to make sure two proteins are in contact
    #use MinMover
    #randomize_upstream = RigidBodyRandomizeMover(pose, dock_jump, partner_upstream)
    #randomize_downstream = RigidBodyRandomizeMover(pose, dock_jump, partner_downstream)
    dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
    #spin = RigidBodySpinMover(dock_jump)
    #trans_mover = RigidBodyTransMover(pose, dock_jump)
    #slide_into_contact = DockingSlideIntoContact(dock_jump)

    movemap = MoveMap()
    movemap.set_jump(dock_jump, True)
    minmover = MinMover()
    minmover.movemap(movemap)
    minmover.score_function(scorefxn_high_min)

    perturb = SequenceMover()
    #perturb.add_mover(randomize_upstream)
    #perturb.add_mover(randomize_downstream)
    perturb.add_mover(dock_pert)
    #perturb.add_mover(spin)
    #perturb.add_mover(trans_mover)
    #perturb.add_mover(slide_into_contact)
    perturb.add_mover(to_fullatom)
    perturb.add_mover(recover_sidechains)
    perturb.add_mover(minmover)

    dock_prot = DockingProtocol()
    dock_prot.set_movable_jumps(Vector1([1]))
    dock_prot.set_lowres_scorefxn(scorefxn_low)
    dock_prot.set_highres_scorefxn(scorefxn_high)
    dock_prot.set_partners('ABCDEFPTMUV_O')

    jd = PyJobDistributor(job_output, jobs, scorefxn_high)
    temp_pose = Pose()
    temp_pose.assign(pose)
    to_fullatom.apply(temp_pose)
    recover_sidechains.apply(temp_pose)
    jd.native_pose = temp_pose

    while not jd.job_complete:
        try:
            test_pose.assign(pose)
            perturb.apply(test_pose)
            dock_prot.apply(test_pose)
            to_fullatom.apply(test_pose)
            jd.output_decoy(test_pose)
        except RuntimeError:
            print 'ERROR:NAN occured in H-bonding calculations'
main()
