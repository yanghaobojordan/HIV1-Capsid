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
    pose=pose_from_pdb("Folding_WT_8.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_WT_8_bind_unbind.pdb'
    newfile=open("dg.txt", "w")
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_P1411Y.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_P1411Y_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_S1412P.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_S1412P_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_G1413W.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_G1413W_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_V1414W.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_V1414W_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_F1415G.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_F1415G_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_T1416R.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_T1416R_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_F1417G.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_F1417G_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_G1418Y.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_G1418Y_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_P1411M.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_P1411M_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_S1412M.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_S1412M_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_G1413M.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_G1413M_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_V1414I.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_V1414I_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_F1415M.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_F1415M_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_T1416M.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_T1416M_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_F1417Y.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_F1417Y_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    #####################################################
    pose=pose_from_pdb("Folding_Output_G1418A.pdb")
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_G1418A_bind_unbind.pdb'
    bind_score=float(scorefxn(pose))
    setup_foldtree(pose, 'ABCDEF_OPMTUV', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    unbind_score=float(scorefxn(pose))
    ddg_bind = bind_score-unbind_score
    newfile.write(str(ddg_bind))
    newfile.write('\n')
    pose.dump_pdb(dumpfile)
    newfile.close()
main()
