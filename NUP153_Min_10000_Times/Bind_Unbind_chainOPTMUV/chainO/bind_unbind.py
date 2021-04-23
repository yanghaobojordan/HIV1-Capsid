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
    filename = sys.argv[1]
    pose=pose_from_pdb(filename)
    scorefxn = get_fa_scorefxn()
    dumpfile = 'Folding_Output_G1418Y_bind_unbind.pdb'
    newfile=open("Folding_Output_G1418Y_bind_unbind.txt", "w")
    newfile.write(str(scorefxn(pose)))
    newfile.write('\n')
    setup_foldtree(pose, 'ABCDEFPMTUV_O', Vector1([1]))
    trans_mover = RigidBodyTransMover(pose, 1)
    trans_mover.step_size(100)
    trans_mover.apply(pose)
    newfile.write(str(scorefxn(pose)))
    newfile.close()
    pose.dump_pdb(dumpfile)
main()
