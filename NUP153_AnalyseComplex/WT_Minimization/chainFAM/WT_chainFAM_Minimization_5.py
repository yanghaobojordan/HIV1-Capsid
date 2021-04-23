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
    filename=sys.argv[1]
    pose=pose_from_pdb(filename)
    test=Pose()
    test.assign(pose)
    scorefxn=get_fa_scorefxn()

    dumpfile = 'WT_chainFAM_Minimization_5.pdb'
    txtfile = 'WT_chainFAM_Minimization_5.txt'
    newfile = open(txtfile, "w")
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')
    kT = 1
    mc = MonteCarlo(test, scorefxn, kT)
    min_mover = MinMover()
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
 
    min_mover.movemap(mm)
    min_mover.score_function(scorefxn)
    min_mover.min_type("dfpmin")
    min_mover.tolerance(0.001)
    task_pack=standard_packer_task(test)
    task_pack.restrict_to_repacking()
    task_pack.or_include_current(True)
    pack_mover=PackRotamersMover(scorefxn, task_pack)
    
    for i in range(20):
        pack_mover.apply(test)
        mc.boltzmann(test)
        newfile.write(str(i))
        newfile.write(' ')
        newfile.write(str(scorefxn(test)))
        newfile.write(' ')
        newfile.write(str(CA_rmsd(pose, test)))
        newfile.write('\n')
    mc.recover_low(test)
    print ('Repacking Complete')
    print ('Lowest Score ', scorefxn(test))
    print (mc.show_scores())
    print (mc.show_counters())
    print (mc.show_state())
    newfile.write('Repacking Complete')
    newfile.write('    ')
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')

    for i in range(10000):
        min_mover.apply(test)
        mc.boltzmann(test)
        newfile.write(str(i))
        newfile.write(' ')
        newfile.write(str(scorefxn(test)))
        newfile.write(' ')
        newfile.write(str(CA_rmsd(pose, test)))
        newfile.write('\n')
    mc.recover_low(test)
    print ('Minimization Complete')
    print ('Lowest Score ', scorefxn(test))
    print (mc.show_scores())
    print (mc.show_counters())
    print (mc.show_state())
    newfile.write('Minimization Complete')
    newfile.write('    ')
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')
    newfile.write('RMSD')
    newfile.write('    ')
    newfile.write(str(CA_rmsd(pose, test)))
    newfile.write('\n')
    newfile.close()
    test.dump_pdb(dumpfile)
     
main()
