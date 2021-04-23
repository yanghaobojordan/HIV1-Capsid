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
    scorefxn=get_fa_scorefxn()
    
    MC(pose, scorefxn)

def MC(pose, scorefxn):
    test=Pose()
    test.assign(pose)
    dumpfile = 'Minimize_WT_chainV.pdb'
    txtfile = 'Minimize_WT_chainV.txt'
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

    smallmover=SmallMover(mm, kT, 1) #1 is the number of moves
    #smallmover.angle_max(7)
    shearmover=ShearMover(mm, kT, 1) #1 is the number of moves
    #shearmover.angle_max(7)
    
    task_pack = standard_packer_task(test)
    task_pack.restrict_to_repacking()
    task_pack.or_include_current(True) 
    pack_mover=PackRotamersMover(scorefxn, task_pack)

    combined_mover = SequenceMover()
    combined_mover.add_mover(smallmover)
    combined_mover.add_mover(shearmover)
    combined_mover.add_mover(min_mover)
    trial_mover = TrialMover(combined_mover, mc)

    for i in range (20):
        pack_mover.apply(test)
        mc.boltzmann(test)
        newfile.write(str(i))
        newfile.write(' ')
        newfile.write(str(scorefxn(test)))
        newfile.write(' ')
        newfile.write(str(CA_rmsd(pose, test)))
        newfile.write('\n')
    mc.recover_low(test)
    newfile.write('Repacking Complete')
    newfile.write('    ')
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')
        
    for i in range(5000):
        trial_mover.apply(test)
        #mc.boltzmann(test)
        #print scorefxn(test), i
        newfile.write(str(scorefxn(test)))
        newfile.write(' ')
        newfile.write(str(i))
        newfile.write(' ')
        newfile.write(str(CA_rmsd(pose, test)))
        newfile.write('\n')
    
    mc.recover_low(test)
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')
    newfile.write(str(CA_rmsd(pose, test)))
    newfile.close()
    test.dump_pdb(dumpfile)
    print('Lowest Score ', scorefxn(test))
    print("Number of Acceptances: ", trial_mover.num_accepts())
    print("Acceptance Rate: ", trial_mover.acceptance_rate())

main()
