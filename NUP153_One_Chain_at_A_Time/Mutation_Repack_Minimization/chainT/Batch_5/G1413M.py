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
    
    mutate_residue(pose, pose.pdb_info().pdb2pose('T', 1413), "M")
    MC(pose, scorefxn, "M")

def MC(pose, scorefxn, mutant):
    test=Pose()
    test.assign(pose)
    dumpfile = 'Folding_Output_G1413'+str(mutant)+'.pdb'
    txtfile = 'Folding_Output_G1413'+str(mutant)+'.txt'
    moveList= 'Folding_Output_G1413'+str(mutant)+'_MoveList.txt'
    move_list_file=open(moveList, 'w')
    newfile = open(txtfile, "w")
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')
    kT = 1
    mc = MonteCarlo(test, scorefxn, kT)

    count=0
    move_list=[]

    residue=int(test.pdb_info().pdb2pose('T', 1413))
    residue=test.residue(residue).xyz("CA")
    for i in range(1, test.total_residue()+1):
        i_residue=test.residue(i).xyz("CA")
        if (residue-i_residue).norm()<10:
            move_list.append(i)
            count +=1

    move_list_file.write(str(count))
    move_list_file.write('\n')
    for i in move_list:
        move_list_file.write(str(pose.pdb_info().pose2pdb(i)))
        move_list_file.write('    ')
        move_list_file.write(pose.residue(i).name())
        move_list_file.write('\n')        
    move_list_file.close()

    min_mover = MinMover()
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    for i in move_list:
        mm.set_bb(i, True)
        mm.set_chi(i, True)
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
    task_pack.temporarily_fix_everything()
    for i in move_list:
        task_pack.temporarily_set_pack_residue(i,True)
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
