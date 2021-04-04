
from pyrosetta import *
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta import PyMOLMover
from pyrosetta.rosetta import *
from pyrosetta.toolbox.py_jobdistributor import output_scorefile
import numpy as np
import json
import time

import MutateInterface

def fp_dock_init(addl_flags=None):
    """
    Use to init parameters for the flexpep run. The ones below are currently working reasonably. 
    Best guide to the params that I have found is here:
    
    https://new.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock
    
    """
    opts = "-pep_refine -ex1 -ex2aro -use_input_sc -ignore_unrecognized_res -mute all"
    if addl_flags:
        opts += " "+addl_flags
    init(opts)

# copied from D090_Ala_scan.py in the example scripts for PyRosetta
def calc_binding_energy(pose, scorefxn, to_repack, cutoff = 8.0):
    scorefxn = create_score_function(scorefxn)
    
    # create a copy of the pose for manipulation
    test_pose = Pose()
    test_pose.assign(pose)

    # setup packer options
    # the sidechain conformations of residues "near the interface", defined as
    #    within  <cutoff>  Angstroms of an interface residue, may change and
    #    must be repacked, if all residues are repacked, aberrant sidechain
    #    conformations near the interface, but independent of complex
    #    interactions, will be repacked for the mutant and wild-type structures
    #    preventing them from adding noise to the score difference
    # this method of setting up a PackerTask is different from packer_task.py
    tf = standard_task_factory()    # create a TaskFactory
    tf.push_back(core.pack.task.operation.RestrictToRepacking())    # restrict it to repacking

    # this object contains repacking options, instead of turning the residues
    #    "On" or "Off" directly, this will create an object for these options
    #    and assign it to the TaskFactory
    prevent_repacking = core.pack.task.operation.PreventRepacking()
    
    for residue in to_repack:
        # the "center" (nbr_atom) of the mutant residue, for distance calculation
        center = test_pose.residue(residue).nbr_atom_xyz()
        for i in range(1, test_pose.total_residue() + 1):
            # the .distance_squared method is (a little) lighter than .norm
            # if the residue is further than <cutoff> Angstroms away, do not repack
            if center.distance_squared(
                    test_pose.residue(i).nbr_atom_xyz()) > cutoff**2:
                prevent_repacking.include_residue(i)

    # apply these settings to the TaskFactory
    tf.push_back(prevent_repacking)

    # setup a PackRotamersMover
    # changed by colin
    packer = protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)

    #### create a Mover for performing translation
    #### RigidBodyTransMover is SUPPOSED to translate docking partners of a
    ####    pose based on an axis and magnitude
    #### test it using the PyMOLMover, it does not perform a simple translation
    ####    I also observed a "Hbond Tripped" error when packing after applying
    ####    the Mover, it appears to store inf and NaN values into hbonds
    #transmover = RigidBodyTransMover()
    # calc_interaction_energy separates the chains by 500.0 Angstroms,
    #    so does this Mover
    # if using this Mover, the step_size MUST be a float
    # if this setting is left to default, it will move the proteins
    #    VERY far apart
    #transmover.step_size( 5.0 )

    # repack the test_pose
    packer.apply(test_pose)

    # score this structure
    before = scorefxn(test_pose)

    # separate the docking partners
    #### since RigidBodyTransMover DOES NOT WORK, it is not used
    #transmover.apply(test_pose)

    # here are two methods for applying a translation onto a pose structure
    # both require an xyzVector
    xyz = rosetta.numeric.xyzVector_double_t()    # a Vector for coordinates
    xyz.x = 500.0    # arbitrary separation magnitude, in the x direction
    xyz.y = 0.0    #...I didn't have this and it defaulted to 1e251...?
    xyz.z = 0.0    #...btw thats like 1e225 light years,
                   #    over 5e245 yrs at Warp Factor 9.999 (thanks M. Pacella)

    #### here is a hacky method for translating the downstream partner of a
    #    pose protein-protein complex (must by two-body!)
    chain2starts = len(pose.chain_sequence(1)) + 1
    for r in range(chain2starts, test_pose.total_residue() + 1):
        for a in range(1, test_pose.residue(r).natoms() + 1):
            test_pose.residue(r).set_xyz(a,
                test_pose.residue(r).xyz(a) + xyz)

    # here is an elegant way to do it, it assumes that jump number 1
    #    defines the docking partners "connectivity"
    # the pose.jump method returns a jump object CREATED from the pose jump
    #    data, the pose itself does not own a Jump object, thus you can use
    #    Jump methods, such as pose.jump(1).set_translation, however the object
    #    has not been properly constructed for manipulation, thus performing
    #    a change does not cause any problems, but is not permanently applied
    #translate = test_pose.jump( 1 )    # copy this information explicitly
    # adjust its translation via vector addition
    #translate.set_translation( translate.get_translation() + xyz )
    #test_pose.set_jump( 1 , translate )
    # as explained above, this call will NOT work
    #test_pose.jump(1).set_translation( test_pose.get_translation() + xyz )

    # repack the test_pose after separation
    packer.apply(test_pose)

    # return the change in score
    return before - scorefxn(test_pose)

# Also copied from D090_Ala_scan.py in the example scripts for PyRosetta
def ddG(wt_pose, mut_pose, to_repack, sf='docking', repack_distance=8, replicate_runs=20):
    
    a = []
    
    for _ in range(replicate_runs):
        wt_score = calc_binding_energy(wt_pose, sf,
            to_repack, repack_distance)
        mut_score = calc_binding_energy(mut_pose, sf,
            to_repack, repack_distance)
        #### the method calc_interaction_energy separates an input pose by
        ####    500 Angstroms along the jump defined in a Vector1 of jump numbers
        ####    for movable jumps, a ScoreFunction must also be provided
        #### if setup_foldtree has not been applied, calc_interaction_energy may be
        ####    wrong (since the jumps may be wrong)
        #wt_score = calc_interaction_energy(wt, scorefxn, movable_jumps)
        #mut_score = calc_interaction_energy(mutant, scorefxn, movable_jumps)
        ddg = mut_score - wt_score
        a.append(ddg)
    
    return np.mean(a)

def dG(pose, to_repack, sf='docking', repack_distance=8, replicate_runs=10):
    a = []
    
    for _ in range(replicate_runs):
        mut_score = calc_binding_energy(pose, sf,
            to_repack, repack_distance)
        a.append(mut_score)
    
    return np.mean(a)

def pep_run(decoy_name, n_decoys, pose, sf='docking', pymol_ip_addr=None):
    """
    Tested score function is 'docking' not sure if this is the best. Is there a way to run this in a 
    distributed fashion? Each decoy takes ~3 minutes on my machine.
    """
    if pymol_ip_addr:
        pmm = PyMOLMover(pymol_ip_addr, 65000) #enter the IP that pymol displays on startup
        pmm.keep_history(True)
    
    # Score function and starting PDB
    sf = create_score_function(sf) # no idea what this sf is...
    #pose = pose_from_pdb(pdb)

    # Creating FlexPepDock protocol using init options
    fpdock = FlexPepDockingProtocol()

    jd = PyJobDistributor(decoy_name, n_decoys, sf)
    while not jd.job_complete:
        pp = Pose()
        pp.assign(pose)
        fpdock.apply(pp)
        if pymol_ip_addr:
            pmm.apply(pp)
        jd.output_decoy(pp) # this will output a PDB file, which is not really necessary.

def pep_multicore_run(input_path, destination_path, filenames, nstructs):
    """
    Spreads a run onto as many cores as filenames. Good for looking at many peptides/mutants at once.
    This might need to get changed to input poses at some point?
    """
    import time
    from multiprocessing import Pool
    
    timestr = time.strftime("%Y%m%d-%H%M%S")

    with Pool() as p:
        work = [(destination_path+filename[:-4]+"-"+timestr+"_"+str(n), nstructs, None, None, pose_from_pdb(input_path+filename), 'score12') for n, filename in enumerate(filenames)]
        p.starmap(pep_run, work)
        
def fpdock(pose):
    fpdock = FlexPepDockingProtocol() # maybe some day I'll figure out how to add a scorefunction flag here??
    fpdock.apply(pose)
    
def add_score(pose, name, score):
    pose.scores[name] = score
    
def dump_scores(pose, path, i):
    name = pose.pdb_info().name().split('/')[2]
    fullpath = path+name+"_"+str(i)+".json"
    with open(fullpath, 'w') as json_file:
        json.dump(dict(pose.scores), json_file)
    
def dock_peptide(pdb_file, dump_path, scorefunction='docking', i='debug'):
    
    fp_dock_init("-score:weights "+scorefunction)
    
    pose = pose_from_pdb(pdb_file)
    
    fpdock(pose)
    scorefxn = create_score_function(scorefunction)
    s1 = scorefxn(pose)
    interface = MutateInterface.interface_res(pose, 0, False) # this will be different when mut has deletions
    d = dG(pose, interface, sf=scorefunction, replicate_runs=8)
    s2 = scorefxn(pose)
    add_score(pose, 'A_B binding energy', d)
    dump_scores(pose, dump_path, i)   
