
from pyrosetta.rosetta import *
from pyrosetta import *
import ScreenPeptides
import MutateInterface

def dask_testmut(pdb_filepath):
    from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
    import ScreenPeptides
    import MutateInterface
    import random
    
    ScreenPeptides.fp_dock_init()

    pose = pose_from_pdb(pdb_filepath)
    position = random.choice(MutateInterface.interface_res(pose))
    aa = MutateInterface.choose_AA()
    toolbox.mutants.mutate_residue(pose, position, aa, pack_radius=8)
    
    fpdock = FlexPepDockingProtocol()
    fpdock.apply(pose)
    
    return pose.scores

def ddG_peptide(args):
    """
    Run FpDock and then calc ddG.
    
    First arg is path to pdb file. Second arg is path to native pdb file.
    """
    from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
    import ScreenPeptides
    
    pdb_filepath, pdb_native_filepath = args
    
    ScreenPeptides.fp_dock_init()

    pose = pose_from_pdb(pdb_filepath)
    native_pose = pose_from_pdb(pdb_native_filepath)
    
    fpdock = FlexPepDockingProtocol()
    fpdock.apply(pose)
    
    ddG = ScreenPeptides.ddG(native_pose, pose, 151)
    
    pose.scores['ddG'] = ddG
    pose.scores['fname'] = pose.pdb_info().name()
    
    return pose.scores
