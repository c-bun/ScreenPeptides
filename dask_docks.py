
from pyrosetta.rosetta import *
from pyrosetta import *

def dask_dock(**kwargs):
    from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
    
    import pyrosetta.distributed.io as io # Local import
    import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts # Local import
    
    pose = io.pose_from_file(kwargs["filepath"])
    #sf = create_score_function('docking') # no idea what this sf is...
    # Creating FlexPepDock protocol using init options
    fpdock = FlexPepDockingProtocol()
    fpdock.apply(pose)
    
    return pose
