# Dev path

- DONE - add design of protein after peptide is docked (Packer task?)
    - is there a way to make sure it only mutates two positions? Find interface, choose two random positions, make resfile, pack with resfile. This should be doable looking at the AlaScan script and the PyRosetta notebooks.
    - design only protein side, but allow packing of peptide.
    
### Solution:

- ```pyrosetta.toolbox.mutants.mutate_residue(pack_or_pose, mutant_position, mutant_aa, pack_radius=0.0, pack_scorefxn=None)```
- Will mutate selected residue and repack. Only works for one res at a time, but could just stack them?
- ```pyrosetta.rosetta.protocols.interface.select_interface_residues(pose: pyrosetta.rosetta.core.pose.Pose, interface: str, interface_distance: int) → pyrosetta.rosetta.utility.vector1_bool```
- Will return the interface residues for two interacting chains (A_B). Not sure what form the data comes back in, but could use to determine the residues on the protein side to mutate? Might be good to look at second sphere too though?
    
- add ability to run decoys in parallel. (Need new seed every time?)
    - take a look at the PyRosetta notebook for this.
    - is there a way to set up Dask so that the computer is recognised as a cluster?
    - how are seeds handled?
    
- add checking new protein-peptide pair against potential orthogonal pair. (testing the negative pair)
    - would this require way too much CPU time?
    - what is the minimum number of decoys to get a good fpdock?
    
- add checking ddG of binding with respect to WT configuration. Is this better than just looking at the total_energy in the FP dock result?
    - should probably look to see if this is correlated with kd first...