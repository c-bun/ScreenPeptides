# Dev path

## TODO

- would using "design" be better than point mutants?
    - convert interface list to resfile
    - fpdock (25 times?), then design, then fpdock (25 times?)?
    
- add checking new protein-peptide pair against potential orthogonal pair. (testing the negative pair)
    - would this require way too much CPU time?
    - what is the minimum number of decoys to get a good fpdock? Maybe ~25?
    
## DONE
    
- DONE - add checking ddG of binding with respect to WT configuration. Is this better than just looking at the total_energy in the FP dock result?
    - should probably look to see if this is correlated with kd first...
    
- DONE - add ability to run decoys in parallel. (Need new seed every time?)
    - take a look at the PyRosetta notebook for this.
    - is there a way to set up Dask so that the computer is recognised as a cluster?
    - how are seeds handled?
    
- DONE - add design of protein after peptide is docked (Packer task?)
    - is there a way to make sure it only mutates two positions? Find interface, choose two random positions, make resfile, pack with resfile. This should be doable looking at the AlaScan script and the PyRosetta notebooks.
    - design only protein side, but allow packing of peptide.