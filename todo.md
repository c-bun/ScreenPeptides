# Dev path

- add design of protein after peptide is docked (Packer task?)
    - is there a way to make sure it only mutates two positions? Find interface, choose two random positions, make resfile, pack with resfile. This should be doable looking at the AlaScan script and the PyRosetta notebooks.
    - design only protein side, but allow packing of peptide.
    
- add ability to run decoys in parallel. (Need new seed every time?)
    - take a look at the PyRosetta notebook for this.
    - is there a way to set up Dask so that the computer is recognised as a cluster?
    - how are seeds handled?
    
- add checking new protein-peptide pair against potential orthogonal pair. (testing the negative pair)
    - would this require way too much CPU time?
    - what is the minimum number of decoys to get a good fpdock?