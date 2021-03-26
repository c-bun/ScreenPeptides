
# put it all together!!

import pyrosetta.rosetta.protocols as protocols
import random
    
def interface_res(pose, interface_distance=8, exclude_chain=True):
    
    interface = protocols.interface.select_interface_residues(pose, 'A_B', interface_distance)
    
    if exclude_chain_B:
        B_chain_start = nbit.chain_end(1)
        for i, res in enumerate(interface):
            if i > B_chain_start:
                interface[i] = 0
                
    # convert to a list of residue numbers
    interface_res = []
    for i, res in enumerate(interface):
        if res == 1:
            interface_res.append(i+1)
    
    return interface_res
#    return random.choice(interface_res)

def choose_AA(aas=["A", "V", "I", "L", "L", "F", "Y", "W", "S", "T", "N", "Q", "G", "R", "H", "K", "D", "E"]):
    return random.choice(aas)

def point_mutation_and_dock(pose, ndecoys):
    position = random.choice(interface_res(pose))
    aa = choose_AA()
    toolbox.mutants.mutate_residue(pose, position, aa, pack_radius=8)
    
    pep_run("./decoys/"+str(position)+aa, ndecoys, pose)
