from pyrosetta import *
from Bio.PDB import PDBParser, PDBIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio import Seq
from Bio import pairwise2
import math
import Bio
import pandas
import array


def get_pose_from_pdb_with_chain(path, chain):
    p = PDBParser()
    struct = p.get_structure('TEST', path)
    c = struct[0][chain]
    io = PDBIO()
    io.set_structure(c)
    # Yuck - we have to save in PDB state
    io.save('/tmp/mypdb.pdb')
    pose = pose_from_pdb('/tmp/mypdb.pdb')
    os.remove('/tmp/mypdb.pdb')
    return pose


def pose_structure_df(pose, display_residues=[]):
    """
    Extracts and displays various structural properties of the input  <pose>
        and its  <display_residues>  including:
            -PDB numbering
            -chain identification
            -sequence

    """
    # store the pose's number of residues, example Python syntax
    nres = pose.total_residue()

    # 1. obtain the pose's sequence
    sequence = pose.sequence()

    # 2. obtain a list of PDB numbering and icode as a single string
    pdb_info = pose.pdb_info()
    PDB_nums = [(str(pdb_info.number(i)) + pdb_info.icode(i)).strip()
                for i in range(1, nres + 1)]
    # 3. obtains a list of the chains organized by residue
    chains = [pdb_info.chain(i) for i in range(1, nres + 1)]
    # 4. extracts a list of the unique chain IDs
    unique_chains = []
    for c in chains:
        if c not in unique_chains:
            unique_chains.append(c)

    phis = [pose.phi(i) for i in range(1, nres + 1)]
    psis = [pose.psi(i) for i in range(1, nres + 1)]
    omegas = [pose.omega(i) for i in range(1, nres + 1)]

    # Secondrry structure
    DSSP = rosetta.protocols.moves.DsspMover()
    DSSP.apply(pose)    # populates the pose's Pose.secstruct
    ss = pose.secstruct()

    return_list = []
    for i in range(1, nres + 1):
        return_list.append(
            {'PDB': PDB_nums[i-1],
             'Chain': chains[i-1],
             'Pose': i,
             'Residue': sequence[i-1],
             'Phi': phis[i-1],
             'Psi': psis[i-1],
             'Omega': omegas[i-1],
             'SS': ss[i-1],
             'Resi3': pose.residue(i).name()})

    df = pandas.DataFrame(return_list)
    return df.set_index('Pose')

def compress_file(file_):
    '''
        Zips a file in place
    '''
    import gzip
    import shutil
    import os
    print("Compressing File {} to {}".format(file_, file_+'.gz'))
    with open(file_, 'rb') as f_in, gzip.open(file_+'.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(file_)

def score_pose_to_df(input_pose, score_function='ref2015'):
    '''
        Add score information to pose dataframe
    '''
    ###Pose structure dataframe from 
    pose_df = pose_structure_df(input_pose)
    ref2015_sf = create_score_function(score_function)
    ref2015_sf(input_pose)
    energies = input_pose.energies()
    residue_energies = [energies.residue_total_energy(i) for i in range(1, input_pose.total_residue() + 1)]
    pose_df['residue_energy'] = residue_energies

    weights = [pyrosetta.rosetta.core.scoring.ScoreType(s)
        for s in range(1, int(
            pyrosetta.rosetta.core.scoring.end_of_score_type_enumeration) + 1)
        if ref2015_sf.weights()[pyrosetta.rosetta.core.scoring.ScoreType(s)]]

    per_residue_unwiehgts = energies.residue_total_energies
    per_residue_weighted = []
    for residue_index in range(1,input_pose.total_residue() + 1):
        sum_ = 0.0
        for weight in weights:
        #print(weight)
            entry = {'Pose': residue_index,
             'score_type': str(weight).split('.')[-1],
         'score': per_residue_unwiehgts(residue_index)[weight] * ref2015_sf.weights()[weight]}
            per_residue_weighted.append(
            entry)
    per_res_df = pandas.DataFrame(per_residue_weighted)
    scored_info_df = pose_df.join(per_res_df.pivot('Pose','score_type','score'))
    return scored_info_df

def get_sphere_sasa(input_pose):
    '''return a list of all residues per residue sasa'''
    neighbor_counts = []
    p = input_pose
    num_residues = p.total_residue()
    for res_target in range(1,num_residues+1):
        neighbor_count = 0.0
        CB_atom = "CB"
        ##If glycine, then use 1HA 
        if p.residue(res_target).type().name1() == 'G':
            CB_atom = "1HA"
        for res_neighbor in range(1, num_residues+1):
            ##Dont measure if res_neighbor and target is self
            if res_neighbor == res_target:
                continue
            else:
                neighbor_atom = "CB"
                if p.residue(res_neighbor).type().name1() == 'G':
                    neighbor_atom = "1HA"
            distance = p.residue(res_target).xyz(CB_atom).distance(
                p.residue(res_neighbor).xyz(neighbor_atom))
            neighbor_count += 1.0/(1.0 + math.exp(1.0*(distance-9.0)))
        neighbor_counts.append(neighbor_count)
    return neighbor_counts



if __name__ == "__main__":
    init()
    p = pose_from_pdb('2ny7_clean.pdb')
    print(pose_structure_df(p))
