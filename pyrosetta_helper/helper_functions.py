from pyrosetta import *
from Bio.PDB import PDBParser, PDBIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio import Seq
from Bio import pairwise2
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

def get_per_residue_energies_df(p,scorefxn):
    for_pandas = []
    scorefxn(p)
    score_type_dataframe = {}
    #pandasDataFrame()
    weights_map = scorefxn.weights()
    energy = p.energies()
    
    for s in range(1, int(rosetta.core.scoring.end_of_score_type_enumeration)+ 1):
        score_type = rosetta.core.scoring.ScoreType(s)
        score_weight = weights_map[score_type]
        if score_weight:
           score_type_dataframe[score_type] = score_weight
 
    for resi in range(1,p.total_residue()+1):
        raw_score_total = 0.0
        total_score = 0.0
        for score_term in score_type_dataframe:
            st = score_term
            w = score_type_dataframe[score_term]
            score = energy.residue_total_energies(resi)[score_term]
            wt_score = w * score
            for_pandas.append({
             'Pose':resi,
             'ScoreType':str(st).split('.')[-1],
             'Weight': w,
             'RawScore': score,
             'Score':wt_score})
            raw_score_total += score
            total_score += wt_score
        for_pandas.append({
             'Pose':resi,
             'ScoreType':'total',
             'Weight': 0,
             'RawScore': raw_score_total,
             'Score':total_score})

    df = pandas.DataFrame(for_pandas)[['Pose','ScoreType','Weight','RawScore','Score']]
    return df

def compress_file(file_):
    import gzip
    import shutil
    import os
    print("Compressing File {} to {}".format(file_, file_+'.gz'))
    with open(file_, 'rb') as f_in, gzip.open(file_+'.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(file_)



if __name__ == "__main__":
    init()
    p = pose_from_pdb('2ny7_clean.pdb')
    print(pose_structure_df(p))