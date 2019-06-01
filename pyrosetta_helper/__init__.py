from pyrosetta import *
from Bio import Seq, pairwise2
from Bio.PDB import PDBParser, PDBIO
from Bio.SubsMat import MatrixInfo as matlist
from itertools import combinations
import math
import pandas
import array
import gzip
import shutil
import os


#INIT Package, so need this
from .GenKic import GenKic
from .BluePrintManager import BluePrintManager
from .BluePrintManager import BluePrintEntity

__all__ = ["GenKic","get_one_to_three","get_pose_from_pdb_with_chain",
"get_pose","pose_structure_df","compress_file","score_pose_to_df","score_interface_to_df",
"get_sphere_sasa","reverse_pose","slice_pose","join_poses","find_current_disulfides"]
#from . import DeNovoInterface
one_to_three = {
    'A': 'ALA',
    'C': 'CYS',
    'D': 'ASP',
    'E': 'GLU',
    'F': 'PHE',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'K': 'LYS',
    'L': 'LEU',
    'M': 'MET',
    'N': 'ASN',
    'P': 'PRO',
    'Q': 'GLN',
    'R': 'ARG',
    'S': 'SER',
    'T': 'THR',
    'V': 'VAL',
    'Y': 'TYR',
    'W': 'TRP'
}

def get_one_to_three(res):
    """Return the 3 letter code of a 1 letter amino acid
    
    Arguments:
        res {str} -- AA
    
    Returns:
        str -- three letter
    """
    return one_to_three[res]

def get_pose(pdb):
    """Get Pose Object from PDB Path
    
    Arguments:
        pdb {path} -- path of PDB file
    """
    return pose_from_pdb(pdb)

def get_pose_from_pdb_with_chain(path, chain):
    """Given a PDB file path, return a Pose Object of just the Chain
    
    Arguments:
        path {filepath} -- PDB File Path
        chain {str} -- Chain you want
        
    Returns:
        pyrosetta.Pose -- Pose Object
    """
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
    """Given a Pose Object, return a DataFrame of all info
    
    Arguments:
        pose {Pose Object}
    
    Keyword Arguments:
        display_residues {only list these residues given residue numbers} -- [residue numbers] (default: {[]})
    
    Returns:
        pandas.DataFrame
    """
    # store the pose's number of residues, example Python syntax
    nres = pose.total_residue()

    # 1. obtain the pose's sequence
    sequence = pose.sequence()

    # 2. obtain a list of PDB numbering and icode as a single string
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        print('setting PDB Info')
        #If its none for some reason, lets set it manually
        pdb_info = rosetta.core.pose.PDBInfo(pose)
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
    am = rosetta.core.sequence.ABEGOManager()
    abego_string = am.get_abego_string(am.get_symbols(pose))
    
    return_list = []
    for i in range(1, nres + 1):
        return_list.append(
            {'Pose Name': pdb_info.name(),
             'PDB': PDB_nums[i-1],
             'Chain': chains[i-1],
             'Pose': i,
             'Residue': sequence[i-1],
             'Phi': phis[i-1],
             'Psi': psis[i-1],
             'Omega': omegas[i-1],
             'SS': ss[i-1],
             'ABEGO':abego_string[i-1],
             'Resi3': pose.residue(i).name()})

    df = pandas.DataFrame(return_list)
    return df.set_index('Pose')

def compress_file(file_):
    """Compress File
    
    Arguments:
        file_ {path} -- file path
    """
    print("Compressing File {} to {}".format(file_, file_+'.gz'))
    with open(file_, 'rb') as f_in, gzip.open(file_+'.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(file_)

def score_pose_to_df(input_pose, score_function=''):
    """[summary]
    
    Arguments:
        input_pose {[type]} -- [description]
    
    Keyword Arguments:
        score_function {str} -- [description] (default: {''})
    
    Returns:
        [type] -- [description]
    """
    if isinstance(input_pose, str):
        input_pose = pose_from_file(input_pose)
    pose_df = pose_structure_df(input_pose)
    if not score_function:
        ref2015_sf = create_score_function('ref2015')
        score_function = ref2015_sf
    
    score_function(input_pose)
    energies = input_pose.energies()
    residue_energies = [energies.residue_total_energy(
        i) for i in range(1, input_pose.total_residue() + 1)]
    pose_df['residue_energy'] = residue_energies

    weights = [pyrosetta.rosetta.core.scoring.ScoreType(s)
               for s in range(1, int(
                   pyrosetta.rosetta.core.scoring.end_of_score_type_enumeration) + 1)
               if score_function.weights()[pyrosetta.rosetta.core.scoring.ScoreType(s)]]

    per_residue_unwiehgts = energies.residue_total_energies
    per_residue_weighted = []
    for residue_index in range(1, input_pose.total_residue() + 1):
        for weight in weights:
            entry = {'Pose': residue_index,
                     'score_type': str(weight).split('.')[-1],
                     'score': per_residue_unwiehgts(residue_index)[weight] * score_function.weights()[weight]}
            per_residue_weighted.append(
                entry)
    per_res_df = pandas.DataFrame(per_residue_weighted)
    scored_info_df = pose_df.join(
        per_res_df.pivot('Pose', 'score_type', 'score'))
    return scored_info_df

def score_interface_to_df(input_pose, interface_string, score_function='ref2015'):
    """[summary]
    
    Arguments:
        input_pose {[type]} -- [description]
        interface_string {[type]} -- [description]
    
    Keyword Arguments:
        score_function {str} -- [description] (default: {'ref2015'})
    
    Returns:
        [type] -- [description]
    """
    # Pose structure dataframe from
    ref2015_sf = create_score_function(score_function)
    if isinstance(input_pose, str):
        input_pose = pyrosetta.pose_from_pdb(input_pose)
    ref2015_sf(input_pose)
    df = score_pose_to_df(input_pose)

    # Lets use the IAMover
    iam = rosetta.protocols.analysis.InterfaceAnalyzerMover
    ia = iam(interface_string, False, ref2015_sf, True, True, True, False)
    ia.apply(input_pose)

    ri = ia.get_all_per_residue_data()
    df['InterFaceResidue'] = [i for i in ri.interface_residues]
    df['SeparatedSASA'] = [i for i in ri.separated_sasa]
    df['ComplexedSASA'] = [i for i in ri.complexed_sasa]
    df['dSASA'] = [i for i in ri.dSASA]
    df['dSASA_sc'] = [i for i in ri.dSASA_sc]
    df['dhSASA'] = [i for i in ri.dhSASA]
    df['dhSASA_sc'] = [i for i in ri.dhSASA_sc]
    df['dhSASA_rel_by_charge'] = [i for i in ri.dhSASA_rel_by_charge]
    df['SASA'] = [i for i in ri.SASA]
    df['dSASA_fraction'] = [i for i in ri.dSASA_fraction]
    df['separated_energy'] = [i for i in ri.separated_energy]
    df['complexed_energy'] = [i for i in ri.complexed_energy]
    df['ddG'] = [i for i in ri.dG]

    return df

def get_sphere_sasa(input_pose):
    """[summary]
    
    Arguments:
        input_pose {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    neighbor_counts = []
    p = input_pose
    num_residues = p.total_residue()
    for res_target in range(1, num_residues+1):
        neighbor_count = 0.0
        CB_atom = "CB"
        # If glycine, then use 1HA
        if p.residue(res_target).has('CB'):
            CB_atom = 'CB'
        elif p.residue(res_target).has('1HA'):
            CB_atom = '1HA'
        else:
            print("Target residue {} does not have CB or 1HA".format(
                p.residue(res_target)))
            continue
        # if p.residue(res_target).type().name1() == 'G':
        #    CB_atom = "1HA"
        for res_neighbor in range(1, num_residues+1):
            # Dont measure if res_neighbor and target is self
            if res_neighbor == res_target:
                continue
            else:
                if p.residue(res_neighbor).has('CB'):
                    neighbor_atom = "CB"
                elif p.residue(res_neighbor).has('1HA'):
                    neighbor_atom = "1HA"
                else:
                    print("Neighbor residue {} does not have CB or 1HA".format(
                        p.residue(res_neighbor)))
                    continue
                # p.residue(res_neighbor).type().name1() == 'G':
            distance = p.residue(res_target).xyz(CB_atom).distance(
                p.residue(res_neighbor).xyz(neighbor_atom))
            neighbor_count += 1.0/(1.0 + math.exp(1.0*(distance-9.0)))
        neighbor_counts.append(neighbor_count)
    return neighbor_counts

def reverse_pose(p):
    """[summary]
    
    Arguments:
        p {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    reversed_pose = Pose()
    reversed_pose.pdb_info(rosetta.core.pose.PDBInfo(reversed_pose))
    reversed_pose.pdb_info().name('Reverse'+p.pdb_info().name())

    for i in range(p.total_residue(),0,-1):
        reversed_pose.append_residue_by_bond(p.residue(i))
    return reversed_pose

def slice_pose(p,start,end,retain_pdb=True):

    sliced = Pose()
    disulfides = []
    if end > p.size() or start > p.size():
        return "end/start slice is longer than total lenght of pose {} {}".format(start,end) 
    for dis_index, i in enumerate(range(start,end+1),start=1):
        #print(p.residue(i).name())
        if p.residue(i).name() == 'CYS:disulfide':
            disulfides.append(dis_index)
        sliced.append_residue_by_bond(p.residue(i))


    pdb_info = rosetta.core.pose.PDBInfo(sliced)
    if retain_pdb:
        pdb_info.copy(p.pdb_info(),start,end,1)
    sliced.pdb_info(pdb_info)
    ##We could check if the disulfide mate is in here, but I think it's a little to complicated
    if disulfides:
        for disulfide in disulfides:
            print('unsetting {}'.format(disulfide))
            mr = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
            mr.set_target(disulfide)
            mr.set_res_name('CYS')
            mr.apply(sliced)
    return sliced

def join_poses(pose_A, pose_B):
    
    '''Take two poses and join them'''
    p = Pose()
    p.assign(pose_A)
    for index,residue in enumerate(range(1,pose_B.total_residue()+1),start=p.size()):
        p.append_residue_by_bond(pose_B.residue(residue))
    return p
     

def find_current_disulfides(p):
    """[summary]
    
    Arguments:
        p {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    disulfds = []   
    pairs = []
    for r in range(1,p.size()+1):
        if p.residue(r).type().is_disulfide_bonded():
            disulfds.append(r)
    for tupe in combinations(disulfds,2):
        if p.residue(tupe[0]).is_bonded(p.residue(tupe[1])):
            pairs.append(tupe)
    return pairs

def fill_pose_with_pdb(p,chain='A'):
    """Fill pose with PDB info
    
    Arguments:
        p {Pose} -- PyRosetta Pose Options 
    
    Keyword Arguments:
        chain {str} -- chain character (default: {'A'})
    """
    pdb_info = p.pdb_info()
    for pdb in range(1,p.total_residue()+1):
        pdb_info.set_resinfo(pdb,'A',pdb) 
