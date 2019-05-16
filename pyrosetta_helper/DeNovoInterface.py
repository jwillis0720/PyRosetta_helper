from itertools import combinations,combinations_with_replacement
import pyrosetta
import glob
import pandas

class DistanceSelector():
    
    def __init__(self,chain, residues=[]):
        
        self.chain = chain
        
        #Convert to list incase its a generator
        self.residues = [int(i) for i in residues]
        
    def get_chain(self):
        
        return self.chain
        
    def get_residues(self):
        
        return self.residues
    
    def __str__(self):
        return "Chain {} wit Residues {}".format(self.chain, self.residues)

def permutate_distance_constraints(distance_df, ds, dist='dist'):
    sub_df = []
    for combo in combinations(ds.get_residues(),2):
        chain_1 = ds.get_chain()
        chain_2 = ds.get_chain()
        pdb_1 = combo[0]
        pdb_2 = combo[1]
        row = distance_df[(distance_df['chain1'] == chain_1) &\
                               (distance_df['chain2'] == chain_2) &\
                               (distance_df['pdb1'] == str(pdb_1)) &\
                               (distance_df['pdb2'] == str(pdb_2))]
        sub_df.append(row)
    return pandas.concat(sub_df)

def get_all_distances(p):
    sequence = p.sequence()
    num_residues = p.total_residue()
    all_distances = []
    pdb_info = p.pdb_info()
    PDB_nums = [(str(pdb_info.number(i)) + pdb_info.icode(i)).strip()
                for i in range(1, num_residues + 1)]

    three_letter_seqs = [str(p.residue(i).name3()) for i in range(1,num_residues + 1)]
    chains = [pdb_info.chain(i) for i in range(1, num_residues + 1)]
    for res_1 in range(1,num_residues+1):
        for res_2 in range(res_1+1, num_residues+1):
            for atom_1, atom_2 in [i for i in combinations_with_replacement(['C','CA','O','N','CB'],2)]:
                if p.residue(res_1).name3() == "GLY" and atom_1 == 'CB':
                    continue
                if p.residue(res_2).name3() == 'GLY' and atom_2 == 'CB':
                    continue
                distance = p.residue(res_1).xyz(atom_1).distance(
                    p.residue(res_2).xyz(atom_2))
                all_distances.append({
                    'atom_1': atom_1,
                    'atom_2': atom_2,
                    'res1':sequence[res_1-1],
                    'res1_3':three_letter_seqs[res_1-1],
                    'p1':res_1,
                    'pdb1':PDB_nums[res_1-1],
                    'chain1':chains[res_1-1],
                    'res2':sequence[res_2-1],
                    'res2_3':three_letter_seqs[res_2-1],
                    'p2':res_2,
                    'pdb2':PDB_nums[res_2-1],
                    'chain2':chains[res_2-1],
                    'dist':distance
                })
    return all_distances
    
def get_distance_df(glob_path):
    distances = []
    for i in glob.glob(glob_path):
        p = pyrosetta.pose_from_pdb(i)
        l = get_all_distances(p)
        df = pandas.DataFrame(l)
        df['INPUT_PDB'] = i
        distances.append(df)
    
    ##If there is more than one
    distances_df = pandas.concat(distances).groupby(
        ['atom_1','atom_2','p1','p2','res1','res2','res1_3','res2_3','pdb1','pdb2','chain1','chain2']).aggregate(
        {'dist':['mean','min','max','std']}).reset_index()
    return distances_df

def get_restraint_df(FIM, distance_df, chain):
    
    ###Get FIM CA DISTANCE Permutations
    all_FIM_distances = permutate_distance_constraints(
                distance_df, DistanceSelector(chain,FIM))
         
    return all_FIM_distances