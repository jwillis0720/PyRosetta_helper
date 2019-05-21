from pyrosetta.rosetta.protocols.generalized_kinematic_closure import GeneralizedKIC
from pyrosetta import create_score_functions
class GenKic():
    
    def __init__(self, loop_residues):
        """VK Mulligans insanely customizable GenKic algorithm adapted to PyRosetta
        
        Arguments:
            loop_residues {list} -- list of pose residues that will be rebuilt in the GenKic Protocol
        
        it is imperative that you read V. Mulligan tutorials to have any idea what the hell is going on. Currently this is behind
        the RosettaCommons Firwall
        https://github.com/RosettaCommons/demos/blob/master/tutorials/GeneralizedKIC/generalized_kinematic_closure_1.md
        
        """
        ##Make an Instance
        self.gk_instance = GeneralizedKIC()
        
        ##Pose Numbering of residues to consider in the GK. Can be a loop. or Whatever you want actually
        self.loop_residues = loop_residues
        
        for res_num in self.loop_residues:
            self.gk_instance.add_loop_residue(res_num)
        
        #The residues to pivot around GK, Uses, again see part1 of VK tutorial 
        self.pivot_residues = [self.loop_residues[0], 
                               self.loop_residues[int(len(self.loop_residues)/2)],
                               self.loop_residues[-1]]


        ###Here are some basics we can set in the constructor but also have the option to set somewhere else
        self.closure_attempts = 100000
        self.min_solutions = 1
        self.scorefxn = create_score_function('ref2015')
        
        ##Defaults
        #1. Make sure PIVOT atoms obey rama space
        self.filter_pivot = True
        
        #2. Make the loop not bump into itself, i'm not sure why it would ever be off
        self.filter_loop_bump = True
        
        #Every dihedral perterbation will be boltzman dist around this value
        self.dihedral_perturb_value = 3.0
        
    
    def set_scorefxn(self,s):
        self.scorefxn = s
       
    def get_gk_selectors(self):
        return selectors
    
    def get_gk_perturbors(self):
        return perturbers
        
    def set_pivot_residues(self,pivot_residues):
        self.pivot_residues = pivot_residues
    
    def set_dihedral_pertub_value(self,value):
        self.dihedral_perturb_value = value
        
    def set_selector_type(self,selector_type):
        self.selector_type = selector_type
        
    def set_omega_angles(self,omega_value=180):
        print('Setting Omega Angles')
        time.sleep(1)
        self.gk_instance.add_perturber('set_dihedral')
        for res_num in self.loop_residues[:-1]:
            omega_atoms = vector1_core_id_NamedAtomID()
            omega_atoms.append(pyrosetta.rosetta.core.id.NamedAtomID('C',res_num))
            omega_atoms.append(pyrosetta.rosetta.core.id.NamedAtomID('N',res_num+1))
            self.gk_instance.add_atomset_to_perturber_atomset_list(omega_atoms)
        self.gk_instance.add_value_to_perturber_value_list(omega_value)


    def set_closure_attempts(self,closure_attempts):
        self.closure_attempts = closure_attempts
        # Try N times

    def set_stop_after_n(self,n_sol):
        self.min_solutions = n_sol
        
    def set_selector_score_functions(self,scorefxn):
        self.scorefxn = scorefxn
        
    def set_filter_pivot(self,boolean):
        self.filter_pivot = boolean
        
    def close_normal_bond(self,r1,r2):
        self.gk_instance.close_bond(
            r1, 'C',
            r2, 'N',
            r1, 'CA', 
            r2, 'CA',
            1.32829, #ideal bond length
            116.2, #ideal bond angle 1 N-C-CA
            121.7, #ideal bond angle 2 CA-N-C
            180.,
            False, 
            False)
        
    def set_atom_pair_filter(self,r1,r2,a1,a2,d,gt=False):
        self.gk_instance.add_filter('atom_pair_distance')
        self.gk_instance.add_filter_parameter("atom1", a1)
        self.gk_instance.add_filter_parameter("atom2",a2)
        self.gk_instance.add_filter_parameter("res1", r1)
        self.gk_instance.add_filter_parameter("res2", r2)
        self.gk_instance.add_filter_parameter("distance", d)
        if gt:
            self.gk_instance.add_filter_parameter("greater_than", True)
        
    def set_rama_prepro_check(self,r,e=1.0):
        if r not in self.rama_residues:
            self.gk_instance.add_filter('rama_prepro_check')
            self.gk_instance.set_filter_resnum(r)
            self.gk_instance.set_filter_rama_cutoff_energy(e)
            self.rama_residues.append(r)
        else:
            print('Residue {} already being checked by rama prepro filter'.format(r))
        
    def set_backbone_bin(self,r,b,bin_params_file='ABEGO'):
        '''
        Designation that indicates a residue's position in Ramachandran space (A = right-handed alpha or 310 helix; B = right-handed beta strands and extended conformations; E = left-handed beta strands; G = left-handed helices) and cis         omega angles (O). See citation here.
        '''
        if r not in self.rama_bin_residues:
            self.gk_instance.add_filter('backbone_bin')
            self.gk_instance.set_filter_resnum(r)
            self.gk_instance.load_filter_bin_params(bin_params_file)
            self.gk_instance.set_filter_bin(b)
        else:
            print('Residue {} already being checked by backbone bin filter'.format(r))
      
    def add_residue_to_loop(self,r):
       self.gk_instance.add_loop_residue(r)
       
        
    def add_residues_to_perturb_rama(self,r):    
        self.gk_instance.add_perturber('randomize_backbone_by_rama_prepro')
        self.gk_instance.add_residue_to_perturber_residue_list(r)
     
    def add_residue_to_set_backbone_bin(self,r,b):
        
        '''Randomly select mainchain torsions from within a mainchain torsion bin (effect="set_backbone_bin")
This perturber takes a user-specified bin and bin transitions probability file, and randomly chooses mainchain torsions for specified residues from within that bin. The perturber has two additional input options: a bin transitions probability file (bin_params_file="filename.bin_params") and a bin (bin="binname"), where the bin must match a bin named in the bin transitions probability file. Note that the selection of torsion angles from within a bin is based on the sub-bin distribution specified in the bin transitions probability file; it can be uniform or based on the Ramachandran distribution for an alpha-amino acid. (See the file type documentation for details on this.)'''
        self.gk_instance.add_perturber('set_backbone_bin')
        self.gk_instance.add_residue_to_perturber_residue_list(r)
        self.gk_instance.load_perturber_bin_params('ABEGO')
        self.gk_instance.set_perturber_bin(b)
        
    
    
    
    def add_residue_to_perturb_dihedral(self,r):
        self.gk_instance.add_perturber('perturb_dihedral')
        phi_atomlist = vector1_core_id_NamedAtomID()
        psi_atomlist = vector1_core_id_NamedAtomID()    
        ###PHI
        ##C(-1)-N-Cα-C
        phi_atomlist.append(
            pyrosetta.rosetta.core.id.NamedAtomID('N',r))
        phi_atomlist.append(
            pyrosetta.rosetta.core.id.NamedAtomID('CA',r))

        self.gk_instance.add_atomset_to_perturber_atomset_list(phi_atomlist)
        self.gk_instance.add_value_to_perturber_value_list(self.dihedral_perturb_value)

        print("Set Residue {} PHI in Perturber by {}".format(r,self.dihedral_perturb_value))

        ###PSI
        ##N-Cα-C-N(+1)
        psi_atomlist.append(
            pyrosetta.rosetta.core.id.NamedAtomID('CA',r))
        psi_atomlist.append(
            pyrosetta.rosetta.core.id.NamedAtomID('C',r))

        self.gk_instance.add_atomset_to_perturber_atomset_list(psi_atomlist)
        self.gk_instance.add_value_to_perturber_value_list(self.dihedral_perturb_value)
        print("Set Residue {} Psi in Perturber by {}".format(r,self.dihedral_perturb_value))

    
    def get_instance(self):        
        '''Get the instance back with everything set'''
        #Selector
        self.gk_instance.set_selector_type(self.selector_type)
        

        #Set closure attempts
        self.gk_instance.set_closure_attempts(self.closure_attempts)
        
        #Set N solutions
        self.gk_instance.set_min_solution_count(self.min_solutions)
        
        #Set selector scorefunction
        self.gk_instance.set_selector_scorefunction(self.scorefxn)
                        

        ##ADD PIVOT ATOMS
        self.gk_instance.set_pivot_atoms(
            self.pivot_residues[0],'CA',
            self.pivot_residues[1],'CA',
            self.pivot_residues[2],'CA')

        #Rama prepro check on pivot_residues
        if self.filter_pivot:  
            for p in self.pivot_residues:
                self.set_rama_prepro_check(p)
    
        #Self Loop Bump
        if self.filter_loop_bump:  
            self.gk_instance.add_filter('loop_bump_check')
        
        return self.gk_instance