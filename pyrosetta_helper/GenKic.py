from pyrosetta.rosetta.protocols.generalized_kinematic_closure import GeneralizedKIC

class GenKic():
    
    def __init__(self, loop_residues):
        '''           
        loop_resiues = list of residues involved in GK
       
        '''
        
        self.loop_residues = loop_residues
        self.pivot_residues = [self.loop_residues[0], 
                               self.loop_residues[int(len(self.loop_residues)/2)],
                               self.loop_residues[-1]]

        self.gk_instance = GeneralizedKIC()
        self.selector_type = "lowest_energy_selector"
        self.perturber_type = "randomize_backbone_by_rama_prepro"
        self.closure_attempts = 100000
        self.min_solutions = 1
        self.scorefxn = create_score_function('ref2015')
        
        ##Defaults
        #1. Make sure PIVOT ATAMS ober RAMA
        self.filter_pivot = True
        self.filter_loop_bump = True
        self.dihedral_perturb_value = 3.0
       
    def set_pivot_residues(self,pivot_residues):
        
        self.pivot_residues = pivot_residues
    
    def set_dihedral_pertub_value(self,value):
        self.dihedral_perturb_value = value
    
    def close_normal_bond(self,first_residue,second_residue):
        self.gk_instance.close_bond(
            first_residue, 'C',
            second_residue, 'N',
            first_residue, 'C', 
            second_residue, 'N',
                1.32,
                114,
                123,
                180.,False, False)
        
    def set_selector_type(self,selector_type):
        self.selector_type = selector_type
        
    def set_perturber_type(self,perturber_type):
        self.perturber_type = perturber_type


    def set_closure_attempts(self,closure_attempts):
        self.closure_attempts = closure_attempts
        # Try N times

    def set_stop_after_n(self,n_sol):
        self.min_solutions = n_sol
        
    def set_selector_score_functions(self,scorefxn):
        self.scorefxn = scorefxn
        
    def set_filter_pivot(self,boolean):
        self.filter_pivot = boolean

    def get_instance(self):
        '''Get the instance back with everything set'''
        
        #Selector
        self.gk_instance.set_selector_type(self.selector_type)
        
        #Set perturber
        self.gk_instance.add_perturber(self.perturber_type)
        
        #Set closure attempts
        self.gk_instance.set_closure_attempts(self.closure_attempts)
        
        #Set N solutions
        self.gk_instance.set_min_solution_count(self.min_solutions)
        
        #Set selector scorefunction
        self.gk_instance.set_selector_scorefunction(self.scorefxn)
        
        ##ADD LOOP to PERTURBER AND LoopSet
        for res_num in self.loop_residues:
            self.gk_instance.add_loop_residue(res_num)
            if self.perturber_type == 'randomize_backbone_by_rama_prepro':
                self.gk_instance.add_residue_to_perturber_residue_list(res_num)
                
            if self.perturber_type == 'perturb_dihedral':
                for res_num in self.loop_residues:
                    atomlist_1,atomlist_2 = vector1_core_id_NamedAtomID(), vector1_core_id_NamedAtomID()
                    atomlist_1.append(pyrosetta.rosetta.core.id.NamedAtomID('N',res_num))
                    atomlist_1.append(pyrosetta.rosetta.core.id.NamedAtomID('CA',res_num))
                    atomlist_2.append(pyrosetta.rosetta.core.id.NamedAtomID('CA',res_num))
                    atomlist_2.append(pyrosetta.rosetta.core.id.NamedAtomID('C',res_num))
                    self.gk_instance.add_atomset_to_perturber_atomset_list(atomlist_1)
                    self.gk_instance.add_value_to_perturber_value_list(self.dihedral_perturb_value) 
                    self.gk_instance.add_atomset_to_perturber_atomset_list(atomlist_2)
                    self.gk_instance.add_value_to_perturber_value_list(self.dihedral_perturb_value) 

        ##ADD PIVOT ATOMS
        self.gk_instance.set_pivot_atoms(
            self.pivot_residues[0],'CA',
            self.pivot_residues[1],'CA',
            self.pivot_residues[2],'CA')

        #Rama prepro check on pivot_residues
        if self.filter_pivot:  
            for p in self.pivot_residues:
                self.gk_instance.add_filter('rama_prepro_check')
                self.gk_instance.set_filter_resnum(p)
                self.gk_instance.set_filter_rama_cutoff_energy(2.0)
    
        if self.filter_loop_bump:  
            self.gk_instance.add_filter('loop_bump_check')
        
        
        return self.gk_instance
        

    
def rebuild_with_GenKic(working_pose, native_pose,from_residue=0,to_residue=0,repeats=5):
    '''rebuild whole peptide with GenKic'''    

    # Go through rebuilt_loop and add those inddexes
    rebuilt_loop = []
    if from_residue == 0:
        #Have to start at second residue
        from_residue = 2
    if to_residue == 0:
        to_residue = working_pose.size()-1
    
    for res_num in range(from_residue, to_residue+1):
        rebuilt_loop.append(res_num)
    print(rebuilt_loop)

    pm.keep_history(True)
    
    #Setup GK
    gen_kic_object = GenKic(rebuilt_loop)
    gen_kic_object.set_perturber_type('perturb_dihedral')
    gen_kic_object.set_dihedral_pertub_value = 5.0
    gen_kic_object_instance = gen_kic_object.get_instance()

    
    
    #pyrosetta.rosetta.core.scoring.CA_rmsd(working_pose,)
    scorefxn = create_score_function('ref2015')
    poses = []
    for x in range(1,repeats):
        while True:
            copy_pose = Pose()
            copy_pose.assign(working_pose)
            gen_kic_object_instance.apply(copy_pose)
            rmsd = pyrosetta.rosetta.core.scoring.CA_rmsd(working_pose,copy_pose)
            if rmsd < 2.5:
                fd = fast_design('site1.resfile',working_pose)
                fd.apply(copy_pose)
                scorefxn(copy_pose)
                poses[scorefxn(copy_pose)] = copy_pose
                break
                #print("Yaya!!!")
                #pm.apply(copy_pose)
                #return copy_pose
                
                
    winning_pose = Pose()
    winning_pose.assign(poses[sorted(poses)[0]])        
    return winning_pose    

def connect_disembodied_SSE_with_GK(embodiment_pose_1, 
                                    embodiment_pose_2, 
                                    connect_with,repeats=10,design=False):
    
    pose_a = Pose()
    pose_a.assign(embodiment_pose_1)

    pose_b = Pose()
    pose_b.assign(embodiment_pose_2)

    # Setup CHEMICAL MANAGER TO MAKE NEW RESIDUES
    chm = pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set('fa_standard')

    def rtn_residue(x): return pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
        rts.name_map(x))

    # Residue objects of the connecting loop
    connecting_loop_objects = [rtn_residue(hf.get_one_to_three(i))
                       for i in list(connect_with)]

    # Will keep track of indexing of rebuilt loop
    rebuilt_loop = []

    # Get last residue postion on the first pose
    last_residue_on_c_terminus = pose_a.size()
    
    # First residue to rebuilt
    rebuilt_loop.append(last_residue_on_c_terminus)
    pose_a.set_omega(last_residue_on_c_terminus, 180.1)

    # Iterate through connecting loop and connect it to pose 1
    for resi in connecting_loop_objects:
        pose_a.append_residue_by_bond(resi, True)
        # Add to the index since we added it
        last_residue_on_c_terminus += 1
        rebuilt_loop.append(last_residue_on_c_terminus)
        # And set that omega angle to 180
        pose_a.set_omega(last_residue_on_c_terminus, 180.)

    # Iterate through pose 2 and connect to the C term of the loop we just added
    for residue_index in range(1, pose_b.total_residue()+1):
        pose_a.append_residue_by_bond(
            pose_b.residue(residue_index))

    # Since we are adding a pose, we don't have to rebuild it with GENKIC. 
    #But we should add the Nterm of POSE2
    rebuilt_loop.append(last_residue_on_c_terminus+1)
    #Setup GK
    gen_kic_object = GenKic(rebuilt_loop)
    gen_kic_object.close_normal_bond(rebuilt_loop[-2],rebuilt_loop[-1])
    gk_instance = gen_kic_object.get_instance()
    
    native_pose = Pose()
    native_pose.assign(pose_a)
    
    pm.keep_history(True)
    
    scorefxn = create_score_function('ref2015')
    poses = {}
    for x in range(0,repeats):
        copy_pose = Pose()
        copy_pose.assign(native_pose)
        fill_pose_with_pdb(copy_pose)
        gk_instance.apply(copy_pose)
        
        if design:
            pm.apply(copy_pose)
            with open('/tmp/Temp.resfile','w') as f:
                f.write('NATAA\nEX 1 EX 2\nUSE_INPUT_SC\nstart\n')
                for line in rebuilt_loop[1:-1]:
                    print('designing res {}'.format(line))
                    f.write('{} A NOTAA CM\n'.format(line))
            fd = fast_design('/tmp/Temp.resfile',native=copy_pose,rounds=1)
            fd.apply(copy_pose)
            copy_pose.pdb_info().name('Design_{}'.format(x))
        scorefxn(copy_pose)
        poses[scorefxn(copy_pose)] = copy_pose
        pm.apply(copy_pose)
    winning_pose = Pose()
    winning_pose.assign(poses[sorted(poses)[0]])
    
    return winning_pose
            
    
def relax_mover():
    # Easy Fast Relax Mover
    score_high = pyrosetta.create_score_function('ref2015')
    fr = pyrosetta.rosetta.protocols.relax.FastRelax(score_high, standard_repeats=3)
    fr.constrain_relax_to_start_coords(True)
    return fr
    #fr.apply(p)

def fast_design(resfile,native="",rounds=3):
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These three are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))


    #Score Function
    sf = pyrosetta.create_score_function('ref2015')

    #Fast Relax Mover
    fr = pyrosetta.rosetta.protocols.relax.FastRelax(sf,rounds)

    if native:
        fr.set_native_pose(native)
        fr.constrain_relax_to_native_coords(True)
    
    #Set task factory
    fr.set_task_factory(tf)
    return fr
