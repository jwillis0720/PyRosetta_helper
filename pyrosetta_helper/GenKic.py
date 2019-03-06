
from pyrosetta_helper import helper_functions as hf
import pyrosetta
from pyrosetta import (
    init,
    pose_from_file,
    Pose,
    get_fa_scorefxn,
    get_score_function,
    create_score_function,
    pose_from_sequence,
    PyMOLMover)


def get_bb_only_scorefxn():
    backbone_sfxn = create_score_function('empty')
    backbone_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 0.1)
    backbone_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 0.2)
    backbone_sfxn.set_weight(pyrosetta.rosetta.core.scoring.hbond_sr_bb, 2.0)
    backbone_sfxn.set_weight(pyrosetta.rosetta.core.scoring.hbond_lr_bb, 2.0)
    backbone_sfxn.set_weight(pyrosetta.rosetta.core.scoring.rama_prepro, 0.45)
    backbone_sfxn.set_weight(pyrosetta.rosetta.core.scoring.omega, 0.4)
    backbone_sfxn.set_weight(pyrosetta.rosetta.core.scoring.p_aa_pp, 0.6)
    return backbone_sfxn


def connect_disembodied_SSE_with_GK(embodiment_pose_1, embodiment_pose_2, connect_with, trial_name='GKTrial', **kwargs):
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
    connecting_loop = [rtn_residue(hf.get_one_to_three(i))
                       for i in list(connect_with)]

    # Will keep track of indexing of rebuilt loop
    rebuilt_loop = []

    # Get last residue postion on the first pose
    last_residue_on_c_terminus = pose_a.size()
    # First residue to rebuilt
    rebuilt_loop.append(last_residue_on_c_terminus)
    pose_a.set_omega(last_residue_on_c_terminus, 180.1)

    # Iterate through connecting loop and connect it to pose 1
    for resi in connecting_loop:
        pose_a.append_residue_by_bond(resi, True)
        # Add to the index since we added it
        last_residue_on_c_terminus += 1
        rebuilt_loop.append(last_residue_on_c_terminus)
        # And set that omega angle to 180
        pose_a.set_omega(last_residue_on_c_terminus, 180.)

    # Iterate through pose 2 and connect to the C term of the loop we just added
    for residue_index in range(1, pose_b.total_residue()):
        pose_a.append_residue_by_bond(
            pose_b.residue(residue_index))

    # Since we are adding a pose, we don't have to rebuild it with GENKIC. But we should add the Nterm of POSE2
    rebuilt_loop.append(last_residue_on_c_terminus+1)

    print(rebuilt_loop)
    # use GK
    copy_pose = Pose()
    copy_pose.assign(pose_a)
    copy_pose.pdb_info().name(trial_name)

    # Now we just use GenKIC
    # SETUP GK
    gk = pyrosetta.rosetta.protocols.generalized_kinematic_closure.GeneralizedKIC()

    # key words in arguments:
    GK_KEYS = kwargs.keys()
    # We will select the lowest energy loop
    if 'selector_type' in GK_KEYS:
        gk.set_selector_type(kwargs.get('selector_type'))
    else:
        gk.set_selector_type('lowest_energy_selector')

    # Use a backbone only score_function
    if 'selector_scorefunction' in GK_KEYS:
        gk.set_selector_scorefunction(kwargs.get('selector_scorefunction'))
    else:
        gk.set_selector_scorefunction(get_bb_only_scorefxn())

    # Try N times
    if 'closure_attempts' in GK_KEYS:
        gk.set_closure_attempts(kwargs.get('closure_attempts'))
    else:
        gk.set_closure_attempts(100000)

    # Start after N solutions found
    if 'min_solution' in GK_KEYS:
        gk.set_min_solution_count(kwargs.get('min_solution'))
    else:
        gk.set_min_solution_count(1000)

    # We want random perterbations in rama space since we don't know what this loop should look like
    gk.add_perturber('randomize_backbone_by_rama_prepro')

    # Go through rebuilt_loop and add those inddexes
    for res_num in rebuilt_loop:
        gk.add_loop_residue(res_num)
        gk.add_residue_to_perturber_residue_list(res_num)

    # For pivot atoms, lets use first and last and midway
    pivot_atoms = [
        rebuilt_loop[0],
        rebuilt_loop[int(len(rebuilt_loop)/2)],
        rebuilt_loop[-1]]
    gk.set_pivot_atoms(
        rebuilt_loop[0],
        'CA',
        rebuilt_loop[int(len(rebuilt_loop)/2)],
        'CA',
        rebuilt_loop[-1],
        'CA')

    # FILTER 1 - Make sure pivot atoms don't violate rama space
    for p in pivot_atoms:
        gk.add_filter('rama_prepro_check')
        gk.set_filter_resnum(p)
        gk.set_filter_rama_cutoff_energy(2.0)

    # FIlter 2 - Make sure the loop we add does not clash with itself
    gk.add_filter('loop_bump_check')

    # Filter 3 - The two termini residues are off SSE elements, lets enforce that they stay that way
    if 'filter_n_term_ABBA' in GK_KEYS:
        gk.add_filter('backbone_bin')
        gk.set_filter_resnum(rebuilt_loop[0])
        gk.load_filter_bin_params('ABBA')
        gk.set_filter_bin(kwargs.get('filter_n_term_ABBA'))

    # Filter 4 - The two termini residues are off SSE elements, lets enforce that they stay that way
    if 'filter_c_term_ABBA' in GK_KEYS:
        gk.add_filter('backbone_bin')
        gk.set_filter_resnum(rebuilt_loop[-1])
        gk.load_filter_bin_params('ABBA')
        gk.set_filter_bin(kwargs.get('filter_c_term_ABBA'))

    # Since we built from N-C, we need to connect the last residue and n-1 residue
    gk.close_bond(rebuilt_loop[-2], 'C',
                  rebuilt_loop[-1], 'N',
                  # optional params -- use default values
                  rebuilt_loop[-2], 'C', rebuilt_loop[-1], 'N',
                  1.32,
                  114,
                  123,
                  180.,
                  False, False)

    # And finally apply
    gk.apply(copy_pose)
    if gk.get_last_move_status() == pyrosetta.rosetta.protocols.moves.FAIL_RETRY:
        return False
    return copy_pose
