#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_tandem module
"""

import os
import shutil

from rmgpy import settings as rmg_settings
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.thermo.model import HeatCapacityModel

from arc.common import read_yaml_file

from t3.common import DATA_BASE_PATH, EXAMPLES_BASE_PATH, PROJECTS_BASE_PATH
from t3.main import T3, legalize_species_label


t3_minimal = {'options': {'all_core_reactions': False,
                          'all_core_species': False,
                          'collision_violators_thermo': False,
                          'collision_violators_rates': False,
                          'fit_missing_GAV': False,
                          'flux_adapter': 'RMG',
                          'library_name': 'T3',
                          'max_RMG_exceptions_allowed': 10,
                          'max_RMG_walltime': '00:00:05:00',
                          'max_T3_iterations': 2,
                          'max_T3_walltime': None,
                          'max_rmg_iterations': None,
                          'max_rmg_processes': None,
                          'profiles_adapter': 'RMG'},
              'sensitivity': {'ME_methods': ['CSE', 'MSC'],
                              'SA_threshold': 0.01,
                              'adapter': 'RMG',
                              'atol': 1e-06,
                              'global_observables': None,
                              'pdep_SA_thershold': 0.001,
                              'rtol': 0.0001,
                              'P_list': None,
                              'T_list': None,
                              'top_SA_reactions': 10,
                              'top_SA_species': 10},
              'uncertainty': None}

rmg_minimal = {'database': {'kinetics_depositories': 'default',
                            'kinetics_estimator': 'rate rules',
                            'kinetics_families': 'default',
                            'kinetics_libraries': [],
                            'seed_mechanisms': [],
                            'thermo_libraries': ['primaryThermoLibrary'],
                            'transport_libraries': ['OneDMinN2',
                                                    'PrimaryTransportLibrary',
                                                    'NOx2018',
                                                    'GRI-Mech']},
               'model': {'atol': 1e-16,
                         'branching_index': None,
                         'branching_ratio_max': None,
                         'core_tolerance': [0.01, 0.001],
                         'dynamics_time_scale': None,
                         'filter_reactions': False,
                         'filter_threshold': 100000000.0,
                         'ignore_overall_flux_criterion': None,
                         'max_num_objs_per_iter': 1,
                         'max_num_species': None,
                         'maximum_edge_species': None,
                         'min_core_size_for_prune': None,
                         'min_species_exist_iterations_for_prune': None,
                         'rtol': 1e-08,
                         'terminate_at_max_objects': False,
                         'tolerance_branch_reaction_to_core': None,
                         'tolerance_interrupt_simulation': [0.01, 0.001],
                         'tolerance_keep_in_edge': None,
                         'tolerance_move_edge_reaction_to_core': None,
                         'tolerance_move_edge_reaction_to_core_interrupt': None,
                         'tolerance_move_edge_reaction_to_surface': None,
                         'tolerance_move_edge_reaction_to_surface_interrupt': None,
                         'tolerance_move_surface_reaction_to_core': None,
                         'tolerance_move_surface_species_to_core': None,
                         'tolerance_thermo_keep_species_in_edge': None},
               'options': None,
               'pdep': None,
               'reactors': [{'P': 1.0,
                             'T': 1000.0,
                             'conditions_per_iteration': 12,
                             'termination_conversion': {'H2': 0.9},
                             'termination_rate_ratio': None,
                             'termination_time': 1000000.0,
                             'type': 'gas batch constant T P'}],
               'species': [{'SA_observable': False,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0.67,
                            'constant': False,
                            'inchi': None,
                            'label': 'H2',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[H][H]',
                            'solvent': False},
                           {'SA_observable': False,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0.33,
                            'constant': False,
                            'inchi': None,
                            'label': 'O2',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[O][O]',
                            'solvent': False},
                           {'SA_observable': True,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0,
                            'constant': False,
                            'inchi': None,
                            'label': 'H',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[H]',
                            'solvent': False},
                           {'SA_observable': True,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0,
                            'constant': False,
                            'inchi': None,
                            'label': 'OH',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[OH]',
                            'solvent': False}],
               'species_constraints': None}

qm_minimal = {'adapter': 'ARC',
              'job_types': {'conformers': True,
                            'fine': False,
                            'freq': True,
                            'opt': True,
                            'rotors': False,
                            'sp': True},
              'level_of_theory': 'b3lyp/6-31g(d,p)'}

project_directory = os.path.join(PROJECTS_BASE_PATH, 'test_minimal_delete_after_usage')


def run_minimal() -> T3:
    """A helper function for running the minimal example"""
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    input_dict['verbose'] = 10
    input_dict['project_directory'] = project_directory
    t3 = T3(**input_dict)
    return t3


def test_args_and_attributes():
    """Test passing args and assigning attributes in T3"""
    run_minimal()
    assert os.path.isfile(os.path.join(project_directory, 't3.log'))
    assert not os.path.isdir(os.path.join(project_directory, 'log_archive'))

    t3 = run_minimal()
    assert os.path.isfile(os.path.join(project_directory, 't3.log'))
    assert os.path.isdir(os.path.join(project_directory, 'log_archive'))

    assert t3.project == 'T3_minimal'
    assert t3.project_directory == os.path.join(PROJECTS_BASE_PATH, 'test_minimal_delete_after_usage')
    assert t3.verbose == 10

    assert t3.rmg_exceptions_counter == 0
    assert t3.iteration == 0
    assert t3.thermo_lib_base_path == os.path.join(rmg_settings['database.directory'], 'thermo', 'libraries')
    assert t3.kinetics_lib_base_path == os.path.join(rmg_settings['database.directory'], 'kinetics', 'libraries')
    assert t3.executed_networks == list()
    assert t3.t3 == t3_minimal
    assert t3.rmg == rmg_minimal
    assert t3.qm == qm_minimal
    shutil.rmtree(project_directory)


def test_as_dict():
    """Test T3.as_dict()"""
    t3 = run_minimal()
    assert t3.as_dict() == {'project': 'T3_minimal',
                            'project_directory': project_directory,
                            'qm': qm_minimal,
                            'rmg': rmg_minimal,
                            't3': t3_minimal,
                            'verbose': 10}
    shutil.rmtree(project_directory)


def test_write_input_file():
    """Test automatic input file writing"""
    t3 = run_minimal()
    t3.write_t3_input_file()
    assert os.path.isfile(os.path.join(project_directory, 'auto_saved_input.yml'))
    with open(os.path.join(project_directory, 'auto_saved_input.yml'), 'r') as f:
        assert f.readline() == 'project: T3_minimal\n'
    shutil.rmtree(project_directory)


def test_check_arc_args():
    """Test the check_arc_args() method"""
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    input_dict['verbose'] = 10
    input_dict['project_directory'] = project_directory
    input_dict['qm'] = {'adapter': 'ARC',
                        'unsupported_ARC_arg': 'value',
                        'bac_type': 'm',
                        }
    t3 = T3(**input_dict)
    assert t3.qm['adapter'] == 'ARC'
    assert t3.qm['bac_type'] == 'm'
    assert 'unsupported_ARC_arg' not in t3.qm


def test_run_arc():
    """Test executing ARC"""
    t3 = run_minimal()
    t3.update_paths(iteration=1)
    t3.run_arc(arc_kwargs=t3.qm)
    with open(os.path.join(project_directory, 'iteration_1', 'ARC', 'arc.log'), 'r') as f:
        lines = f.readlines()
        for line in ['Starting project T3\n',
                     'Geometry optimization: b3lyp/6-31g(d,p), software: gaussian (dft)\n',
                     'All jobs terminated. Summary for project T3:\n',
                     'Total execution time: 00:00:00\n',
                     ]:
            assert line in lines
    shutil.rmtree(project_directory)


def test_restart():
    """Test that the restart() method deduces the correct status of a project"""
    restart_base_path = os.path.join(DATA_BASE_PATH, 'restart')

    # empty folders are not saved in git, add them if they don't already exist
    empty_dirs = [os.path.join(restart_base_path, 'r0'),
                  os.path.join(restart_base_path, 'r1', 'iteration_1')]
    for empty_dir in empty_dirs:
        if not os.path.isdir(empty_dir):
            os.makedirs(empty_dir)

    # empty project directory
    # results in iteration=0, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r0'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (0, True)

    # empty 'iteration_1' folder in project directory
    # results in iteration=1, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r1'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (1, True)

    # 'iteration_2' folder with an 'RMG.log' indicating a non-converged job
    # results in iteration=2, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r2'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (2, True)

    # 'iteration_3' folder with an 'RMG.log' indicating a converged job
    # results in iteration=3, run_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r3'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (3, False)

    # 'iteration_4' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # results in iteration=4, run_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r4'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (4, False)

    # 'iteration_5' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # results in iteration=5, run_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r5'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (5, False)

    # 'iteration_6' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # and an ARC 'restart.yml' file
    # results in a complete ARC run, iteration=6+1=7, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r6'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    t3.species = {0: {
                'RMG label': 'Imipramine_1_peroxy',
                'Chemkin label': 'Imipramine_1_peroxy',
                'QM label': 'Imipramine_1_peroxy_0',
                'object': Species(smiles='C'),
                'reasons': ['reason'],
                'converged': None,
                'iteration': 2,
            }}
    t3.dump_species()
    assert t3.restart() == (7, True)
    t3.process_arc_run()
    assert t3.species[0]['converged'] is True
    with open(os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'arc.log'), 'r') as f:
        lines = f.readlines()
        assert 'Starting project T3_ARC_restart_test\n' in lines
        assert 'All jobs terminated. Summary for project T3_ARC_restart_test:\n' in lines

    # 'iteration_7' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a converged job
    # results in iteration=7+1=8, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r7'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (8, True)

    # restore r6 log file
    with open(os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'arc.log'), 'w') as f:
        f.writelines("""Dummy ARC log file\n\n""")

    # delete log files
    for i in range(10):
        directory = os.path.join(restart_base_path, f'r{i}')
        if os.path.isdir(directory):
            log_file = os.path.join(directory, 't3.log')
            if os.path.isfile(log_file):
                os.remove(log_file)
            log_archive = os.path.join(directory, 'log_archive')
            if os.path.isdir(log_archive):
                shutil.rmtree(log_archive)
    shutil.rmtree(os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'output'))
    shutil.rmtree(os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'log_and_restart_archive'))
    os.remove(os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'T3_ARC_restart_test.info'))


def test_legalize_species_label():
    """Test the legalize_species_label() function"""
    species = Species(smiles='C', label='CH4')
    legalize_species_label(species=species)
    assert species.label == 'CH4'

    species = Species(smiles='C#C', label='C#C')
    legalize_species_label(species=species)
    assert species.label == 'C2H2'

    species = Species(smiles='C=CC', label='S(2398)')
    legalize_species_label(species=species)
    assert species.label == 'C3H6'


def test_determine_species_based_on_sa():
    """Test determining species to calculate based on sensitivity analysis"""
    t3 = run_minimal()
    t3.paths['chem annotated'] = os.path.join(DATA_BASE_PATH, 'archive', 'iteration_1', 'chemkin', 'chem_annotated.inp')
    t3.paths['species dict'] = os.path.join(DATA_BASE_PATH, 'archive', 'iteration_1', 'chemkin', 'species_dictionary.txt')
    t3.paths['SA solver'] = os.path.join(DATA_BASE_PATH, 'archive', 'iteration_1', 'sa', 'solver')
    t3.paths['PDep SA'] = os.path.join(DATA_BASE_PATH, 'archive', 'iteration_1', 'PDep_SA')
    t3.qm['species'] = list()
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    species_keys = t3.determine_species_based_on_sa(rmg_species, rmg_reactions)
    assert species_keys == [0, 1, 2, 3, 4, 5, 6]



# spc1 = Species().from_smiles('CC')
# spc1.thermo = HeatCapacityModel()
# spc1.thermo.comment = 'Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH)' \
#                       ' + group(Cs-CsHHH) + radical(RCCJ)'
#
# spc2 = Species().from_smiles('CC')
# spc2.thermo = HeatCapacityModel()
# spc2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(CH3)'
#
# spc3 = Species().from_smiles('CCO')
# spc3.thermo = HeatCapacityModel()
# spc3.thermo.comment = 'Thermo library: primaryThermoLibrary'

# arc_input_file_path = os.path.join(BASE_PATH, 'tandem_1.yml')
# rmg_input_file_path = os.path.join(BASE_PATH, 'rmg', 'input.py')
# 
# spc1 = Species().from_smiles('CC')
# spc1.thermo = HeatCapacityModel()
# spc1.thermo.comment = 'Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) ' \
#                       '+ group(Cs-CsHHH) + radical(RCCJ)'
# 
# spc2 = Species().from_smiles('CC')
# spc2.thermo = HeatCapacityModel()
# spc2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(CH3)'
# 
# spc3 = Species().from_smiles('CCO')
# spc3.thermo = HeatCapacityModel()
# spc3.thermo.comment = 'Thermo library: primaryThermoLibrary'
# 



# def test_run_rmg():
#     """Test the ability to run RMG from T3"""
#     run_rmg_path = os.path.join(BASE_PATH, 'rmg', 'run_rmg')
#     if os.path.isdir(run_rmg_path):
#         shutil.rmtree(run_rmg_path)
#     os.makedirs(run_rmg_path)
#     rmg_input_file_path = os.path.join(run_rmg_path, 'input.py')
#     shutil.copyfile(src=os.path.join(BASE_PATH, 'rmg', 'input.py'), dst=rmg_input_file_path)
#     main.run_rmg(input_file=rmg_input_file_path,
#                  output_directory=run_rmg_path,
#                  kwargs={'restart': '',
#                          'walltime': '00:00:00:00',
#                          'maxproc': 1,
#                          'kineticsdatastore': False},
#                  arguments={'max RMG walltime': '00:00:01:00',
#                             'max RMG exceptions allowed': 0},
#                  tolerance=0.01,
#                  thermo_library=None,
#                  verbose=False)
#     with open(os.path.join(run_rmg_path, ''RMG.log''), 'r') as f:
#         line = f.readline()
#     assert 'RMG execution initiated' in line
#     shutil.rmtree(run_rmg_path)
#
#
# def test_run_arc():
#     """Test the ability to run ARC from T3"""
#     main.run_arc(input_dict={'level_of_theory': 'b3lyp/6-31g**//b3lyp/6-31g**',
#                              'calc_freq_factor': False},
#                  run_directory=BASE_PATH,
#                  species_to_calc=list(),
#                  verbose=False)
#     with open(os.path.join(BASE_PATH, 'ARC', ''arc.log''), 'r') as f:
#         line = f.readline()
#     assert 'ARC execution initiated' in line
#     shutil.rmtree(os.path.join(BASE_PATH, 'ARC'))
#
#
# def test_parse_arc_input_file():
#     """Test parsing an ARC input file and extracting T3 parameters"""
#     arguments, input_dict = main.parse_arc_input_file(os.path.join(BASE_PATH, 'tandem_1.yml'),
#                                                       has_sa=True,
#                                                       has_pdep=False,
#                                                       verbose=False)
#     expected_arguments = {'SA observables': [{'label': 'OH', 'smiles': '[OH]'}],
#                           'SA method': 'RMG',
#                           'SA threshold': 0.001,
#                           'SA species': 10,
#                           'SA reactions': 10,
#                           'SA pdep threshold': 0.1,
#                           'collision violators': True,
#                           'all core species': True,
#                           'RMG tolerances': [0.1, 0.01],
#                           'max tandem iterations': 10,
#                           'max RMG exceptions allowed': 10,
#                           'max RMG walltime': '01:00:00:00'}
#     expected_input_dict = {'level_of_theory': 'b3lyp/6-31g**//b3lyp/6-31g**',
#                            'job_types': {'rotors': False, 'conformers': True, 'fine': False, 'freq': True,
#                                          'opt': True, 'sp': True, 'onedmin': False, 'orbitals': False},
#                            'allow_nonisomorphic_2d': True}
#     assert arguments == expected_arguments
#     assert input_dict == expected_input_dict
#
#
# def test_set_legal_species_labels():
#     """Test setting legal species labels"""
#     # test two species with the same formula
#     species_to_calc = [Species(label='i-C3H7', smiles='C[CH]C'),
#                        Species(label='n-C3H7', smiles='[CH2]CC')]
#     updated_species_to_calc = main.set_legal_species_labels(species_to_calc=species_to_calc, all_species=list())
#     updated_labels = [spc.label for spc in updated_species_to_calc]
#     assert updated_labels == ['C3H7_0', 'C3H7_1']
#
#     # test having a species with the same formula in all_species
#     species_to_calc = [Species(label='C4H9a', smiles='[CH2]CCC'),
#                        Species(label='C4H9b', smiles='C[CH]CC'),
#                        Species(label='NH3', smiles='N')]
#     all_species = [Species(label='C4H9_0', smiles='C[C](C)(C)'),
#                    Species(label='H2O', smiles='O')]
#     updated_species_to_calc = main.set_legal_species_labels(species_to_calc=species_to_calc, all_species=all_species)
#     updated_labels = [spc.label for spc in updated_species_to_calc]
#     assert updated_labels == ['C4H9_1', 'C4H9_2', 'H3N_0']
#
#
# def test_get_species_label_by_structure():
#     """Test getting the species label from a list by its structure"""
#     adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
# 2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
# 3 H u0 p0 c0 {1,S}
# 4 H u0 p0 c0 {1,S}
# 5 H u0 p0 c0 {1,S}
# 6 H u0 p0 c0 {2,S}
# 7 H u0 p0 c0 {2,S}
# 8 H u0 p0 c0 {2,S}"""
#     spc_list = [Species(label='propane').from_smiles('CCC'), Species(label='ethane').from_smiles('CC')]
#     label = main.get_species_label_by_structure(adj=adj, species_list=spc_list)
#     assert label == 'ethane'
#
#
# def test_get_species_by_label():
#     """Test getting a species from a list by its label"""
#     species_to_calc = [Species(label='C4H9a', smiles='[CH2]CCC'),
#                        Species(label='C4H9b', smiles='C[CH]CC'),
#                        Species(label='NH3', smiles='N')]
#     spc = main.get_species_by_label(label='C4H9a', species_list=species_to_calc)
#     print(type(spc))
#     assert isinstance(spc, Species)
#     print(spc.label)
#     assert spc.label == 'C4H9a'
#
#     spc = main.get_species_by_label(label='x', species_list=species_to_calc)
#     print(spc)
#     assert spc is None
#
#
# def test_species_not_in_list():
#     """Test determining whether a species is NOT in a list of species"""
#     # calling set_legal_species_labels() sets the global species_labels_dict
#     species_to_calc = [Species(label='C4H9a', smiles='[CH2]CCC'),
#                        Species(label='C4H9b', smiles='C[CH]CC'),
#                        Species(label='NH3', smiles='N')]
#     updated_species_to_calc = main.set_legal_species_labels(species_to_calc=species_to_calc, all_species=list())
#     is_species_in_list = main.species_not_in_list('C4H9_1', updated_species_to_calc)
#     assert not is_species_in_list
#     is_species_in_list = main.species_not_in_list('C4H9_5', updated_species_to_calc)
#     assert is_species_in_list
#
#
# def test_get_reaction_by_index():
#     """Test getting a reaction from a list by its index"""
#     reaction_list = [Reaction(index=0,
#                               reactants=[Species().from_smiles('[CH2]CC')],
#                               products=[Species().from_smiles('C[CH]C')])]
#     rxn = main.get_reaction_by_index(index=0, reaction_list=reaction_list)
#     assert isinstance(rxn, Reaction)
#     assert rxn.index == 0
#     assert str(rxn) == '[CH2]CC <=> C[CH]C'
#
#     rxn = main.get_reaction_by_index(index=5, reaction_list=reaction_list)
#     assert rxn is None
#
#
# def test_calc_based_on_thermo_comment():
#     """Test which species are selected for calculation based on their thermo comment"""
#     assert main.calc_based_on_thermo_comment(spc1)
#     assert main.calc_based_on_thermo_comment(spc2)
#     assert not main.calc_based_on_thermo_comment(spc3)
#
#
# def test_has_high_uncertainty():
#     """Test determining whether a species thermo should be calculated"""
#     should_species_be_calculated = main.has_high_uncertainty(
#         species=spc1, unconverged_species=list(), species_to_calc=dict())
#     assert should_species_be_calculated
#
#     should_species_be_calculated = main.has_high_uncertainty(
#         species=spc1, unconverged_species=[spc1.copy()], species_to_calc=dict())
#     assert not should_species_be_calculated
#
#     should_species_be_calculated = main.has_high_uncertainty(
#         species=spc1, unconverged_species=list(), species_to_calc={'label': {'spc': spc2}})
#     assert not should_species_be_calculated
#
#     should_species_be_calculated = main.has_high_uncertainty(
#         species=spc1, unconverged_species=list(), species_to_calc={'label': {'spc': spc3}})
#     assert should_species_be_calculated
#
#
# def test_load_species_and_reactions_from_chemkin_file():
#     """Test loading species and reactions from a Chemkin file"""
#     run_directory = os.path.join(BASE_PATH, 'iteration_0')
#     rmg_species, rmg_reactions = main.load_species_and_reactions_from_chemkin_file(run_directory, verbose=False)
#     assert len(rmg_species) == 27
#     assert len(rmg_reactions) == 227
#     assert isinstance(rmg_species[0], Species)
#     assert isinstance(rmg_reactions[0], Reaction)
#
# # def test_determine_species_to_calculate():
# #     """Test that we correctly determine the species to be calculated from an RMG job"""
# #     # TODO: add an actual case with SA and PDep and coll violators
# #     pass
#
#
# def test_determine_species_based_on_sensitivity():
#     """Test determining species to calculate based on sensitivity analysis"""
#     run_directory = os.path.join(BASE_PATH, 'iteration_1')
#     arguments = main.parse_arc_input_file(input_file_path=os.path.join(BASE_PATH, 'tandem_1.yml'),
#                                           has_sa=True,
#                                           has_pdep=False,
#                                           verbose=False)[0]
#     rmg_species, rmg_reactions = main.load_species_and_reactions_from_chemkin_file(run_directory, verbose=False)
#     unconverged_species = list()
#     species_to_calc = main.determine_species_based_on_sensitivity(run_directory,
#                                                                   arguments,
#                                                                   rmg_species,
#                                                                   rmg_reactions,
#                                                                   unconverged_species,
#                                                                   iteration=1,
#                                                                   executed_networks=list(),
#                                                                   verbose=False)[0]
#     assert len(list(species_to_calc.values())) == 7
#     species_to_calc_str = main.dict_to_str(species_to_calc)
#     expected_species_to_calc_str = """ethane(1):
#   spc: ethane(1)
#   reason: observable
# C2H4(17):
#   spc: C=C(17)
#   reason: (iteration 1) participates in the 2nd most sensitive reaction for ethane(1): C=C(17) + H(5) <=> C[CH2](4)
# C2H5(4):
#   spc: C[CH2](4)
#   reason: (iteration 1) participates in the 2nd most sensitive reaction for ethane(1): C=C(17) + H(5) <=> C[CH2](4)
# HO2(8):
#   spc: [O]O(8)
#   reason: (iteration 1) participates in the 5th most sensitive reaction for ethane(1): H(5) + O2(2) <=> [O]O(8)
# C2H4O(41):
#   spc: C=CO(41)
#   reason: (iteration 1) participates in the 6th most sensitive reaction for ethane(1): C=C(17) + O(T)(14) <=> C=CO(41)
# OO(34):
#   spc: OO(34)
#   reason: (iteration 1) participates in the 8th most sensitive reaction for ethane(1): OH(D)(33) + OH(D)(33) <=> OO(34)
# CH3(3):
#   spc: [CH3](3)
#   reason: (iteration 1) the 5th most sensitive species thermo for ethane(1)
# """
#     assert species_to_calc_str == expected_species_to_calc_str
#
#
# #     def test_determine_species_from_pdep_network():
# #         """"""
# # use ests/data/tandem/sa_coefficients.yml with '+' in species names
# #         ch2_adj = """1 C u0 p1 c0 {2,S} {3,S}
# # 2 H u0 p0 c0 {1,S}
# # 3 H u0 p0 c0 {1,S}"""
# #         rmg_species = [Species(label='CO[O](9)').from_smiles('CO[O]'),
# #                        Species(label='[CH2]OO(10)').from_smiles('[CH2]OO'),
# #                        Species(label='O2(2)').from_smiles('[O][O]'),
# #                        Species(label='CH3_0(5)').from_smiles('[CH3]'),
# #                        Species(label='OH(D)(27)').from_smiles('[OH]'),
# #                        Species(label='C=O(26)').from_smiles('C=O'),
# #                        Species(label='[O]CO(28)').from_smiles('[O]CO'),
# #                        Species(label='[O]O(8)').from_smiles('[O]O'),
# #                        Species(label='CH2(S)(3)').from_adjacency_list(ch2_adj),
# #                        Species(label='N2').from_smiles('N#N'),
# #                        Species(label='Ar').from_smiles('[Ar]'),
# #                        Species(label='He').from_smiles('[He]'),
# #                        Species(label='Ne').from_smiles('[Ne]')]
# #         pdep_rxn = PDepReaction(index=0,
# #                                 reactants=[Species(label='[O]O(8)').from_smiles('[O]O'),
# #                                            Species(label='CH2(S)(3)').from_adjacency_list(ch2_adj)],
# #                                 products=[Species(label='CO[O](9)').from_smiles('CO[O]')],
# #                                 network=PDepNetwork(index=27))
# #         pdep_rxns_to_explore = [(pdep_rxn, 1, 'CH2(S)(3)')]
# #         species_to_calc, executed_networks = \
# #             main.determine_species_from_pdep_network(
# #                 run_directory=os.path.join(BASE_PATH, 'iteration_0'),
# #                 pdep_rxns_to_explore=pdep_rxns_to_explore,
# #                 unconverged_species=list(),
# #                 species_to_calc=dict(),
# #                 iteration=1,
# #                 threshold=0.1,
# #                 executed_networks=list(),
# #                 rmg_species=rmg_species,
# #                 verbose=False)
# #         # print(species_to_calc)
# #         # print(executed_networks)
# #         # raise
#
#
# def test_modify_pdep_network_file():
#     """Test modifying an Arkane P-dep network input file"""
#     pdep_sa_path = os.path.join(BASE_PATH, 'iteration_0', 'pdep_sa')
#     if os.path.isdir(pdep_sa_path):
#         shutil.rmtree(pdep_sa_path)
#
#     input_file_path, output_file_path, isomer_labels = \
#         main.modify_pdep_network_file(run_directory=os.path.join(BASE_PATH, 'iteration_0'),
#                                       network_name='network27_2',
#                                       method='CSE')
#     assert 'iteration_0/pdep_sa/network27_2/CSE/input.py' in input_file_path
#     assert 'iteration_0/pdep_sa/network27_2/CSE/sensitivity/sa_coefficients.yml' in output_file_path
#     assert isomer_labels == ('CO[O](9)', '[CH2]OO(10)')
#
#     with open(input_file_path, 'r') as f:
#         lines = f.readlines()
#     sensitivity_conditions, cse = False, False
#     for line in lines:
#         if "    sensitivity_conditions = [[(300, 'K'), (0.01, 'bar')]," in line:
#             sensitivity_conditions = True
#         elif "    method = 'chemically-significant eigenvalues'," in line:
#             cse = True
#     assert sensitivity_conditions
#     assert cse
#
#     input_file_path, output_file_path, isomer_labels = \
#         main.modify_pdep_network_file(run_directory=os.path.join(BASE_PATH, 'iteration_0'),
#                                       network_name='network3_1',
#                                       method='MSC')
#     assert 'iteration_0/pdep_sa/network3_1/MSC/input.py' in input_file_path
#     assert 'iteration_0/pdep_sa/network3_1/MSC/sensitivity/sa_coefficients.yml' in output_file_path
#     assert isomer_labels == ('ethane(1)',)
#
#     with open(input_file_path, 'r') as f:
#         lines = f.readlines()
#     sensitivity_conditions, msc, wrong_msc = False, False, False
#     for line in lines:
#         if "    sensitivity_conditions = [[(300, 'K'), (0.01, 'bar')]," in line:
#             sensitivity_conditions = True
#         elif "    method = 'modified strong collision'," in line:
#             msc = True
#         elif "modified strong collision" in line:
#             # this string must not appear in any other line
#             wrong_msc = True
#     assert sensitivity_conditions
#     assert msc
#     assert not wrong_msc
#
#     shutil.rmtree(pdep_sa_path)
#
#
# def test_determine_species_based_on_collision_violators():
#     """Test determining species to calculate based on collision rate violating reactions"""
#     run_directory = os.path.join(BASE_PATH, 'iteration_2')
#     rmg_species, rmg_reactions = main.load_species_and_reactions_from_chemkin_file(run_directory, verbose=False)
#     unconverged_species = list()
#     species_to_calc = main.determine_species_based_on_collision_violators(run_directory,
#                                                                           rmg_species,
#                                                                           unconverged_species,
#                                                                           verbose=False)
#     assert len(list(species_to_calc.values())) == 1
#     species_to_calc_str = main.dict_to_str(species_to_calc)
#     expected_species_to_calc_str = """C3H6(19):
#   spc: [CH2]C[CH2](19)
#   reason: species participates in a collision rate violating reaction, C3H6(19)+H(4)=C3H7_0(15)
# """
#     assert species_to_calc_str == expected_species_to_calc_str
#
#
# def test_add_rmg_libraries():
#     """Test adding an RMG library to the RMG database repository"""
#     rmg_thermo_lib_1_path = os.path.join(settings['database.directory'], 'thermo', 'libraries', 't3_thermo.py')
#     rmg_thermo_lib_2_path = os.path.join(settings['database.directory'], 'thermo', 'libraries', 't3_thermo_0.py')
#     libraries_path = os.path.join(BASE_PATH, 'iteration_0')
#     local_context = {'ThermoData': ThermoData, 'Wilhoit': Wilhoit, 'NASAPolynomial': NASAPolynomial, 'NASA': NASA}
#
#     if os.path.isfile(rmg_thermo_lib_1_path):
#         os.remove(rmg_thermo_lib_1_path)
#     if os.path.isfile(rmg_thermo_lib_2_path):
#         os.remove(rmg_thermo_lib_2_path)
#
#     # test adding a library for the fist time
#     library_name = main.add_rmg_libraries(run_directory=libraries_path, library_name=None, verbose=False)
#     assert library_name == 't3_thermo'
#     thermo_lib = ThermoLibrary()
#     thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
#     assert len(list(thermo_lib.entries.values())) == 10
#
#     # test adding a library for the fist time with a different name ('t3_thermo' is occupied)
#     library_name = main.add_rmg_libraries(run_directory=libraries_path, library_name=None, verbose=False)
#     assert library_name == 't3_thermo_0'
#     thermo_lib = ThermoLibrary()
#     thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
#     assert len(list(thermo_lib.entries.values())) == 10
#
#     # test appending entries to an existing library
#     libraries_path = os.path.join(BASE_PATH, 'iteration_1')
#     library_name = main.add_rmg_libraries(run_directory=libraries_path, library_name='t3_thermo', verbose=False)
#     assert library_name == 't3_thermo'
#     thermo_lib = ThermoLibrary()
#     thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
#     assert len(list(thermo_lib.entries.values())) == 11  # extended with one additional entry
#
#     os.remove(rmg_thermo_lib_1_path)
#     os.remove(rmg_thermo_lib_2_path)
#
#
# def test_get_unconverged_species():
#     """Test attaining a list of unconverged species from the ARC project info file"""
#     labels = ['C3H8_0', 'C3H7_1', 'C3H6_0', 'C3H6_1', 'C3H5_1',
#               'C3H5_2', 'C3H6_2', 'C3H4_0', 'C3H3_0', 'C3H4_1', 'C3H4_2']
#     all_species = [Species(label=label) for label in labels]
#     unconverged_species = main.get_unconverged_species(run_directory=os.path.join(BASE_PATH, 'iteration_0'),
#                                                        all_species=all_species,
#                                                        log_species=False,
#                                                        verbose=False)
#     assert len(unconverged_species) == 1
#     assert unconverged_species[0].label == 'C3H4_1'
