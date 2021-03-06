#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_tandem module
"""

import os
import shutil
from typing import Optional

from rmgpy import settings as rmg_settings
from rmgpy.rmg.pdep import PDepNetwork, PDepReaction
from rmgpy.species import Species
from rmgpy.thermo import NASA

from arc.common import read_yaml_file

from t3.common import DATA_BASE_PATH, EXAMPLES_BASE_PATH, PROJECTS_BASE_PATH
from t3.main import (T3,
                     legalize_species_label,
                     get_species_by_label,
                     get_reaction_by_index,
                     get_species_label_by_structure)
from t3.utils.writer import write_rmg_input_file


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
                              'pdep_SA_threshold': 0.001,
                              'rtol': 0.0001,
                              'P_list': None,
                              'T_list': None,
                              'top_SA_reactions': 10,
                              'top_SA_species': 10},
              'uncertainty': None,
              }

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
               'species_constraints': None,
               }

qm_minimal = {'adapter': 'ARC',
              'job_types': {'conformers': True,
                            'fine': False,
                            'freq': True,
                            'opt': True,
                            'rotors': False,
                            'sp': True},
              'level_of_theory': 'b3lyp/6-31g(d,p)',
              'reactions': [],
              'species': [],
              }

test_minimal_project_directory = os.path.join(PROJECTS_BASE_PATH, 'test_minimal_delete_after_usage')
restart_base_path = os.path.join(DATA_BASE_PATH, 'restart')
dump_species_path = os.path.join(DATA_BASE_PATH, 'test_dump_species')


def setup_module():
    """
    Setup.
    Useful for rerunning these tests after a failed test during development.
    """
    if os.path.isdir(test_minimal_project_directory):
        shutil.rmtree(test_minimal_project_directory)


def run_minimal(project: Optional[str] = None,
                project_directory: Optional[str] = None,
                iteration: Optional[int] = None,
                set_paths: bool = False,
                ) -> T3:
    """A helper function for running the minimal example"""
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    input_dict['verbose'] = 10
    input_dict['project_directory'] = project_directory or test_minimal_project_directory
    if project is not None:
        input_dict['project'] = project
    t3 = T3(**input_dict)
    t3.iteration = iteration or 0
    if set_paths:
        t3.set_paths()
    return t3


def test_args_and_attributes():
    """Test passing args and assigning attributes in T3"""
    run_minimal()
    assert os.path.isfile(os.path.join(test_minimal_project_directory, 't3.log'))
    assert not os.path.isdir(os.path.join(test_minimal_project_directory, 'log_archive'))

    t3 = run_minimal()
    assert os.path.isfile(os.path.join(test_minimal_project_directory, 't3.log'))
    assert os.path.isdir(os.path.join(test_minimal_project_directory, 'log_archive'))

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
    shutil.rmtree(test_minimal_project_directory)


def test_as_dict():
    """Test T3.as_dict()"""
    t3 = run_minimal()
    assert t3.as_dict() == {'project': 'T3_minimal',
                            'project_directory': test_minimal_project_directory,
                            'qm': qm_minimal,
                            'rmg': rmg_minimal,
                            't3': t3_minimal,
                            'verbose': 10}
    shutil.rmtree(test_minimal_project_directory)


def test_write_t3_input_file():
    """Test automatically writing a T3 input file"""
    t3 = run_minimal()
    t3.write_t3_input_file()
    assert os.path.isfile(os.path.join(test_minimal_project_directory, 'auto_saved_input.yml'))
    with open(os.path.join(test_minimal_project_directory, 'auto_saved_input.yml'), 'r') as f:
        assert f.readline() == 'project: T3_minimal\n'
    shutil.rmtree(test_minimal_project_directory)


def test_set_paths():
    """Test updating self.paths"""
    t3 = run_minimal(iteration=1, set_paths=True)
    paths = {'ARC': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC',
             'ARC info': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/T3.info',
             'ARC input': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/input.yml',
             'ARC kinetics lib': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/output/RMG '
                                 'libraries/kinetics',
             'ARC log': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/arc.log',
             'ARC restart': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/restart.yml',
             'ARC thermo lib': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/output/RMG '
                               'libraries/thermo/T3.py',
             'PDep SA': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/PDep_SA',
             'RMG': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG',
             'RMG PDep': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/pdep',
             'RMG coll vio': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/collision_rate_violators.log',
             'RMG input': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/input.py',
             'RMG log': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/RMG.log',
             'SA': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA',
             'SA input': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/input.py',
             'SA solver': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/solver',
             'chem annotated': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/chemkin/chem_annotated.inp',
             'iteration': 'T3/Projects/test_minimal_delete_after_usage/iteration_1',
             'species dict': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/chemkin/'
                             'species_dictionary.txt'}
    for key, path in t3.paths.items():
        assert paths[key] in path


def test_restart():
    """Test that the restart() method deduces the correct status of a project"""
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


def test_check_arc_args():
    """Test the check_arc_args() method"""
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    input_dict['verbose'] = 10
    input_dict['project_directory'] = test_minimal_project_directory
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
    t3 = run_minimal(iteration=1, set_paths=True)
    t3.run_arc(arc_kwargs=t3.qm)
    with open(t3.paths['ARC log'], 'r') as f:
        lines = f.readlines()
    for line in ['Starting project T3\n',
                 'Geometry optimization: b3lyp/6-31g(d,p), software: gaussian (dft)\n',
                 'All jobs terminated. Summary for project T3:\n',
                 'Total execution time: 00:00:00\n',
                 ]:
        assert line in lines
    shutil.rmtree(test_minimal_project_directory)


def test_process_arc_run():
    """Tests processing an ARC run and copying over a thermo library to the RMG-database repository"""
    t3 = run_minimal(project='T3',
                     project_directory=os.path.join(DATA_BASE_PATH, 'process_arc'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.species = {0: {'RMG label': 'imipramine_ol_2_ket_4',
                      'Chemkin label': 'imipramine_ol_2_ket_4',
                      'QM label': 'imipramine_ol_2_ket_4_0',
                      'object': Species(smiles='C'),
                      'reasons': ['reason 1', 'reason 2'],
                      'converged': None,
                      'iteration': 1},
                  1: {'RMG label': 'imipramine_ol_2_ket_5',
                      'Chemkin label': 'imipramine_ol_2_ket_5',
                      'QM label': 'imipramine_ol_2_ket_5_1',
                      'object': Species(smiles='CC'),
                      'reasons': ['reason 3'],
                      'converged': None,
                      'iteration': 1},
                  }
    t3.process_arc_run()
    assert t3.species[0]['converged'] is True
    assert t3.species[1]['converged'] is False
    thermo_lib_path = os.path.join(t3.thermo_lib_base_path, 'T3.py')
    assert os.path.isfile(thermo_lib_path)
    with open(thermo_lib_path, 'r') as f:
        lines = f.readlines()
    for line in ['name = "T3"\n',
                 "Species imipramine_ol_2_ket_4 (run time: 1 day, 8:24:38)\n",
                 '    label = "imipramine_ol_2_ket_4",\n',
                 "        E0 = (-171.078,'kJ/mol'),\n"]:
        assert line in lines
    os.remove(thermo_lib_path)


def test_get_current_rmg_tol():
    """Test getting the correct RMG tolerances"""
    t3 = run_minimal()
    t3.rmg['model']['core_tolerance'] = [0.1, 0.05, 0.01, 0.001]
    t3.iteration = 1
    assert t3.get_current_rmg_tol() == 0.1
    t3.iteration = 2
    assert t3.get_current_rmg_tol() == 0.05
    t3.iteration = 3
    assert t3.get_current_rmg_tol() == 0.01
    t3.iteration = 4
    assert t3.get_current_rmg_tol() == 0.001
    t3.iteration = 5
    assert t3.get_current_rmg_tol() == 0.001
    t3.iteration = 6
    assert t3.get_current_rmg_tol() == 0.001
    t3.iteration = 238
    assert t3.get_current_rmg_tol() == 0.001


def test_run_rmg():
    """Test the ability to run RMG from T3"""
    t3 = run_minimal(iteration=1, set_paths=True)
    write_rmg_input_file(
        kwargs=t3.rmg,
        iteration=t3.iteration,
        path=t3.paths['RMG input'],
        walltime=t3.t3['options']['max_RMG_walltime'],
    )
    t3.run_rmg()
    with open(t3.paths['RMG input'], 'r') as f:
        lines = f.readlines()
    for line in ["    thermoLibraries=['primaryThermoLibrary'],\n",
                 "simulator(atol=1e-16, rtol=1e-08)\n",
                 ]:
        assert line in lines
    with open(t3.paths['RMG log'], 'r') as f:
        lines = f.readlines()
    for line in ["    thermoLibraries=['primaryThermoLibrary'],\n",
                 "simulator(atol=1e-16, rtol=1e-08)\n",
                 "No collision rate violators found in the model's core.\n",
                 "MODEL GENERATION COMPLETED\n",
                 "The final model core has 12 species and 18 reactions\n",
                 ]:
        assert line in lines
    assert os.path.isfile(t3.paths['chem annotated'])
    assert os.path.isfile(t3.paths['species dict'])
    shutil.rmtree(test_minimal_project_directory)


def test_determine_species_to_calculate():
    """Test determining the species to be calculated"""

    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_species'))

    # 1. no calculations required
    t3.iteration = 1
    t3.set_paths()
    t3.t3['options']['all_core_species'] = True
    additional_calcs_required = t3.determine_species_to_calculate()
    assert not additional_calcs_required

    # 2. All core species
    t3.iteration = 2
    t3.set_paths()
    t3.t3['options']['all_core_species'] = True
    additional_calcs_required = t3.determine_species_to_calculate()
    assert additional_calcs_required
    assert len(list(t3.species.keys())) == 3
    assert all([species_dict['reasons'] == ['All core species'] for species_dict in t3.species.values()])
    assert all([species_dict['RMG label'] in ['OH', 'HO2', 'H2O2'] for species_dict in t3.species.values()])
    assert all([species_dict['QM label'] in ['OH_0', 'HO2_1', 'H2O2_2'] for species_dict in t3.species.values()])

    # 3. collision violators
    t3.iteration = 3
    t3.set_paths()
    t3.species = dict()
    t3.t3['options']['all_core_species'] = False
    t3.t3['options']['collision_violators_thermo'] = True
    additional_calcs_required = t3.determine_species_to_calculate()
    assert additional_calcs_required
    assert len(list(t3.species.keys())) == 20
    assert all(['Species participates in collision rate violating reaction:' in species_dict['reasons'][0]
                for species_dict in t3.species.values() if species_dict['RMG label'] not in ['H', 'OH']])

    # 4. SA observables
    assert t3.species[0]['RMG label'] == 'H'
    assert t3.species[0]['reasons'] == ['SA observable']
    assert t3.species[1]['RMG label'] == 'OH'
    assert t3.species[1]['reasons'] == ['SA observable']


def test_determine_species_based_on_sa():
    """Test determining species to calculate based on sensitivity analysis"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    species_keys = t3.determine_species_based_on_sa()
    assert species_keys == [0, 1]


def test_determine_species_from_pdep_network():
    """Test determining species from pdep network"""
    t3 = run_minimal(project_directory = os.path.join(DATA_BASE_PATH, 'pdep_network'),
                     iteration=1,
                     set_paths=True,
                     )
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    # focus on reaction 2 in network4_2.py, whose species correspond to the indices below
    # reactants = ['H(34)', 'C4ene(26)'],
    # products = ['C4rad(5)'],
    pdep_rxn = PDepReaction(index=1,
                            reactants=[t3.rmg_species[35],
                                       t3.rmg_species[27]],
                            products=[t3.rmg_species[6]],
                            network=PDepNetwork(index=4))
    pdep_rxns_to_explore = [(pdep_rxn, 2, t3.rmg_species[6].label)]
    species_keys = t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
    assert len(species_keys) == 1
    shutil.rmtree(t3.paths['PDep SA'])


def test_determine_species_based_on_collision_violators():
    """Test determining species to calculate based on collision rate violating reactions"""
    t3 = run_minimal()
    t3.paths['RMG coll vio'] = os.path.join(DATA_BASE_PATH, 'collision_rate_violators', 'collision_rate_violators.log')
    t3.paths['chem annotated'] = os.path.join(DATA_BASE_PATH, 'collision_rate_violators', 'chem_annotated.inp')
    t3.paths['species dict'] = os.path.join(DATA_BASE_PATH, 'collision_rate_violators', 'species_dictionary.txt')
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    species_to_calc = t3.determine_species_based_on_collision_violators()
    assert len(species_to_calc) == 18
    expected_species_to_calc = [
        'C7H13(920)',
        'C6H9(1933)',
        'C6H8(2025)',
        'C6H9(2035)',
        'C6H8(2027)',
        'S(1752)',
        'S(11767)',
        'S(11972)',
        'S(17233)',
        'S(16488)',
        'S(16530)',
        'S(25139)',
        'S(16448)',
        'S(16972)',
        'S(13229)',
        'C6H8(8657)',
        'S(26357)',
        'S(25149)'
    ]
    for index in species_to_calc:
        assert t3.species[index]['Chemkin label'] == expected_species_to_calc[index]


def test_trsh_rmg_tol():
    """Test troubleshooting the RMG tolerance"""
    t3 = run_minimal()
    t3.t3['options']['max_T3_iterations'] = 10

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 1
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.05, 0.001, 0.0001]

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 2
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.1, 0.001, 0.0001]

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 6
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.1, 0.001, 0.0001, 0.00005, 0.00005, 0.00005]

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 12
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.1, 0.001, 0.0001]


def test_species_requires_refinement():
    """Test whether a species thermo requires refinement"""
    t3 = run_minimal()
    spc_1 = Species(smiles='C')
    spc_1.thermo = NASA()
    spc_1.thermo.comment = 'Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-CO) + ' \
                           'missing(O2d-CO) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + ' \
                           'group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)'
    spc_2 = Species(smiles='C')
    spc_2.thermo = NASA()
    spc_2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(HOOj)'
    spc_3 = Species(smiles='C')
    spc_3.thermo = NASA()
    spc_3.thermo.comment = 'Thermo library: primaryThermoLibrary'
    assert t3.species_requires_refinement(spc_1) is True
    assert t3.species_requires_refinement(spc_2) is True
    assert t3.species_requires_refinement(spc_3) is False


def test_get_species_key():
    """Test checking whether a species already exists in self.species and getting its key"""
    t3 = run_minimal(project_directory = os.path.join(DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.t3['options']['all_core_species'] = True
    t3.determine_species_to_calculate()

    # 1. by species
    assert t3.get_species_key(species=Species(smiles='[OH]')) == 0
    assert t3.get_species_key(species=Species(smiles='O[O]')) == 1
    assert t3.get_species_key(species=Species(smiles='OO')) == 2
    assert t3.get_species_key(species=Species(smiles='O')) is None

    # 2. by label
    t3.species = {5: {'QM label': 'O2'}}
    key = t3.get_species_key(label='O2')
    assert key == 5


def test_load_species_and_reactions_from_chemkin_file():
    """Test loading RMG species and reactions from a Chemkin file"""
    t3 = run_minimal(project_directory = os.path.join(DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    assert len(rmg_species) == 12
    assert len(rmg_reactions) == 18
    assert rmg_species[0].label == 'Ar'
    assert rmg_species[10].label == 'H2O'
    assert str(rmg_reactions[0]) == 'H(3) + H(3) <=> H2(1)'
    assert str(rmg_reactions[10]) == 'OH(4) + H2O2(9) <=> HO2(6) + H2O(7)'


def test_add_species():
    """Test adding a species to self.species and to self.qm['species']"""
    t3 = run_minimal(project_directory = os.path.join(DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.t3['options']['all_core_species'] = True
    t3.determine_species_to_calculate()
    spc_1 = Species(label='OH', smiles='[OH]')
    spc_2 = Species(label='hydrazine', smiles='NN')

    assert t3.get_species_key(species=spc_1) == 0
    assert t3.species[0]['RMG label'] == 'OH'
    assert t3.species[0]['reasons'] == ['All core species']

    t3.add_species(species=spc_1, reasons='Some other reason')
    assert t3.get_species_key(species=spc_1) == 0
    assert t3.species[0]['RMG label'] == 'OH'
    assert t3.species[0]['reasons'] == ['All core species', 'Some other reason']

    assert t3.get_species_key(species=spc_2) is None

    t3.add_species(species=spc_2, reasons=['R1', 'R2'])
    assert t3.get_species_key(species=spc_2) == 3
    assert t3.species[3]['RMG label'] == 'hydrazine'
    assert t3.species[3]['reasons'] == ['R1', 'R2']


def test_dump_species():
    """Test dump species for restart purposes"""
    # create an empty `iteration_5` directory
    if not os.path.isdir(os.path.join(dump_species_path, 'iteration_5')):
        os.makedirs(os.path.join(dump_species_path, 'iteration_5'))
    t3 = T3(project='test_dump_species',
            project_directory=dump_species_path,
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
    assert os.path.isfile(os.path.join(dump_species_path, 't3.log'))
    assert os.path.isfile(os.path.join(dump_species_path, 'species.yml'))
    assert t3.restart() == (5, True)


def test_load_species():
    """Test loading the dumped species dictionary from `test_dump_species()` above"""
    t3 = T3(project='test_dump_species',
            project_directory=dump_species_path,
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    t3.load_species()
    assert t3.species[0]['Chemkin label'] == 'Imipramine_1_peroxy'
    assert t3.species[0]['QM label'] == 'Imipramine_1_peroxy_0'


# main functions:


def test_get_species_by_label():
    """Test getting species by label"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    label = 'H2O'
    species = get_species_by_label(label, rmg_species)
    assert species.label == label
    assert species.index == 7


def test_get_reaction_by_index():
    """Test getting reaction by index"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    index = 5
    reaction = get_reaction_by_index(index, rmg_reactions)
    assert reaction.reactants[0].label == 'H'
    assert reaction.reactants[1].label == '[O]O'
    assert reaction.products[0].label == 'OO'


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


def test_get_species_label_by_structure():
    """Test getting the species label from a list by its structure"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    adj = """1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}"""
    label = get_species_label_by_structure(adj, rmg_species)
    assert label == 'H2O'


def teardown_module():
    """teardown any state that was previously setup with a setup_module method."""
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
    # delete project folders
    for directory in [test_minimal_project_directory,
                      dump_species_path,
                      os.path.join(DATA_BASE_PATH, 'minimal_data', 'log_archive')]:
        if os.path.isdir(directory):
            shutil.rmtree(directory)
