#!/usr/bin/env python3
# encoding: utf-8

"""
schema test module
"""

import pytest
from pydantic import ValidationError

from t3.schema import T3Options, T3Sensitivity, T3Uncertainty

# define a long quote of length 444 characters to test constraints on string length
quote = 'Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et ' \
        'dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ' \
        'ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu ' \
        'fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt ' \
        'mollit anim id est laborum'

def test_T3Options_schema():
    """Test creating an instance of T3Options"""
    # test that a T3Options object can be instantiated properly
    t3_options = T3Options(flux_adapter='RMG',
                           profiles_adapter='RMG',
                           collision_violator_thermo=False,
                           all_core_species=False,
                           all_core_reactions=False,
                           fit_missing_GAV=False,
                           max_T3_iterations=5,
                           max_RMG_exceptions_allowed=5,
                           max_RMG_walltime='00:02:00:00',
                           max_T3_walltime='01:00:00:00',
                           library_name='T3_library'
                           )
    assert t3_options.flux_adapter == 'RMG'
    assert t3_options.profiles_adapter == 'RMG'
    assert t3_options.collision_violator_thermo is False
    assert t3_options.all_core_species is False
    assert t3_options.all_core_reactions is False
    assert t3_options.fit_missing_GAV is False
    assert t3_options.max_T3_iterations == 5
    assert t3_options.max_RMG_exceptions_allowed == 5
    assert t3_options.max_RMG_walltime == '00:02:00:00'
    assert t3_options.max_T3_walltime == '01:00:00:00'
    assert t3_options.library_name == 'T3_library'

    with pytest.raises(ValidationError):
        # check that flux_adapter is constrained to at most 255 characters
        T3Options(flux_adapter=quote)

    with pytest.raises(ValidationError):
        # check that profiles_adapter is constrained to at most 255 characters
        T3Options(profiles_adapter=quote)

    with pytest.raises(ValidationError):
        # check that max_T3_iterations is > 0
        T3Options(max_T3_iterations=0)

    with pytest.raises(ValidationError):
        # check that max_RMG_exceptions_allowed is >= 0
        T3Options(max_RMG_exceptions_allowed=-1)

    with pytest.raises(ValidationError):
        # check that max_RMG_walltime is in the expected format
        T3Options(max_RMG_walltime='05-02:00:00')

    with pytest.raises(ValidationError):
        # check that max_T3_walltime is in the expected format
        T3Options(max_T3_walltime='05-02:00:00')

    with pytest.raises(ValidationError):
        # check that library_name is constrained to at most 255 characters
        T3Options(library_name=quote)


def test_T3Sensitivity_schema():
    """Test creating an instance of T3Sensitivity"""
    # test that a T3Sensitivity object can be instantiated properly
    t3_Sensitivity = T3Sensitivity(adapter=None,
                                   atol=1e-6,
                                   rtol=1e-4,
                                   global_observables=None,
                                   SA_threshold=0.01,
                                   pdep_SA_thershold=0.001,
                                   ME_methods=['CSE', 'MSC'],
                                   top_SA_species=10,
                                   top_SA_reactions=10
                                   )
    assert t3_Sensitivity.adapter is None
    assert t3_Sensitivity.atol == 1e-6
    assert t3_Sensitivity.rtol == 1e-4
    assert t3_Sensitivity.global_observables is None
    assert t3_Sensitivity.SA_threshold == 0.01
    assert t3_Sensitivity.pdep_SA_thershold == 0.001
    assert t3_Sensitivity.ME_methods == ['CSE', 'MSC']
    assert t3_Sensitivity.top_SA_species == 10
    assert t3_Sensitivity.top_SA_reactions == 10

    with pytest.raises(ValidationError):
        # check that adapter is constrained to at most 255 characters
        T3Sensitivity(adapter=quote)

    with pytest.raises(ValidationError):
        # check that atol is constrained to > 0
        T3Sensitivity(atol=0)

    with pytest.raises(ValidationError):
        # check that atol is constrained to < 1e-1
        T3Sensitivity(atol=1)

    with pytest.raises(ValidationError):
        # check that rtol is constrained to > 0
        T3Sensitivity(rtol=0)

    with pytest.raises(ValidationError):
        # check that rtol is constrained to < 1e-1
        T3Sensitivity(rtol=1)

    with pytest.raises(ValidationError):
        # check that entries in global_observables have at least 2 characters
        T3Sensitivity(global_observables='a')

    with pytest.raises(ValidationError):
        # check that entries in global_observables have at most 3 characters
        T3Sensitivity(global_observables='abcd')

    with pytest.raises(ValidationError):
        # check that SA_threshold is constrained to > 0
        T3Sensitivity(SA_threshold=0)

    with pytest.raises(ValidationError):
        # check that SA_threshold is constrained to < 0.5
        T3Sensitivity(SA_threshold=1)

    with pytest.raises(ValidationError):
        # check that pdep_SA_thershold is constrained to > 0
        T3Sensitivity(pdep_SA_thershold=0)

    with pytest.raises(ValidationError):
        # check that pdep_SA_thershold is constrained to < 0.5
        T3Sensitivity(pdep_SA_thershold=1)

    with pytest.raises(ValidationError):
        # check that entries in ME_methods have at least 2 characters
        T3Sensitivity(ME_methods='a')

    with pytest.raises(ValidationError):
        # check that entries in ME_methods have at most 3 characters
        T3Sensitivity(ME_methods='abcd')

    with pytest.raises(ValidationError):
        # check that top_SA_species is constrained to >= 0
        T3Sensitivity(top_SA_species=-1)

    with pytest.raises(ValidationError):
        # check that top_SA_reactions is constrained to >= 0
        T3Sensitivity(top_SA_reactions=-1)


def test_T3Uncertainty():
    """Test creating an instance of T3Uncertainty"""
    # test that a T3Uncertainty object can be instantiated properly
    t3_Uncertainty = T3Uncertainty(adapter=None,
                                   local_analysis=False,
                                   global_analysis=False,
                                   correlated=True,
                                   local_number=10,
                                   global_number=5,
                                   termination_time=None,
                                   PCE_run_time=1800,
                                   PCE_error_tolerance=None,
                                   PCE_max_evals=None,
                                   logx=False
                                   )
    assert t3_Uncertainty.adapter is None
    assert t3_Uncertainty.local_analysis is False
    assert t3_Uncertainty.global_analysis is False
    assert t3_Uncertainty.correlated is True
    assert t3_Uncertainty.local_number == 10
    assert t3_Uncertainty.global_number == 5
    assert t3_Uncertainty.termination_time is None
    assert t3_Uncertainty.PCE_run_time == 1800
    assert t3_Uncertainty.PCE_error_tolerance is None
    assert t3_Uncertainty.PCE_max_evals is None
    assert t3_Uncertainty.logx is False

    with pytest.raises(ValidationError):
        # check that adapter is constrained to at most 255 characters
        T3Uncertainty(adapter=quote)

    with pytest.raises(ValidationError):
        # check that local_number is constrained to > 0
        T3Uncertainty(local_number=0)

    with pytest.raises(ValidationError):
        # check that global_number is constrained to > 0
        T3Uncertainty(global_number=0)

    with pytest.raises(ValidationError):
        # check that termination_time is in the expected format
        T3Uncertainty(termination_time='05-02:00:00')

    with pytest.raises(ValidationError):
        # check that PCE_run_time is constrained to > 0
        T3Uncertainty(PCE_run_time=0)

    with pytest.raises(ValidationError):
        # check that PCE_max_evals is constrained to > 0
        T3Uncertainty(PCE_max_evals=0)
