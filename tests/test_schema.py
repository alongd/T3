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
    T3_Options = T3Options(flux_adapter='RMG',
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
    assert T3_Options.flux_adapter == 'RMG'
    assert T3_Options.profiles_adapter == 'RMG'
    assert T3_Options.collision_violator_thermo == False
    assert T3_Options.all_core_species == False
    assert T3_Options.all_core_reactions == False
    assert T3_Options.fit_missing_GAV == False
    assert T3_Options.max_T3_iterations == 5
    assert T3_Options.max_RMG_exceptions_allowed == 5
    assert T3_Options.max_RMG_walltime == '00:02:00:00'
    assert T3_Options.max_T3_walltime == '01:00:00:00'
    assert T3_Options.library_name == 'T3_library'

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
    T3_Sensitivity = T3Sensitivity(adapter=None,
                                   atol=1e-6,
                                   rtol=1e-4,
                                   global_observables=None,
                                   SA_threshold=0.01,
                                   pdep_SA_thershold=0.001,
                                   ME_methods=['CSE', 'MSC'],
                                   top_SA_species=10,
                                   top_SA_reactions=10
                                   )
    assert T3_Sensitivity.adapter == None
    assert T3_Sensitivity.atol == 1e-6
    assert T3_Sensitivity.rtol == 1e-4
    assert T3_Sensitivity.global_observables == None
    assert T3_Sensitivity.SA_threshold == 0.01
    assert T3_Sensitivity.pdep_SA_thershold == 0.001
    assert T3_Sensitivity.ME_methods == ['CSE', 'MSC']
    assert T3_Sensitivity.top_SA_species == 10
    assert T3_Sensitivity.top_SA_reactions == 10

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
    T3_Uncertainty = T3Uncertainty(adapter=None,
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
    assert T3_Uncertainty.adapter == None
    assert T3_Uncertainty.local_analysis == False
    assert T3_Uncertainty.global_analysis == False
    assert T3_Uncertainty.correlated == True
    assert T3_Uncertainty.local_number == 10
    assert T3_Uncertainty.global_number == 5
    assert T3_Uncertainty.termination_time == None
    assert T3_Uncertainty.PCE_run_time == 1800
    assert T3_Uncertainty.PCE_error_tolerance == None
    assert T3_Uncertainty.PCE_max_evals == None
    assert T3_Uncertainty.logx == False

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
