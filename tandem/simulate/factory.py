"""
A module for generating simulate adapters.
"""

from typing import Optional, Type

from .adapter import SimulateAdapter

_registered_simulate_adapters = {}


def register_simulate_adapter(simulator: str,
                              simulate_class: Type[SimulateAdapter],
                              ) -> None:
    """
    A register for the simulate adapters.

    Args:
        simulator: A string representation for a simulate adapter.
        simulate_class: The simulate adapter class.

    Raises:
        TypeError: If ``simulate_class`` is not a ``SimulateAdapter`` instance.
    """
    if not issubclass(simulate_class, SimulateAdapter):
        raise TypeError(f'simulate_class is not a SimulateAdapter, got {simulate_class} which is a {type(simulate_class)}')
    _registered_simulate_adapters[simulator] = simulate_class


def simulate_factory(simulate_method: str,
                     run_directory: str,
                     atol: Optional[float] = 1e-6,
                     rtol: Optional[float] = 1e-4,
                     observable_list: Optional[list] = list(),
                     sa_threshold: Optional[float] = 0.001,
                     sa_atol: Optional[float] = 1e-6,
                     sa_rtol: Optional[float] = 1e-4,
                     verbose: Optional[bool] = True,
                     ) -> Type[SimulateAdapter]:
    """
    A factory generating the simulate adapter corresponding to ``simulate_adapter``.

    Args:
        simulate_method (str): The software to use: RMG, RMS, or Cantera.
        run_directory (str): A path to the RMG-ARC iteration directory.
        atol (Optional[float]): The absolute tolerance used when running a simulation.
        rtol (Optional[float]): The relative tolerance used when running a simulation.
        observable_list (Optional[list]): Species used for SA. Entries are dictionaries of 'label'
                                          and structure (either 'smiles' or 'adj').
        sa_threshold (Optional[float]): The sensitivity threshold to use.
        sa_atol (Optional[float]): The absolute tolerance used when running SA.
        sa_rtol (Optional[float]): The relative tolerance used when running SA.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        Type[SimulateAdapter]: The requested SimulateAdapter child, initialized with the respective arguments.
    """

    if simulate_method.lower() not in ['rmg', 'rms', 'cantera']:
        raise ValueError(f'The "simulate_method" argument must equal "RMG", "RMS", or "Cantera". Got: {simulate_method}'
                         f'Currently only RMG is implemented as a simulate method')  # todo update this message once other adapters are implemented
    simulate_adapter_class = _registered_simulate_adapters[simulate_method](run_directory=run_directory,
                                                                            atol=atol,
                                                                            rtol=rtol,
                                                                            observable_list=observable_list,
                                                                            sa_threshold=sa_threshold,
                                                                            sa_atol=sa_atol,
                                                                            sa_rtol=sa_rtol,
                                                                            verbose=verbose,
                                                                            )

    return simulate_adapter_class
