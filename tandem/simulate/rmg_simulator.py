"""
RMG Simulator Adapter module
Used to run mechanism analysis with RMG
"""

import os
import pandas as pd
import shutil
from typing import Optional, Union

from rmgpy.exceptions import InputError
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.species import Species
from rmgpy.tools.loader import load_rmg_py_job
from rmgpy.tools.simulate import simulate

from tandem.logger import log
from .adapter import SimulateAdapter
from .factory import register_simulate_adapter


class RMGSimulator(SimulateAdapter):
    """
    RMGSimulator is an adapter for the abstract class SimulateAdapter. All output is stored in an "rmg_analysis"
    directory inside the RMG-ARC iteration directory.

    Args:
        run_directory (str): The path to the RMG-ARC iteration directory.
        atol (Optional[float]): The absolute tolerance used when running a simulation.
        rtol (Optional[float]): The relative tolerance used when running a simulation.
        observable_list (Optional[list]): Species used for SA. Entries are dictionaries of 'label' and structure
                                          (either 'smiles' or 'adj'). SA is run only if this argument is specified.
        sa_threshold (Optional[float]): The sensitivity threshold to use.
        sa_atol (Optional[float]): The absolute tolerance used when running SA.
        sa_rtol (Optional[float]): The relative tolerance used when running SA.
        verbose (Optional[bool]): Whether or not to log to file.

    Attributes:
        analysis_directory (str): The path to the "rmg_analysis" directory inside the `run_directory`.
        atol (float): The absolute tolerance used when integrating during an RMG iteration.
        input_file )str): The path to the legacy RMG input file.
        model (str): The path to the annotated chemkin file.
        observable_list (list): Species used for SA. Entries are dictionaries of 'label' and structure
                                (either 'smiles' or 'adj').
        observable_species (list): Species object representations of the species used for SA.
        rmg (RMG class): Representation of an RMG job.
        rmg_input_file (str): The path to the RMG input file that was copied into the "rmg_analysis" directory.
        rtol (float): The relative tolerance used when integrating during an RMG iteration.
        run_directory (str): The path to the RMG-ARC iteration directory.
        sa_atol (float): The absolute tolerance used when running SA.
        sa_rtol (float): The relative tolerance used when running SA.
        sa_path (float): The path to the "rmg_analysis/sa" directory inside the `run_directory`.
        sa_threshold (float): The sensitivity threshold to use.
        species_dict (str): The path to the species dictionary.
        verbose (bool): Whether or not to log to file.
    """

    def __init__(self,
                 run_directory: str,
                 atol: Optional[float] = 1e-6,
                 rtol: Optional[float] = 1e-4,
                 observable_list: Optional[list] = list(),
                 sa_threshold: Optional[float] = 0.001,
                 sa_atol: Optional[float] = 1e-6,
                 sa_rtol: Optional[float] = 1e-4,
                 verbose: Optional[bool] = True,
                 ):

        self.run_directory = run_directory
        self.input_file = os.path.join(os.path.dirname(self.run_directory), 'input.py')
        self.analysis_directory = os.path.join(self.run_directory, 'rmg_analysis')
        self.atol = atol
        self.rtol = rtol
        self.observable_list = observable_list
        self.sa_threshold = sa_threshold
        self.sa_atol = sa_atol
        self.sa_rtol = sa_rtol
        self.verbose = verbose

        # initialize other attributes
        self.model = None
        self.observable_species = list()
        self.rmg = None
        self.rmg_input_file = None
        self.sa_path = None
        self.species_dict = None
        self.spcs = None

        self.set_up()

    def set_up(self):
        """
        Read in the chemkin file, species dictionary, and RMG input file.

        If the user requested SA, simulate the job with SA, generate the output csv files,
        and raise the following errors:
            - ValueError: If a given observable in ``observable_list`` does not have structure (smiles or adjacency list).
            - ValueError: If a given observable ``observable_list`` is not present in the RMG species list.
            - ValueError: If RMG SA is not implemented for the given reactor type.

        If SA is not requested, only simulate the job to obtain species profiles. Output is written to
        a `solver` directory inside `rmg_analysis`.
        """

        self.model = os.path.join(self.run_directory, 'chemkin', 'chem_annotated.inp')
        self.species_dict = os.path.join(self.run_directory, 'chemkin', 'species_dictionary.txt')

        # set up directories
        if not os.path.isdir(self.analysis_directory):
            os.mkdir(self.analysis_directory)
        if len(self.observable_list):
            # store SA results in an SA directory
            self.sa_path = os.path.join(self.analysis_directory, 'sa')
            self.rmg_input_file = os.path.join(self.sa_path, 'input.py')
            if not os.path.isdir(self.sa_path):
                os.mkdir(self.sa_path)
        else:
            # store regular simulation results in solver directory that RMG creates automatically
            self.rmg_input_file = os.path.join(self.analysis_directory, 'input.py')

        # must copy the input file since load_rmg_py_job creates many directories based on the input file's location
        if not os.path.isfile(self.rmg_input_file):
            shutil.copyfile(src=self.input_file, dst=self.rmg_input_file)

        # create rmg object that all methods can access
        self.rmg = load_rmg_py_job(input_file=self.rmg_input_file,
                                   chemkin_file=self.model,
                                   species_dict=self.species_dict,
                                   generate_images=True,
                                   use_chemkin_names=False,
                                   check_duplicates=False,
                                   )
        if self.rmg is None:
            log('The RMG adapter did not properly read the rmg input file', 'error', verbose=True)

        # update the rmg object to perform SA if applicable
        if len(self.observable_list):
            log('Running SA using RMG...', verbose=self.verbose)
            self.spcs = self.rmg.reaction_model.core.species
            for observable in self.observable_list:
                for rmg_spc in self.spcs:
                    if 'adj' in observable:
                        observable_spc = Species(label=observable['label']).from_adjacency_list(observable['adj'])
                    elif 'smiles' in observable:
                        observable_spc = Species(label=observable['label']).from_smiles(observable['smiles'])
                    else:
                        raise ValueError(f'All SA observables must have structure (smiles or adj), got: {observable}')
                    if observable_spc.label == rmg_spc.label or observable_spc.is_isomorphic(rmg_spc):
                        self.observable_species.append(rmg_spc)
                        break
                else:
                    raise ValueError(f'Could not find the observable species {observable["label"]} '
                                     f'in the RMG species list.')

            for reaction_system in self.rmg.reaction_systems:
                if isinstance(reaction_system, SimpleReactor):
                    reaction_system.sensitive_species = self.observable_species
                    reaction_system.sensitivity_threshold = self.sa_threshold
                    if hasattr(reaction_system, 'Trange') and reaction_system.Trange is not None:
                        temperature = sum([t.value_si for t in reaction_system.Trange]) / len(reaction_system.Trange)
                    else:
                        temperature = reaction_system.T.value_si
                    reaction_system.sens_conditions['T'] = temperature
                    if hasattr(reaction_system, 'Prange') and reaction_system.Prange is not None:
                        pressure = sum([p.value_si for p in reaction_system.Prange]) / len(reaction_system.Prange)
                    else:
                        pressure = reaction_system.P.value_si
                    reaction_system.sens_conditions['P'] = pressure
                elif isinstance(reaction_system, LiquidReactor):
                    reaction_system.sensitive_species = self.observable_species
                    reaction_system.sensitivity_threshold = self.sa_threshold
                    if hasattr(reaction_system, 'Trange') and reaction_system.Trange is not None:
                        temperature = sum([t.value_si for t in reaction_system.Trange]) / len(reaction_system.Trange)
                    else:
                        temperature = reaction_system.T.value_si
                    reaction_system.sens_conditions['T'] = temperature
                    if hasattr(reaction_system, 'Vrange') and reaction_system.Vrange is not None:
                        volume = sum([v for v in reaction_system.Vrange]) / len(reaction_system.Vrange)
                    else:
                        volume = reaction_system.V
                    reaction_system.sens_conditions['V'] = volume
                else:
                    raise NotImplementedError(f'RMG SA not implemented for Reactor type {type(reaction_system)}.')

        try:
            log('Simulating a job using RMG...', verbose=self.verbose)
            simulate(self.rmg)
        except FileNotFoundError:
            log('The RMG adapter did not properly simulate the job.', 'error', verbose=True)


    def get_sa_coefficients(self):
        """
        Obtain the SA coefficients.

        Returns:
             A SA dictionary, whose structure is given below.
             sa_dict = { 'thermo' :
                                    { observable_1 :
                                                    { species_1 : 1D array with one entry per time point. Each entry
                                                                  is dLn(observable_1) / dG_species_1 in mol / kcal
                                                                  at the respective time.

                                                     continues for all other species in the model i.e species_2, etc.

                                                    }
                                      observable_2 :  etc...
                                    }

                         'kinetcs' :
                                    { observable_1 :
                                                    { 1     :   1D array with one entry per time point. Each entry
                                                                is dLn(observable_1) / dLn(k_1) at the respective time.

                                                     continues for all other reactions in the model i.e. 2, 3, etc.

                                                    }
                                      observable_2 :  etc...
                                    }

                         'time' :  1D array of time points in seconds.
                        }
        """

        def get_species_by_label(label: str,
                                 species_list: list,
                                 ) -> Union[Species, None]:
            """
            Get a species from a list of species by its label. This label includes the index in parenthesis as assigned
            in the chemkin file. Ex: H(3)

            Args:
                label (str): A species' label (species.label).
                species_list (list): Entries are RMG Species objects.

            Raises:
                InputError: If ``label`` is None.

            Returns:
                Union[Species, None]: The corresponding species from the species_list.
                                      Returns ``None`` if no species was found.
            """
            if label is None:
                raise InputError('Got None as label input')
            for spc in species_list:
                if spc.label == label or spc.to_chemkin() == label:
                    return spc
            if '(' in label and ')' in label:
                # try by the RMG species index
                for spc in species_list:
                    if spc.index == int(label.split('(')[-1].split(')')[0]):
                        return spc
            return None


        solver_path = os.path.join(self.sa_path, 'solver')
        if not os.path.exists(solver_path):
            log("Could not find the path to RMG's solver output folder.", level='error', verbose=self.verbose)
            return None

        sa_files = list()
        for file_ in os.listdir(solver_path):
            if 'sensitivity' in file_ and file_.endswith(".csv"):
                sa_files.append(file_)

        for sa_file in sa_files:
            # iterate through all SA .csv files in the solver folder
            df = pd.read_csv(os.path.join(solver_path, sa_file))
            sa_dict = {'kinetics': dict(), 'thermo': dict(), 'time': list()}
            for header in df.columns:
                # iterate through all headers in the SA .csv file,
                sa_type = None
                if 'Time' in header:
                    sa_dict['time'] = df[header].values
                elif 'dln[k' in header:
                    sa_type = 'kinetics'
                elif 'dG' in header:
                    sa_type = 'thermo'
                if sa_type is not None:
                    # proceed only if we care about this column

                    # check whether the observable requires calculations:
                    observable_label = header.split('[')[1].split(']')[0]
                    observable = get_species_by_label(observable_label, self.spcs)
                    if observable is None:
                        log(f'Could not identify observable species {observable_label}!',
                            level='error', verbose=self.verbose)

                    # continue with the parameter this column represents
                    observable_label = observable.to_chemkin()
                    if observable_label not in sa_dict[sa_type]:
                        sa_dict[sa_type][observable_label] = dict()

                    # parameter extraction examples:
                    # for species get 'C2H4(8)' from `dln[ethane(1)]/dG[C2H4(8)]`
                    # for reaction, get 8 from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)`
                    parameter = header.split('[')[2].split(']')[0]
                    if sa_type == 'kinetics':
                        parameter = int(parameter[1:])

                    sa_dict[sa_type][observable_label][parameter] = df[header].values

        return sa_dict



register_simulate_adapter("RMG", RMGSimulator)
