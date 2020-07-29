import h5py
from typing import Union
from warnings import warn
from .misc import *


class ResultExtractor:

    def __init__(self,
                 result_file: Union[str, None] = None,
                 simulation_file: Union[str, None] = None,
                 validation_file: Union[str, None] = None,
                 optimization_type: Union[str, None] = None):
        """
        Args:
            result_file:
                Path to h5 file with parPE optimization results
            simulation_file:
                Path to h5 file with parPE simulation results
            validation_file:
                Path to h5 file with parPE validation results
            optimization_type:
                String indicating 'fullbatch' or 'minibatch' optimization type.
                 If a result_file is provided, this will be taken from the
                 settings of the result_file, potentially overriding this input
        """

        # Get final llhs
        self.result_file = result_file
        if result_file is not None:
            self.optimization_type, self.optimizer, self.n_starts, \
                self.final_parameters, self.final_cost_values, \
                self.final_cpu_times = \
                self.process_result_file(self.result_file)

            if optimization_type is not None and \
                    optimization_type != self.optimization_type:
                warn('Optimization type was passed to ResultExtractor but not '
                     'coinciding with the settings from the parPE optimization '
                     'results. Overriding user-defined optimization_type!')

        else:
            if optimization_type is not None:
                self.optimization_type = optimization_type

            # We have no result_file, so we can't set all values...
            self.optimizer = self.n_starts = self.final_parameters = \
                self.final_cost_values = self.final_cpu_times = None

        self.simulation_file = simulation_file
        if simulation_file is not None:
            self.n_starts, self.llhs, self.simulations, self.correlations,\
                self.cpu_times = \
                self.process_simulation_file(self.simulation_file)

        self.validation_file = validation_file
        if validation_file is not None:
            self.process_validation_file(self.validation_file)

    @staticmethod
    def process_result_file(result_file):
        # create return values
        finalCosts = finalParameters = finalCpuTimes = []

        with h5py.File(result_file, 'r') as f:
            # first determine optimizer and optimization type
            optimizer = f['inputData/optimizationOptions'].attrs['optimizer']
            if optimizer == 10:
                optimizationType = 'minibatch'
            else:
                optimizationType = 'fullbatch'

            try:
                # first get the number of starts, then try to iterate over them
                # and extract information
                nStarts = len(f['/multistarts'])
                for iStart in range(nStarts):
                    current_start = f[f'/multistarts/{iStart}/']

                    # extract cpu time, final cost and parameters
                    if 'finalCost' in current_start:
                        finalCosts.append(current_start['finalCost'])
                    if 'finalParameters' in current_start:
                        finalCosts.append(current_start['finalCost'])
                    if 'cpuSec' in current_start:
                        finalCosts.append(current_start['finalCpuTimes'])

            except KeyError:
                nStarts = 0
                warn('No multistarts were found in result_file. Stopping!')

        return optimizationType, optimizer, nStarts, finalCosts, \
               finalParameters, finalCpuTimes

    @staticmethod
    def process_simulation_file(simulation_file):
        """

        Args:
            simulation_file:
                path to simulation file
        Returns:

        """
        #
        try:
            meaurements, simulations, cpuTimes, llhs = \
                readSimulationsFromFile(simulation_file)
            correlations = getCorrTable(meaurements, simulations)
        except:
            meaurements, simulations, llhs = \
                readSimulationsFromFileLegacy(simulation_file)
            cpuTimes = None
            correlations = getCorrTableLegacy([meaurements, simulations])

        # Get the number of multistarts. Needs to be compared to the one from the result file
        with h5py.File(simulation_file, 'r') as f:
            nStarts = len(f['/multistarts'])

        return nStarts, llhs, simulations, correlations, cpuTimes

    @staticmethod
    def process_validation_file(validation_file):
        pass
