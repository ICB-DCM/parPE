"""Reading parPE input files"""
import amici
import h5py
import numpy as np


class DataProvider:
    """Interface to HDF5 input files for parPE parameter estimation

    Attributes:
        filename: hdf5 file to access

    """

    def __init__(self, filename: str):
        self.filename = filename

    def apply_solver_settings(self, solver: 'amici.Solver'):
        """Apply settings to Solver"""

        amici.readSolverSettingsFromHDF5(
            self.filename, solver.get(), '/amiciOptions')

    def apply_model_settings(self, model: 'amici.Model'):
        """Apply settings to Model"""

        amici.readModelDataFromHDF5(
            self.filename, model.get(), '/amiciOptions')

    def get_expdata_for_condition(
            self, model: 'amici.Model', condition_idx: int):
        """Get ExpData for the given simulation condition"""

        edata = amici.ExpData(model.get())

        with h5py.File(self.filename, 'r') as f:
            edata.setTimepoints(self._get_timepoints(f, condition_idx))

            y = np.array(edata.getObservedData())
            y = y.reshape((len(edata.getTimepoints()), model.ny), )
            y[:, :] = self._get_measurements(f, condition_idx)
            edata.setObservedData(y.flatten())

            y = np.array(edata.getObservedDataStdDev())
            y = y.reshape((len(edata.getTimepoints()), model.ny), )
            y[:, :] = self._get_measurement_sigmas(f, condition_idx)
            edata.setObservedDataStdDev(y.flatten())

            preeq_cond_idx, sim_cond_idx, reinit_states_flag = \
                self._get_fixed_par_indices(f, condition_idx)
            edata.fixedParameters = self._get_fixed_parameters(f, sim_cond_idx)
            if preeq_cond_idx > 0:
                edata.fixedParametersPreequilibration = \
                    self._get_fixed_parameters(f, preeq_cond_idx)
            if reinit_states_flag > 0:
                edata.reinitializeFixedParameterInitialStates = True

        return edata

    @staticmethod
    def _get_fixed_parameters(f: h5py.File, idx: int):
        return f['/fixedParameters/k'][:, idx]

    @staticmethod
    def _get_fixed_par_indices(f, simulation_idx):
        preeq_cond_idx, sim_cond_idx, reinit_states_flag = \
            f['/fixedParameters/simulationConditions'][simulation_idx, :]
        return preeq_cond_idx, sim_cond_idx, reinit_states_flag

    @staticmethod
    def _get_measurements(f, simulation_idx):
        return f[f'/measurements/y/{simulation_idx}'][:, :]

    @staticmethod
    def _get_measurement_sigmas(f, simulation_idx):
        return f[f'/measurements/ysigma/{simulation_idx}'][:, :]

    @staticmethod
    def _get_timepoints(f, simulation_idx):
        return f[f'/measurements/t/{simulation_idx}']

    def get_starting_points(self, start_idx=None, simulation_idx=None):
        with h5py.File(self.filename, 'r') as f:
            x_full = f['/optimizationOptions/randomStarts'][:, start_idx]
            if simulation_idx is None:
                return x_full

            mapping = self.get_opt_sim_mapping(simulation_idx)
            x = np.full(shape=(mapping.shape[0], 1), fill_value=np.nan)
            for idx, mapped_idx in enumerate(mapping):
                if idx >= 0:
                    x[idx] = x_full[mapped_idx]
                else:
                    raise AssertionError("Implement overrides")
            return x

    @property
    def num_simulations(self):
        with h5py.File(self.filename, 'r') as f:
            return f['/parameters/optimizationSimulationMapping'].shape[1]

    def get_opt_sim_mapping(self, simulation_idx = None):
        with h5py.File(self.filename, 'r') as f:
            return f['/parameters/optimizationSimulationMapping'][:,
                   simulation_idx]
