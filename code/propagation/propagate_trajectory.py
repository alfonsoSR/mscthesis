# Goal: Define propagation setup for closed-loop estimation
# Subgoal: Find a propagator that is sufficiently good for MEX flyby: ORMM - ROB
# I want an error wrt ephemerides that is consistent with the difference between
# ROB and ORMM, which is under one km over the time-span of the observations
# I will be propagating between 2013-12-27T00:00:00 and 2014-01-01T00:00:00


from tudatpy.astro import frame_conversion as tframe
from proptools import config as pcon, io as pio
import argparse
from pathlib import Path
from tudatpy.dynamics import (
    propagation as tprop,
    simulator as tsim,
)
from tudatpy.interface import spice
import numpy as np
from proptools import propagation as prop
from tudatpy.dynamics import environment_setup as tenvs

Parser = argparse.ArgumentParser()
Parser.add_argument("config_file", help="Path to configuration file")


def extract_simulation_output(
    simulation: tprop.SingleArcSimulationResults, config: pcon.PropSettings
) -> pio.PropagationOutput:

    # Get ephemerides of spacecraft at propagation epochs
    reference_state_j2000 = np.array(
        [
            spice.get_body_cartesian_state_at_epoch(
                target_body_name=config.env.spacecraft,
                observer_body_name=config.center.name,
                reference_frame_name=config.env.global_frame_orientation,
                aberration_corrections="NONE",
                ephemeris_time=epoch,
            )
            for epoch in simulation.state_history_time_object.keys()
        ]
    )

    # Propagation epochs as float
    propagation_epochs = np.array(list(simulation.state_history.keys()))

    # Propagated state vector (Cartesian global frame)
    propagated_state_j2000 = np.array(list(simulation.state_history.values()))

    # Get states in RSW frame (ephemerides)
    reference_state_rsw = np.zeros_like(reference_state_j2000)
    propagated_state_rsw = np.zeros_like(propagated_state_j2000)
    for idx, rstate in enumerate(reference_state_j2000):

        # Calculate rotation matrix
        rotation_matrix = tframe.inertial_to_rsw_rotation_matrix(rstate)

        # Reference state in RSW
        _refpos_rsw = rotation_matrix @ rstate[:3]
        _refvel_rsw = rotation_matrix @ rstate[3:]
        reference_state_rsw[idx] = np.array(
            [_refpos_rsw, _refvel_rsw]
        ).flatten()

        # Propagated state in RSW
        _ppos_rsw = rotation_matrix @ propagated_state_j2000[idx][:3]
        _pvel_rsw = rotation_matrix @ propagated_state_j2000[idx][3:]
        propagated_state_rsw[idx] = np.array([_ppos_rsw, _pvel_rsw]).flatten()

    # Dependent variables
    dependent_variable_values = np.array(
        [x for x in simulation.dependent_variable_history.values()]
    ).T
    if len(dependent_variable_values.shape) == 1:
        dependent_variable_values = dependent_variable_values[None, :]

    dependent_variables: dict[str, np.ndarray] = {}
    cidx: int = 0
    for name, use in config.variables.variables.items():
        if use:
            dependent_variables[name] = dependent_variable_values[cidx]
            cidx += 1

    # Save output
    return pio.PropagationOutput(
        epochs=propagation_epochs,
        cstate_j2000=propagated_state_j2000,
        rstate_j2000=reference_state_j2000,
        cstate_rsw=propagated_state_rsw,
        rstate_rsw=reference_state_rsw,
        dvars=dependent_variables,
        number_of_function_evaluations=simulation.total_number_of_function_evaluations,
        ustate=np.array(list(simulation.unprocessed_state_history.values())).T,
    )


if __name__ == "__main__":

    # Load configuration
    config_path = Path(Parser.parse_args().config_file).resolve()
    if not config_path.exists():
        raise FileNotFoundError("Configuration file does not exist")
    config = pcon.PropSettings(config_path)

    # Load metakernel (Should be next to configuration)
    metakernel = Path(config_path.parent / "metak.tm").resolve()
    try:
        spice.load_kernel(str(metakernel))

        # Define system of bodies from configuration
        bodies = prop.system_of_bodies_from_config(config)

        # Define propagation settings from configuration
        propagator = prop.propagator_settings_from_config(config, bodies)

        # Propagate equations of motion
        simulator = tsim.create_dynamics_simulator(
            bodies=bodies,
            propagator_settings=propagator,
            simulate_dynamics_on_creation=True,
        )
        assert isinstance(simulator, tsim.SingleArcSimulator)
        simulation = simulator.propagation_results
        assert isinstance(simulation, tprop.SingleArcSimulationResults)

        # Extract and save output
        output = extract_simulation_output(simulation, config)
        output.save_to_file(config_path.parent / "results.pkl")

    finally:
        spice.clear_kernels()
