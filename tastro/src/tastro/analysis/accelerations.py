from .. import (
    config as ncon,
    io as nio,
    environment as nenv,
    propagation as nprop,
)
from ..propagation.interface import PropagationSettings
from dataclasses import dataclass
from pathlib import Path
from tudatpy.interface import spice
from tudatpy.dynamics.simulator import (
    create_dynamics_simulator,
    SingleArcSimulator,
)
from tudatpy.dynamics.propagation import AccelerationModel
from tudatpy.dynamics.propagation_setup import (
    propagator as tprops,
    create_acceleration_models,
)
import numpy as np
from nastro import types as nt, graphics as ng
from tudatpy.astro import time_representation as ttime


@dataclass
class AccelerationsInput:

    source: Path


def accelerations(user_input: AccelerationsInput) -> None:

    # Get configuration from user input
    config = ncon.CaseSetup.from_config_file(user_input.source / "configuration.yaml")

    # Define path to metakernel
    metakernel = Path(user_input.source / "metak.tm").absolute()
    try:
        # Load metakernel
        spice.load_kernel(str(metakernel))

        # Create system of bodies
        bodies = nenv.system_of_bodies_from_config(config)

        propagator = nprop.translational_propagator_settings_from_config(config, bodies)

        acceleration_models: dict[str, dict[str, list[AccelerationModel]]] = (
            propagator.accelerations_map
        )
        labels: dict[str, dict[str, list[str]]] = {"MEX": {}}
        for label in labels:

            for ext, accelerations in acceleration_models[label].items():

                labels["MEX"][ext] = []

                fo = config.propagation.accelerations["MEX"].external[ext]
                if fo.gravitational.present:
                    labels["MEX"][ext].append("gravity")
                if fo.aerodynamic.present:
                    labels["MEX"][ext].append("drag")
                if fo.radiation.present:
                    labels["MEX"][ext].append("radiation")
                if fo.relativistic.present:
                    labels["MEX"][ext].append("relativity")

        sol = nio.PropagationOutput.from_config_file(
            user_input.source / "configuration.yaml"
        )

        gen = PropagationSettings("", config)

        acceleration_settings = gen.acceleration_settings()["MEX"]

        results: dict[str, dict[str, nt.CartesianStateDerivative]] = {}

        for external, accelerations in acceleration_settings.items():

            results[external] = {}

            for idx, acceleration in enumerate(accelerations):

                current_settings = {"MEX": {external: [acceleration]}}

                propagator.reset_and_recreate_acceleration_models(
                    current_settings, bodies
                )
                sim: SingleArcSimulator = create_dynamics_simulator(
                    bodies,
                    propagator,
                    simulate_dynamics_on_creation=True,
                )

                results[external][labels["MEX"][external][idx]] = (
                    nt.CartesianStateDerivative(
                        *np.array(
                            [
                                sim.state_derivative_function(ttime.Time(ti), si).T[0]
                                for ti, si in zip(sol.epochs, sol.cstate_j2000)
                            ]
                        ).T
                    )
                )

        legend_setup = ng.PlotSetup(show_axes=False, legend_location="center")

        with ng.Mosaic("aab") as canvas:

            fig = canvas.subplot(ng.PlotSetup(yscale="log"))
            leg = canvas.subplot(legend_setup)

            for planet, data in results.items():

                for label, item in data.items():

                    fig.line(sol.epochs, item.a_mag)
                    leg.line(0, 0, label=f"{planet} {label}")

            fig.__exit__(0, 0, 0)
            leg.__exit__(0, 0, 0)

        # # Define propagator settings from configuration
        # propagator = nprop.translational_propagator_settings_from_config(
        #     config, bodies
        # )

        # sim = create_dynamics_simulator(
        #     bodies, propagator, simulate_dynamics_on_creation=False
        # )
        # assert isinstance(sim, SingleArcSimulator)

        # print(sim.state_derivative_function)

        # # Get acceleration settings from propagator
        # accelerations = tprops.create_acceleration_models(
        #     body_system=bodies,
        #     selected_acceleration_per_body=propagator.acceleration_settings,
        #     bodies_to_propagate=[config.environment.general.spacecraft],
        #     central_bodies=[config.environment.general.center],
        # )

        # for body, body_accelerations in accelerations.items():

        #     print(body_accelerations)

    finally:
        # Clear kernel pool
        spice.clear_kernels()

    return None
