from ...core import SetupBase


class ParametersSetup(SetupBase):

    initial_states: bool
    radiation_pressure_coefficient: bool
    drag_coefficient: bool
    arcwise_drag_coefficient: bool = False
    k2_radiation_coefficient: bool = False
    gm_phobos: bool = False
    open_loop_biases: bool = False
