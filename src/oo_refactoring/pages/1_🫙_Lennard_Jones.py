from base_streamlit_app import BaseStreamlitApp
from parameters import FluidInitialiser
from pathlib import Path


class LennardJonesApp(BaseStreamlitApp):
    def __init__(self, fluid_symbol):
        fluid_initialiser = FluidInitialiser(
            config_path=Path("src/oo_refactoring/fluid_parameters.toml").as_posix()
        )
        super().__init__(fluid_initialiser.get_fluid(fluid_symbol))


app = LennardJonesApp(fluid_symbol="lj")
app.run()
