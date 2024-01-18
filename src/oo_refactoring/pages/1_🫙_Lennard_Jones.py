from oo_refactoring.base_streamlit_app import BaseStreamlitApp
from oo_refactoring.parameters import FluidInitialiser

class LennardJonesApp(BaseStreamlitApp):
    def __init__(self, fluid_symbol):
        fluid_initialiser = FluidInitialiser(config_path="path/to/fluids.toml")
        super().__init__(fluid_initialiser.get_fluid(fluid_symbol))

# Usage
app = LennardJonesApp(fluid_symbol="lj")
app.run()