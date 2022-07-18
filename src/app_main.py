from pathlib import Path
import streamlit as st


from scipy.constants import epsilon_0, elementary_charge
from scipy.integrate import trapz
import scipy.optimize as so

# Other constants - for model of water

q_h = 0.32983 * elementary_charge
q_o = -2 * q_h


# Other constants - Molten KCl (page 65)

# Temperature
# Concentration
# Bulk direct correlation functions
# Short range wall potential
# Surface potential or charge

# Solving for t(z) = ln g + b \phi(z)
# Newton-GMRES algorithm

# so.newton_krylov

so.anderson()

def read_render_markdown_file(markdown_file):
    try:
        md_text = Path(markdown_file).read_text()
        st.markdown(md_text, unsafe_allow_html=False)
        return markdown_file
    except Exception:
        st.error(f"Error with markdown file: {markdown_file}")
        return None


def create_app_header(app_title, subtitle=None):
    st.set_page_config(
        page_title=app_title,
        layout="wide",
        menu_items={
            "About":
            "Created with love & care at DataBooth - www.databooth.com.au"
        },
    )
    st.header(app_title)
    if subtitle is not None:
        st.subheader(subtitle)
    return None


def create_sidebar():
    st.sidebar.text("Inputs:")
    st.sidebar.number_input("Temperature (C):", value=25, min_value=-298, max_value=2000)
    st.sidebar.number_input(label=r"$\psi_0$ (mV):", value=0, min_value=-1000, max_value=1000)
    nPoint = st.sidebar.number_input("nPoint:", value=4096)
    grid_size = st.sidebar.number_input("Grid size (Angstrom):", value=0.006)
    z_cutoff = round((nPoint - 1) * grid_size, 4)
    st.sidebar.text(f"Cutoff (A): {z_cutoff}")

    return None

create_app_header("Charged fluids near interfaces: Integral equation theory", "Michael J. Booth")


create_sidebar()


read_render_markdown_file("docs/intro.md")
read_render_markdown_file("docs/equation_3_1.md")
read_render_markdown_file("docs/equation_3_2.md")