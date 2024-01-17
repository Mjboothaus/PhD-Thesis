import streamlit as st
from .helper_functions import read_render_markdown_file


class StreamlitApp:
    def __init__(self, title, subtitle, about_info, markdown_path):
        self.title = title
        self.subtitle = subtitle
        self.about_info = about_info
        self.markdown_path = markdown_path
        self.configure_page()

    def configure_page(self):
        st.set_page_config(
            page_title=self.title, layout="wide", menu_items={"About": self.about_info}
        )

    def create_app_header(self):
        st.header(self.title)
        if self.subtitle is not None:
            st.subheader(self.subtitle)

    def render_markdown(self):
        read_render_markdown_file(self.markdown_path, output="streamlit")


# Usage
app = StreamlitApp(
    title="Charged fluids near interfaces",
    subtitle="Integral Equation Theory",
    about_info="Created with love & care at DataBooth - www.databooth.com.au",
    markdown_path="docs/app_main.md",
)

app.create_app_header()
app.render_markdown()
