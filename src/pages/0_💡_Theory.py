
import streamlit as st
from helper_functions import read_render_markdown_file
from pathlib import Path

# read_render_markdown_file("docs/intro.md", "streamlit")

read_render_markdown_file("docs/app_theory.md", "streamlit")

# for file in sorted([file.as_posix() for file in Path("docs").glob("eq*.md")]):
#    read_render_markdown_file(file, output="streamlit")

#read_render_markdown_file("docs/equation_3_1.md", "streamlit")
#read_render_markdown_file("docs/equation_3_2.md", "streamlit")
#read_render_markdown_file("docs/equation_4_4.md", "streamlit")


#TODO: Try using SessionState & on_click callback for Next/Prev page in the markdown