# Docs: https://just.systems/man/en/


# Show available recipes
help:
  @just -l

update-st-config:
    streamlit config show > .streamlit/config.toml

app app_name="src/Main.py":
    streamlit run {{app_name}} --server.address=localhost

reqs:
    pdm export --o requirements.txt --without-hashes --prod


# Run bulk-fluid (pyOz) - default 1 component LJ
bulk-fluid-pyoz input_file="lj/nrcg-lj-1comp.in":
	#!/usr/bin/env bash
	cd src/pyoz
	python pyoz.py -i tests/{{input_file}}


# Build and run app.py in a (local) Docker container
container project_name="phd_thesis" server_port="8080": 
    docker build . -t {{project_name}}
    docker run -p {{server_port}}:{{server_port}} {{project_name}}

test:
    pytest

# TODO: Automate deploy to Render.com and/or Railway.app
