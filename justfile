# Docs: https://just.systems/man/en/

project_name := "phd-thesis"
app_py := "src/Main.py"
server_port := "8080"

set dotenv-load

# show available commands
help:
  @just -l

# create the local Python venv (.venv_{{project_name}}) and install requirements(.txt)
setup-python-venv:
	#!/usr/bin/env bash
	pip-compile requirements.in
	python3 -m venv .venv_{{project_name}}
	. .venv_{{project_name}}/bin/activate
	python3 -m pip install --upgrade pip
	pip install -r requirements.txt
	@python -m ipykernel install --user --name .venv_{{project_name}}
	echo -e '\n' source .venv_{{project_name}}/bin/activate '\n'

update-reqs:
    pip-compile requirements.in
    pip install -r requirements.txt --upgrade

rm-python-venv:
	rm -rf .venv/

# run app.py (in Streamlit) locally
run: 
    streamlit run {{app_py}} --server.port={{server_port}} --server.address=localhost
