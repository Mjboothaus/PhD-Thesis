# Docs: https://just.systems/man/en/

project_name := "phd-thesis"

# set dotenv-load

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

rm-python-venv:
	rm -rf .venv/


foo-can-remove-test-only:
	#!/usr/bin/env bash
	set -euxo pipefail
	x=hello
	echo $x
