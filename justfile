# Docs: https://just.systems/man/en/

project_name := "michael-booth-phd-thesis"
app_py := "src/Main.py"
server_port := "8080"
gcp_region := "asia-southeast2"

set dotenv-load

# show available commands
help:
  @just -l


# Create the local Python venv (.venv_{{project_name}}) and install requirements(.txt)

dev-venv:
	#!/usr/bin/env bash
	pip-compile requirements-dev.in
	python3 -m venv .venv_dev_{{project_name}}
	. .venv_dev_{{project_name}}/bin/activate
	# pip install pip==18.1.0
	python3 -m pip install --upgrade pip
	pip install -r requirements.txt
	python -m ipykernel install --user --name .venv_dev_{{project_name}}
	echo -e '\n' source .venv_dev_{{project_name}}/bin/activate '\n'


# Note: no Jupyter or pytest etc in deploy
deploy-venv:
	#!/usr/bin/env bash
	pip-compile requirements-deploy.in -o requirements-deploy.txt
	python3 -m venv .venv_deploy_{{project_name}}
	. .venv_deploy_{{project_name}}/bin/activate
	python3 -m pip install --upgrade pip
	pip install -r requirements-deploy.txt
	echo -e '\n' source .venv_deploy_{{project_name}}/bin/activate '\n'


update-dev-reqs:
	pip-compile requirements-dev.in
	pip install -r requirements-dev.txt --upgrade


update-deploy-reqs:
	pip-compile requirements-deploy.in
	pip install -r requirements-deploy.txt --upgrade

# See custom dvenv command defined in ~/.zshrc

rm-dev-venv:
	#!/usr/bin/env bash
	dvenv
	rm -rf .venv_dev_{{project_name}}


test:
    pytest


# Run app.py (in Streamlit) locally

run: 
    streamlit run {{app_py}} --server.port={{server_port}} --server.address=localhost

# Build and run app.py in a (local) Docker container

container: 
    docker build . -t {{project_name}}
    docker run -p {{server_port}}:{{server_port}} {{project_name}}


# Still in progress - still not STP
gcr-setup:
    gcloud components update
    # gcloud config set region asia-southeast2
    gcloud projects create {{project_name}}
    gcloud beta billing projects link {{project_name}} --billing-account $BILLING_ACCOUNT_GCP
    gcloud services enable cloudbuild.googleapis.com
    gcloud config set project {{project_name}}


# Deploy container to Google Cloud (Cloud Run)

gcr-deploy: 
    gcloud run deploy --source . {{project_name}} --region {{gcp_region}}
    # --image
    # gcloud run deploy {{project_name}} --image [IMAGE]


gcr-list-deployed-url:
    gcloud run services list --platform managed | awk 'NR==2 {print $4}'


gcr-app-disable:   # deleting project does not delete app
    gcloud app versions list


#TODO: 
# gloud init - other stuff?
# gcloud projects list


# Resources:
# - https://stackoverflow.com/questions/59423245/how-to-get-or-generate-deploy-url-for-google-cloud-run-services