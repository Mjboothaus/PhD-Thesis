# Docs: https://just.systems/man/en/

project_name := "michael-booth-phd-thesis"
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

setup-deploy-venv:
	#!/usr/bin/env bash
	pip-compile requirements-deploy.in -o requirements-deploy.txt
	python3 -m venv .venv_deploy_{{project_name}}
	. .venv_deploy_{{project_name}}/bin/activate
	python3 -m pip install --upgrade pip
	pip install -r requirements-deploy.txt
	echo -e '\n' source .venv_deploy_{{project_name}}/bin/activate '\n'

update-reqs:
    #pip-compile requirements.in
    #pip install -r requirements.txt --upgrade

rm-python-venv:
	rm -rf .venv/

# run app.py (in Streamlit) locally
run: 
    streamlit run {{app_py}} --server.port={{server_port}} --server.address=localhost

# build and run app.py in a (local) docker container
run-container: 
    docker build . -t {{project_name}}
    docker run -p {{server_port}}:{{server_port}} {{project_name}}

# in progress
gcloud-setup:
    gcloud components update
    # gcloud config set region asia-southeast2
    gcloud projects create {{project_name}}
    gcloud beta billing projects link {{project_name}} --billing-account $BILLING_ACCOUNT_GCP
    gcloud services enable cloudbuild.googleapis.com
    gcloud config set project {{project_name}}


# deploy container (including app.py) to Google Cloud (App Engine)
#gcloud-deploy-app-engine:
#    # gcloud projects delete {{project_name}}
#    # gcloud projects create {{project_name}}
#    gcloud config set project {{project_name}}

# deploy container to Google Cloud (Cloud Run)
gcloud-deploy-cloud-run: 
    gcloud run deploy --source . {{project_name}}
    # --image
    # gcloud run deploy {{project_name}} --image [IMAGE]


# gloud init - other stuff?


gcloud-view:
    gcloud app browse


gcloud-app-disable:   # deleting project does not delete app
    gcloud app versions list

# gcloud app versions stop {{VERSION.ID}}
# gcloud projects list