# See https://stackoverflow.com/questions/68673221/warning-running-pip-as-the-root-user
# for enhancements to Dockerfile e.g. not running as root & in venv

#FROM python:3.9.12
FROM python:3.10-slim-bullseye
RUN apt-get update && apt-get install -y

# remember to expose the port your app'll be exposed on.
EXPOSE 8080

RUN pip install -U pip

COPY requirements-deploy.txt requirements.txt
RUN pip install -r requirements.txt

# copy into a directory of its own (so it isn't in the toplevel dir)
# RUN mkdir -p /app
COPY docs app/docs
COPY src app/src
COPY data app/data
#COPY output app/output
WORKDIR /app

# run it!
ENTRYPOINT ["streamlit", "run", "src/Main.py", "--server.port=8080", "--server.address=0.0.0.0"]
# ENTRYPOINT ["streamlit", "run", "src/Main.py", "--server.port=8080", "--server.address=0.0.0.0", "--server.enableCORS false", "--server.enableXsrfProtection false"]

# See https://discuss.streamlit.io/t/deploying-streamlit-on-gcp-cloud-run-problem-when-using-new-multipage-app-feature/26316/3