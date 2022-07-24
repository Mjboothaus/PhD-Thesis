FROM python:3.9.12
RUN apt-get update -y

# remember to expose the port your app'll be exposed on.
EXPOSE 8080

RUN pip install -U pip

COPY requirements-deploy.txt requirements.txt
RUN pip install -r requirements.txt

# copy into a directory of its own (so it isn't in the toplevel dir)
RUN mkdir -p /app
COPY docs app/docs
COPY src app/src
COPY data app/data
WORKDIR /app

# run it!
ENTRYPOINT ["streamlit", "run", "src/Main.py", "--server.port=8080", "--server.address=0.0.0.0"]