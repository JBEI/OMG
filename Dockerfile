FROM jupyter/scipy-notebook:latest

COPY ./Pipfile ./Pipfile.lock /home/jovyan/

RUN pip install pipenv \
 && pipenv install --system --deploy

