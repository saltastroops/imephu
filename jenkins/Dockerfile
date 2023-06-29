FROM python:3.10 AS requirements

RUN apt-get update -y

RUN curl -sSL https://install.python-poetry.org | POETRY_HOME=/usr/local python3 -

RUN pip install wheel
RUN pip install allure-pytest bandit black mypy pytest ruff
