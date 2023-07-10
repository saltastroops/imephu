FROM python:3.10 AS requirements

RUN apt-get update -y

RUN curl -sSL https://install.python-poetry.org | POETRY_HOME=/usr/local python3 -
RUN curl https://pyenv.run | bash

# Adapted from https://gist.github.com/jprjr/7667947
ENV HOME /root
ENV PYENV_ROOT $HOME/.pyenv
ENV PATH $PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH

# Create the folder for the virtual Python environments.
# Note: To avoid re-installation of Python and Python packages every time this
# Dockerfile is used, you should map the /venv folder to a volume.
RUN mkdir /venv

# Install all the supported Python versions.
RUN pyenv install 3.8.0
RUN pyenv install 3.9.0
RUN pyenv install 3.10.0
