# ref https://github.com/tebeka/pythonwise/blob/master/docker-miniconda/Dockerfile
FROM ubuntu:18.04
WORKDIR /home/
# System packages 
RUN apt-get update && apt-get install -y curl 

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get -qq update && apt-get -qq install -y curl libcairo2-dev ghostscript texlive-latex-extra texlive-extra-utils vim less 

# Install miniconda to /miniconda
WORKDIR /app

COPY 01-install_conda.sh /app/
RUN ./01-install_conda.sh

ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# Make RUN commands use `bash --login`:
#SHELL ["/bin/bash", "--login"]


COPY 02.1-prepare_yaml.sh 02.2-create_env.sh environment.yaml /app/
RUN ./02.1-prepare_yaml.sh
RUN ./02.2-create_env.sh

RUN conda init bash
RUN echo "conda activate MRutils" >> ~/.bashrc

COPY *post_conda_steps.* /app/
RUN conda init

# some of the commands are rate limited and may fail in the absence of a GITHUB_PAT. This PAT 
# doesn't need to have ANY authorizations, it simply identifies the user who is downloading from 
# public repositories.

# if you want to use this option, you must build the docker and provide a --build-arg GITHUB_PAT 
# where  GITHUB_PAT is an environment variable containing the PAT. so for example
# GITHUB_PAT=$(cat github_pat_file) docker build --build-arg GITHUB_PAT .
# would work if github_pat_file is a local file containing the pat.
# if building in a CI, you'll need to provide this PAT via a "secret" mechanism.

ARG GITHUB_PAT=""

RUN GITHUB_PAT=${GITHUB_PAT} ./03-post_conda_steps.sh

RUN conda env export -n MRutils > /app/linux_docker_env.yml

COPY Dockerfile /app/

ENTRYPOINT ["conda", "run", "-n", "MRutils", "--no-capture-output"]

CMD ["bash"]