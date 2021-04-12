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

COPY 02-create_env.sh environment.yaml /app/
RUN ./02-create_env.sh

RUN conda init bash
RUN echo "conda activate VitaminD_MR" >> ~/.bashrc

COPY *post_conda_steps.* /app/
RUN conda init

ARG GITHUB_PAT=""
RUN echo "${GITHUB_PAT}" > GITHUB_PAT
RUN ./03-post_conda_steps.sh
RUN rm GITHUB_PAT

RUN conda env export -n VitaminD_MR > /app/linux_docker_env.yml

COPY Dockerfile /app/

ENTRYPOINT ["conda", "run", "-n", "VitaminD_MR", "--no-capture-output"]

CMD ["bash"]