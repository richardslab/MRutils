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

COPY install_conda.sh /app/
RUN ./install_conda.sh

ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda


# Make RUN commands use `bash --login`:
#SHELL ["/bin/bash", "--login"]

COPY create_env.sh environment.yaml /app/
RUN ./create_env.sh

RUN conda init bash
RUN echo "conda activate VitaminD_MR" >> ~/.bashrc

COPY post_conda_steps.* /app/
RUN conda init
RUN ./post_conda_steps.sh

RUN conda env export -n VitaminD_MR > /app/linux_env.yml


COPY Dockerfile /app/

ENTRYPOINT ["conda", "run", "-n", "VitaminD_MR", "--no-capture-output", "bash", "-c"]

CMD ["bash"]