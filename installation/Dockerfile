# ref https://github.com/tebeka/pythonwise/blob/master/docker-miniconda/Dockerfile
FROM ubuntu:18.04
WORKDIR /home/
# System packages 
RUN apt-get update && apt-get install -y curl 

ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get -qq update && apt-get -qq install -y curl libcairo2-dev ghostscript texlive-latex-extra texlive-extra-utils vim less 

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

WORKDIR /app

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

COPY create_env.sh environment.yaml /app/
RUN ./create_env.sh

RUN conda init bash
RUN echo "conda activate VitaminD_MR" >> ~/.bashrc

COPY post_conda_steps.* /app/
RUN conda init
RUN ./post_conda_steps.sh

RUN conda env export -n VitaminD_MR > /app/linux_env.yml



COPY Dockerfile /app/
ENTRYPOINT ["conda", "run", "-n", "VitaminD_MR", "--no-capture-output", "/bin/bash", "-c"]
