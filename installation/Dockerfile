# ref https://github.com/tebeka/pythonwise/blob/master/docker-miniconda/Dockerfile
FROM ubuntu:18.04
WORKDIR /home/
# System packages 
RUN apt-get update && apt-get install -y curl

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

WORKDIR /app

COPY create_env.sh environment.yaml /app/
RUN ./create_env.sh

COPY post_conda_steps.* /app/
RUN conda init
RUN ./post_conda_steps.sh
