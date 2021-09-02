FROM continuumio/miniconda
# docker build -t mixture .
# docker run --runtime=nvidia -p 7777:8888 -v ${PWD}:/DP -it jesselee/pytorchls bash
# docker run -v ${dir}:/mixture -it mixture bash
# docker run -v D:\tmp\windows:/mixture -it mixture bash


RUN apt update
RUN apt install -y build-essential
RUN apt install gcc
# RUN apt-get install -y manpages-dev

RUN conda install python=3.7.6
RUN conda install -c r r
RUN conda install scipy

RUN pip install jenkspy
RUN conda install rpy2
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda install pysam

RUN conda install samtools
RUN conda install bowtie2
RUN mkdir -p /usr/share/man/man1
RUN apt install -y default-jre

COPY . /mixture
