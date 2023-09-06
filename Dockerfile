# syntax=docker/dockerfile:1
# tu może rocker zamiast ubuntu
FROM ubuntu:18.04

RUN apt-get update -y && apt-get -y upgrade

RUN apt-get install -y build-essential wget gfortran xorg-dev libpcre3-dev \
		libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev git fort77 libreadline-dev bc libjpeg-dev curl default-jre nano

RUN apt-get install -y python3.6
RUN ln -s /usr/bin/python3.6 /usr/bin/python

WORKDIR /opt

RUN wget https://cran.r-project.org/src/base/R-3/R-3.4.3.tar.gz && \
    tar xzf R-3.4.3.tar.gz && \
    rm R-3.4.3.tar.gz && \
    cd R-3.4.3/ && \
    ./configure && \
    make && \
    make install

RUN wget https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-50.tar.gz && \
    R CMD INSTALL MASS_7.3-50.tar.gz && \
    wget https://cran.r-project.org/src/contrib/Archive/VGAM/VGAM_1.0-5.tar.gz && \
    R CMD INSTALL VGAM_1.0-5.tar.gz && \
    wget https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.11.2.tar.gz && \
    R CMD INSTALL data.table_1.11.2.tar.gz


RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	tar xjf samtools-1.10.tar.bz2 && \
	rm samtools-1.10.tar.bz2 && \
	cd samtools-1.10 && \
	./configure --prefix $(pwd) && \
	make

ENV PATH=${PATH}:/opt/samtools-1.10

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && \
    tar -zxvf bedtools-2.27.1.tar.gz && \
    rm bedtools-2.27.1.tar.gz && \
    cd bedtools2 && \
    make

ENV PATH=${PATH}:/opt/bedtools2/bin

RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make

ENV PATH=${PATH}:/opt/bwa

RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
  tar -vxjf bcftools-1.10.2.tar.bz2 && \
  cd bcftools-1.10.2 && \
  ./configure --prefix=/where/to/install && \
  make && \
  make install

ENV PATH=${PATH}/opt/bcftools-1.10.2/bin:$PATH
WORKDIR /opt
RUN apt-get install -y python3-pip

RUN apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN git clone https://github.com/HuMingLab/MAPS.git

WORKDIR /app

# todo fix package installation via pip3 and requirements
RUN pip3 install --upgrade pip
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
RUN pip3 install -U setuptools
RUN pip3 install cykhash macs3
RUN pip3 install deeptools
RUN pip3 install pairtools
RUN pip3 install plumbum

RUN curl -s https://get.nextflow.io | bash && mv nextflow /opt/

RUN /opt/nextflow
