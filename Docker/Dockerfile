FROM ubuntu:20.04

RUN apt-get update && apt-get install --yes --no-install-recommends \
    wget \
    locales \
    zip \
    unzip \
    bash \
    git 
RUN apt install -y python3-pip


RUN wget -c https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz && \
    tar -xvzf kallisto_linux-v0.44.0.tar.gz && \
    cd kallisto_linux-v0.44.0 && \
    chmod +x kallisto && \
    cp -r kallisto /usr/local/bin
#RUN apt-get install -y kallisto
#=0.44.0
#CMD which kallisto


# Yes, do this twice so it get's cached
RUN pip3 install --upgrade pip
RUN pip3 install RSeQC
RUN pip3 install how_are_we_stranded_here


#RUN pip install numpy==1.11.1
#RUN pip install pandas==0.18.1