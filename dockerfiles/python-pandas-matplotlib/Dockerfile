FROM ubuntu:20.04
MAINTAINER sminot@fredhutch.org

ADD requirements.txt /share/
RUN apt update && \
    apt install -y build-essential python3 python3-pip && \
    cd /share && \
    pip3 install -r requirements.txt
