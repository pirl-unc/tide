FROM ubuntu:20.10
# This image is for running TIDE: https://github.com/liulab-dfci/TIDEpy


ARG recent_commit=898f609af1f5130e8491173bce3a21cfbdd882d1 # most recent as of 20210707

RUN apt-get update && \
    apt-get install -y \
      git \
      nano \
      python3-pip \
      unzip \
      wget && \
    apt-get clean

RUN wget --quiet https://github.com/liulab-dfci/TIDEpy/archive/${recent_commit}.zip && \
    unzip ${recent_commit}.zip && \
    rm ${recent_commit}.zip

RUN mv TIDEpy-${recent_commit} TIDEpy

RUN pip3 install numpy && \
    pip3 install TIDEpy

