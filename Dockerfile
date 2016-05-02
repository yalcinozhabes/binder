FROM andrewosh/binder-base

MAINTAINER Yalcin Ozhabes <yalcinozhabes@gmail.com>

USER root

# Add dependency
RUN apt-get update
RUN apt-get install -y graphviz git
ADD jdftx /home/main/jdftx
ADD .ipython /home/main/.ipython
RUN chmod -R u+w /home/main/.ipython


USER main
ENV PYTHONPATH /home/main/jdftx:$PYTHONPATH

# Install requirements for Python 2
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt

# Install requirements for Python 3
RUN /home/main/anaconda2/envs/python3/bin/pip install -r requirements.txt

#install chemview
RUN /home/main/anaconda2/envs/python3/bin/pip install -r requirements.txt
#RUN git clone https://github.com/gabrielelanaro/chemview
#WORKDIR chemview
# RUN /home/main/anaconda2/envs/python3/bin/pip install .
# ADD chemview /home/main/anaconda2/envs/python3/lib/python3.5/site-packages/chemview/
