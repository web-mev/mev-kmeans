FROM debian:stretch

RUN apt-get update && \
    apt-get install -y python3-dev python3-pip

# Install some Python3 libraries:
RUN mkdir /opt/software
ADD requirements.txt /opt/software/
ADD run_kmeans.py /opt/software/
RUN chmod +x /opt/software/run_kmeans.py
RUN pip3 install -r /opt/software/requirements.txt

ENTRYPOINT ["/opt/software/run_kmeans.py"]
