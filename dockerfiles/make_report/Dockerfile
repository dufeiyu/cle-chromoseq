FROM python:3
MAINTAINER Feiyu Du <fdu@wustl.edu>

LABEL python3 image for make_report task

RUN pip install cyvcf2
RUN mkdir /opt/files/

COPY ../../scripts/make_report.py /usr/local/bin/make_report.py
COPY ../../accessory_files/QCReferenceRanges.json /opt/files/QCReferenceRanges.json

RUN chmod 555 /usr/local/bin/make_report.py
RUN chmod 444 /opt/files/QCReferenceRanges.json
