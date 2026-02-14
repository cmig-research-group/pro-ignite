FROM ubuntu:latest
USER root
SHELL ["/bin/bash", "-c"] 

RUN apt-get update
RUN apt-get install -y sudo
RUN apt-get install unzip
RUN apt-get install -y libxt6

RUN apt-get update && apt-get install dcmtk -y
RUN apt-get update && apt-get install docker.io -y

WORKDIR /pro-ignite
RUN mkdir ./tmp
COPY ./MATLAB_Runtime_R2022a_Update_9_glnxa64.zip ./tmp/
RUN unzip ./tmp/MATLAB_Runtime_R2022a_Update_9_glnxa64.zip -d ./tmp/ 
RUN /pro-ignite/tmp/install -agreeToLicense yes -destinationFolder /pro-ignite -mode silent
RUN rm -r ./tmp

COPY ./compiled/pro_ignite .
COPY ./compiled/run_pro_ignite.sh .

ENTRYPOINT ["/pro-ignite/run_pro_ignite.sh", "/pro-ignite/v912"]
