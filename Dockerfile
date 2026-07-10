FROM containers.mathworks.com/matlab-runtime:r2022a

USER root
SHELL ["/bin/bash", "-c"] 

RUN apt-get update
RUN apt-get install -y sudo
RUN apt-get install -y libxt6

RUN apt-get update && apt-get install dcmtk -y
RUN apt-get update && apt-get install docker.io -y

WORKDIR /pro-ignite

COPY ./compiled/pro_ignite .
COPY ./run_pro_ignite.sh .

ENV AGREE_TO_MATLAB_RUNTIME_LICENSE="yes"

ENTRYPOINT ["/pro-ignite/run_pro_ignite.sh", "/opt/matlabruntime/v912"]
