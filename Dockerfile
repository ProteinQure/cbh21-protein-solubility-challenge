FROM fedora:33

# install python toolchain
RUN sudo dnf install conda python3 python3-devel python3-pip -y

# install build toolchain
RUN sudo dnf install gcc -y

WORKDIR /home/biolib

# setup conda
COPY conda-environment.yaml .
RUN conda env create -f conda-environment.yaml

SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

RUN pip install freesasa atomium temppathlib

# copy all the code
COPY . .
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "myenv", "python3", "predict.py"]
