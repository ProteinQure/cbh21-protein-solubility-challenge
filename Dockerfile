FROM python:3.8-slim
WORKDIR /home/biolib
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY predict.py .
ENTRYPOINT ["python3", "predict.py"]

RUN R -e "install.packages(c('caret', 'doMC', 'tidyverse',''), dependecties = TRUE)" 
