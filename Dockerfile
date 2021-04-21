FROM python:3.8-slim
WORKDIR /home/biolib
RUN pip install numpy biopython
COPY predict.py .
ENTRYPOINT ["python3", "predict.py"]
