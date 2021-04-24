FROM python:3.8-slim
WORKDIR /home/biolib
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY predict.py .
COPY data/test.zip data/
ENTRYPOINT ["python3", "predict.py"]
