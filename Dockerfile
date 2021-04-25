FROM python:3.8
WORKDIR /home/biolib
COPY model.pkl /home/biolib
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY predict.py .
ENTRYPOINT ["python3", "predict.py"]