FROM python:3.8.1-alpine3.11

COPY . /app

WORKDIR /app

RUN pip install pipenv
RUN pipenv install dev

EXPOSE 5000

ENTRYPOINT [ "python" ]

CMD [ "app.py" ]
