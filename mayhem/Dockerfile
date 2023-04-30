FROM debian:bookworm as build 
RUN apt update && apt install build-essential clang -y

COPY . /peanut-gb
WORKDIR /peanut-gb/fuzz
RUN make