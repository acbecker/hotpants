FROM buildpack-deps:focal AS builder
RUN apt update && apt install -y libccfits-dev
RUN git clone https://github.com/acbecker/hotpants.git /workspace
WORKDIR /workspace
RUN git reset --hard 63e126081509a52e96260acdcf2f77b31096c173  # Version 5.1.11
RUN make clean && make

FROM ubuntu:focal
RUN apt update && apt install -y libccfits-dev
COPY --from=builder /workspace/extractkern /workspace/hotpants /workspace/maskim /usr/local/bin/
