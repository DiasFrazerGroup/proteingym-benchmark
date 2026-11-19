FROM python:3.12-slim-bookworm

# install git, curl, certificates
RUN apt-get update && apt-get install -y git curl ca-certificates && rm -rf /var/lib/apt/lists/*

# https://docs.astral.sh/uv/guides/integration/docker/#available-images
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /opt/program

# Copy and install dependencies
COPY README.md README.md
COPY pyproject.toml pyproject.toml
COPY uv.lock uv.lock
COPY src/ src/
RUN uv sync --no-cache

# ProteinMPNN vanilla weights v_48_020 (48 neighbours, 0.20 A training noise; upstream default), pinned to the
# upstream commit and by sha256 (byte-identical to files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_020.pt)
ARG PMPNN_URL=https://github.com/dauparas/ProteinMPNN/raw/8907e6671bfbfc92303b5f79c4b5e6ce47cdef57
RUN mkdir -p /opt/program/data && \
    curl -fL ${PMPNN_URL}/vanilla_model_weights/v_48_020.pt -o /opt/program/data/proteinmpnn_v_48_020.pt && \
    echo "c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd  /opt/program/data/proteinmpnn_v_48_020.pt" | sha256sum -c -

ENTRYPOINT ["uv", "run", "proteinmpnn-window"]
