FROM mambaorg/micromamba:2.0.5

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /opt/genomatch

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/genomatch

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.12 \
    bcftools \
    bcftools-liftover-plugin \
    samtools \
    numpy \
    pandas \
    pysam \
    pytest \
    pyyaml \
    pgenlib \
    pip \
    && micromamba clean --all --yes

RUN python -m pip install --no-deps .

RUN chmod +x /opt/genomatch/src/genomatch/*.py

ENV PATH="/opt/genomatch/src/genomatch:${PATH}" \
    PYTHONUNBUFFERED=1

WORKDIR /work

CMD ["bash"]
