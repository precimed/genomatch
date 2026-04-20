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
    scipy \
    pysam \
    pytest \
    pyyaml \
    pgenlib \
    pip \
    && micromamba clean --all --yes

RUN python -m pip install --no-deps .

ENV PYTHONUNBUFFERED=1
ENV PATH=/opt/conda/bin:${PATH}

WORKDIR /work

CMD ["bash"]
