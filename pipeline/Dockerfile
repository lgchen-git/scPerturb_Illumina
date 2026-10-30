# Locally provisioned, digest-pinned Illumina DRAGEN/PIPseq analysis base.
FROM illumina-pipseq-scrna:4.5.4-1-reports

USER root
COPY requirements-analysis.txt /tmp/requirements-analysis.txt
RUN /opt/conda/envs/illumina_scrna/bin/python -m pip install --no-cache-dir --no-deps -r /tmp/requirements-analysis.txt \
    && rm -f /tmp/requirements-analysis.txt
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONPATH=/project/pipeline
WORKDIR /project
ENTRYPOINT ["/opt/conda/envs/workflow/bin/snakemake"]
