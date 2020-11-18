# see github: alexvpickering/nf-gatk4 for Dockerfile
FROM gatk4
COPY --from=google/deepvariant:0.10.0 / /

RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip

ENV SNPEFF='/gatk/snpEff/snpEff.jar'
ENV SNPGEN='GRCh38.99'

RUN java -jar $SNPEFF download $SNPGEN