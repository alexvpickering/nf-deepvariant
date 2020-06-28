# see github: alexvpickering/nf-gatk4 for Dockerfile
FROM gatk4
COPY --from=google/deepvariant:0.10.0-gpu / /