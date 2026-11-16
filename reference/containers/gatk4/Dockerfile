# GATK 4.6.2.0 for linux/arm64. Same tool versions as sarek 3.10.0's
# gatk4_gcnvkernel image minus gcnvkernel (MKL-only, x86-only; only used for gCNV).
FROM mambaorg/micromamba:1.5.10-noble
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      bioconda::gatk4=4.6.2.0 bioconda::htslib=1.21 bioconda::samtools=1.21 conda-forge::procps-ng \
 && micromamba env export -n base --explicit > /opt/conda/conda-explicit.txt \
 && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
USER root
