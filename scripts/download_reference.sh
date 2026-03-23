#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REF_DIR="${REPO_ROOT}/ref"

GRCH37_FASTA="${REF_DIR}/ucsc/GRCh37/hg19.p13.plusMT.no_alt_analysis_set.fa"
GRCH38_FASTA="${REF_DIR}/ucsc/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHAIN_37_TO_38="${REF_DIR}/chain/hg19ToHg38.over.chain.gz"
CHAIN_38_TO_37="${REF_DIR}/chain/hg38ToHg19.over.chain.gz"

GRCH37_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz"
GRCH38_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
CHAIN_37_TO_38_URL="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_38_TO_37_URL="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz"

log() {
  printf '%s\n' "$*"
}

require_tool() {
  if ! command -v "$1" >/dev/null 2>&1; then
    log "Error: required tool not found: $1"
    exit 1
  fi
}

download_to_stdout() {
  local url="$1"
  if command -v wget >/dev/null 2>&1; then
    wget -O- "$url"
    return
  fi
  if command -v curl >/dev/null 2>&1; then
    curl -LfsS "$url"
    return
  fi
  log "Error: neither wget nor curl is available"
  exit 1
}

download_gzip_if_missing() {
  local output_path="$1"
  local url="$2"
  if [[ -s "${output_path}" ]]; then
    log "Present: ${output_path}"
    return
  fi
  mkdir -p "$(dirname "${output_path}")"
  log "Downloading: ${output_path}"
  local tmp_path="${output_path}.tmp"
  rm -f "${tmp_path}"
  download_to_stdout "$url" | gzip -d > "${tmp_path}"
  mv "${tmp_path}" "${output_path}"
}

download_file_if_missing() {
  local output_path="$1"
  local url="$2"
  if [[ -s "${output_path}" ]]; then
    log "Present: ${output_path}"
    return
  fi
  mkdir -p "$(dirname "${output_path}")"
  log "Downloading: ${output_path}"
  local tmp_path="${output_path}.tmp"
  rm -f "${tmp_path}"
  download_to_stdout "$url" > "${tmp_path}"
  mv "${tmp_path}" "${output_path}"
}

faidx_if_missing() {
  local fasta_path="$1"
  local fai_path="${fasta_path}.fai"
  if [[ -s "${fai_path}" ]]; then
    log "Present: ${fai_path}"
    return
  fi
  log "Indexing: ${fasta_path}"
  samtools faidx "${fasta_path}"
}

main() {
  require_tool gzip
  require_tool samtools
  if ! command -v wget >/dev/null 2>&1 && ! command -v curl >/dev/null 2>&1; then
    log "Error: download_reference.sh requires wget or curl"
    exit 1
  fi

  mkdir -p \
    "${REF_DIR}/ucsc/GRCh37" \
    "${REF_DIR}/ucsc/GRCh38" \
    "${REF_DIR}/chain"

  download_gzip_if_missing "${GRCH37_FASTA}" "${GRCH37_URL}"
  faidx_if_missing "${GRCH37_FASTA}"

  download_gzip_if_missing "${GRCH38_FASTA}" "${GRCH38_URL}"
  faidx_if_missing "${GRCH38_FASTA}"

  download_file_if_missing "${CHAIN_37_TO_38}" "${CHAIN_37_TO_38_URL}"
  download_file_if_missing "${CHAIN_38_TO_37}" "${CHAIN_38_TO_37_URL}"

  log "Required UCSC reference assets are installed under ${REF_DIR}"
}

main "$@"
