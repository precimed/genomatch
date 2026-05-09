#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  normalize_bed_contigs_ncbi.sh BED_DIR

Normalize every *.bed file in BED_DIR to canonical NCBI-style contigs:
  1-22, X, Y, MT

Rows on canonical UCSC-style contigs are renamed:
  chr1..chr22 -> 1..22
  chrX        -> X
  chrY        -> Y
  chrM/chrMT  -> MT

Rows already using 1-22, X, Y, or MT are kept unchanged.
Rows on non-canonical contigs such as *_alt, *_fix, *_random, chrUn, etc.
are removed and recorded in:
  normalize_bed_contigs_ncbi.noncanonical.tsv

Per-file counts are recorded in:
  normalize_bed_contigs_ncbi.report.tsv

The script is idempotent: re-running it on already-normalized files produces
the same BED content and refreshes the same report files.
USAGE
}

log() {
  printf '%s\n' "$*"
}

die() {
  log "Error: $*" >&2
  exit 1
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

[[ "$#" -eq 1 ]] || {
  usage >&2
  exit 2
}

BED_DIR="$1"

[[ -d "${BED_DIR}" ]] || die "BED directory does not exist: ${BED_DIR}"
command -v awk >/dev/null 2>&1 || die "required tool not found: awk"
command -v find >/dev/null 2>&1 || die "required tool not found: find"
command -v cmp >/dev/null 2>&1 || die "required tool not found: cmp"
command -v sort >/dev/null 2>&1 || die "required tool not found: sort"

REPORT="${BED_DIR}/normalize_bed_contigs_ncbi.report.tsv"
NONCANONICAL="${BED_DIR}/normalize_bed_contigs_ncbi.noncanonical.tsv"

printf 'file\tinput_rows\tkept_rows\tremoved_noncanonical_rows\tchanged_contig_rows\n' > "${REPORT}"
printf 'file\tline_number\tcontig\treason\n' > "${NONCANONICAL}"

found=0
updated=0
unchanged=0

while IFS= read -r -d '' bed_path; do
  found=1
  tmp_path="$(mktemp "${bed_path}.tmp.XXXXXX")"

  awk -v report="${REPORT}" -v noncanonical="${NONCANONICAL}" '
    BEGIN {
      OFS = "\t"
      input_rows = 0
      kept_rows = 0
      removed_rows = 0
      changed_rows = 0
    }

    /^[[:space:]]*$/ || /^[[:space:]]*#/ || /^[[:space:]]*track([[:space:]]|$)/ || /^[[:space:]]*browser([[:space:]]|$)/ {
      print
      next
    }

    {
      input_rows++
      original_contig = $1
      normalized_contig = original_contig

      if (original_contig ~ /^chr([1-9]|1[0-9]|2[0-2])$/) {
        normalized_contig = substr(original_contig, 4)
      } else if (original_contig == "chrX") {
        normalized_contig = "X"
      } else if (original_contig == "chrY") {
        normalized_contig = "Y"
      } else if (original_contig == "chrM" || original_contig == "chrMT" || original_contig == "M") {
        normalized_contig = "MT"
      }

      if (normalized_contig ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/) {
        if (normalized_contig != original_contig) {
          changed_rows++
          $1 = normalized_contig
        }
        kept_rows++
        print
      } else {
        removed_rows++
        print FILENAME, FNR, original_contig, "noncanonical_contig" >> noncanonical
      }
    }

    END {
      print FILENAME, input_rows, kept_rows, removed_rows, changed_rows >> report
    }
  ' "${bed_path}" > "${tmp_path}"

  if cmp -s "${bed_path}" "${tmp_path}"; then
    rm -f "${tmp_path}"
    unchanged=$((unchanged + 1))
    log "Unchanged: ${bed_path}"
  else
    mv "${tmp_path}" "${bed_path}"
    updated=$((updated + 1))
    log "Updated:   ${bed_path}"
  fi
done < <(find "${BED_DIR}" -maxdepth 1 -type f -name '*.bed' -print0 | sort -z)

[[ "${found}" -eq 1 ]] || die "no .bed files found in: ${BED_DIR}"

log "Done. Updated ${updated} file(s); ${unchanged} already normalized."
log "Report: ${REPORT}"
log "Noncanonical rows: ${NONCANONICAL}"
