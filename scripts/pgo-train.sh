#!/bin/sh
# Runs an instrumented whittle over every hot path of the training corpus, so the
# profile covers each of them. Every run is single-threaded, so the counters
# record the work rather than the thread scheduling of the machine that trained.
#
# Usage: scripts/pgo-train.sh <whittle-binary> <corpus-directory>
#
# The corpus comes from `cargo run --example gen-training-data`.

set -eu

if [ "$#" -ne 2 ]; then
    printf >&2 'usage: %s <whittle-binary> <corpus-directory>\n' "$0"
    exit 2
fi

whittle=$1
corpus=$2

if [ ! -x "$whittle" ]; then
    printf >&2 'not an executable binary: %s\n' "$whittle"
    exit 1
fi

for file in reads.fastq reads.fastq.gz reads.bam; do
    if [ ! -r "$corpus/$file" ]; then
        printf >&2 'missing corpus file: %s\n' "$corpus/$file"
        printf >&2 'generate it with: cargo run --example gen-training-data -- %s\n' "$corpus"
        exit 1
    fi
done

out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT HUP INT TERM

# Runs one training pass and stops the script on the first failure, naming the
# flags that failed.
run() {
    printf '\ntrain: %s\n' "$*"
    if ! "$whittle" -t 1 "$@"; then
        printf >&2 'training run failed: %s\n' "$*"
        exit 1
    fi
}

# FASTQ to FASTQ, quality trimming and the length filter.
run -i "$corpus/reads.fastq" -o "$out/trimmed.fastq" --qual-trim 10 -l 200

# FASTQ.gz to FASTQ.gz: gzip decode, fixed crop, and the parallel encoder.
run -i "$corpus/reads.fastq.gz" -o "$out/trimmed.fastq.gz" --head-crop 25 --tail-crop 25 -q 8

# Adapter search over the ONT preset, with chimera splitting.
run -i "$corpus/reads.fastq" -o /dev/null --out-format fastq --adapter-preset ont

# BAM to BAM: tag rewriting through a crop, plus a quality filter.
run -i "$corpus/reads.bam" -o "$out/trimmed.bam" --head-crop 30 --tail-crop 30 -q 9 --update-moves

# BAM to BAM with no trimming stage, the filter-only fast path.
run -i "$corpus/reads.bam" -o "$out/filtered.bam" -l 500 -q 9

# BAM to FASTQ.gz, which carries the aux tags into the headers.
run -i "$corpus/reads.bam" -o "$out/converted.fastq.gz" --qual-trim 9

# JSON summary over a splitting quality stage.
run -i "$corpus/reads.bam" -o /dev/null --out-format bam --qual-split 10 \
    --summary-json "$out/summary.json"

printf '\ntraining complete\n'
