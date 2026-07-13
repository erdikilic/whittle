#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
    printf >&2 'usage: %s <commit-message-file>\n' "$0"
    exit 2
fi

message_file=$1
if [ ! -r "$message_file" ]; then
    printf >&2 'commit message file is not readable: %s\n' "$message_file"
    exit 2
fi

subject=$(sed -n '1p' "$message_file")
pattern='^(build|chore|ci|docs|feat|fix|perf|refactor|revert|style|test)(\([a-z0-9][a-z0-9._/-]*\))?!?: [^[:space:]].*$'

if printf '%s\n' "$subject" | LC_ALL=C grep -Eq "$pattern"; then
    exit 0
fi

printf >&2 '\nInvalid commit subject:\n  %s\n\n' "$subject"
printf >&2 'Use Conventional Commit format:\n'
printf >&2 '  <type>[optional scope][optional !]: <description>\n\n'
printf >&2 'Allowed types: build, chore, ci, docs, feat, fix, perf, refactor, revert, style, test\n\n'
printf >&2 'Examples:\n'
printf >&2 '  feat(adapter): add conservative inference\n'
printf >&2 '  fix: handle empty FASTQ input\n'
printf >&2 '  perf!: change the processing pipeline\n\n'
exit 1
