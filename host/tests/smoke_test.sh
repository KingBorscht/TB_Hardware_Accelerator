#!/usr/bin/env sh
set -eu

build_dir="${1:-build/host}"
"${build_dir}/tb_mutation_caller" \
  host/tests/fixtures/references.tsv \
  host/tests/fixtures/catalogue.tsv \
  host/tests/fixtures/alignments.tsv \
  /tmp/tb_mutation_calls.tsv
diff -u host/tests/expected_calls.tsv /tmp/tb_mutation_calls.tsv
