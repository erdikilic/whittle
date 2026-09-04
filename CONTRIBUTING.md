# Contributing

## Building and testing

whittle needs Rust 1.91 or newer; `rust-toolchain.toml` pins the exact stable the
project builds with.

```bash
cargo build --release
cargo test
cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --check
```

There is no external `htslib` to install for a normal build. The `rust-htslib`
dependency is dev-only, used by the decode-equivalence tests described in
[docs/tags.md](docs/tags.md) as an independent `MM`/`ML` oracle.

Some tests are `#[ignore]`d because they need a real uBAM. `WHITTLE_UBAM` names
one:

```bash
WHITTLE_UBAM=/path/to/real.ubam cargo test --test bam_mods_oracle -- --ignored
```

## Generated artifacts

`man/whittle.1` is rendered from the live clap definition, so it is regenerated
and committed after any change to a flag or its help text:

```bash
cargo run --example gen-man
```

`clap_mangen` is a dev-dependency and the generator is an example, so neither
reaches the shipped binary. The release workflow ships the committed page; the
`Man page is current` CI job fails when it is stale.

## Commit messages

Commit subjects follow [Conventional Commits](https://www.conventionalcommits.org/),
for example `feat(adapter): add a trimming mode`, `fix: handle empty input`, or
`perf: improve processing throughput`. The repository's versioned validation
hook is enabled once after cloning:

```bash
git config core.hooksPath .githooks
```

The hook rejects invalid subjects before a commit is created. Git does not enable
repository-provided hooks on its own, for security reasons, so CI validates every
pushed or pull-request commit as the repo-wide backstop. `--no-verify` bypasses
the local hook; branch protection requiring the `commit-message` CI job makes
validation non-optional.

Manual merge commits take a compliant subject instead of Git's default `Merge ...`
message:

```bash
git merge --no-ff feature-branch -m "perf: integrate throughput improvements"
```

## Style

- American English throughout: comments, doc comments, identifiers, and user-facing
  strings.
- No em dashes, en dashes, or ` -- ` in comments, prose, or commit messages.
  Commas, periods, semicolons, colons, or parentheses take their place. Unicode
  arrows, ellipses, and emoji are spelled out: "BAM-to-FASTQ", "...", "yes".
- Comments state what a construct does and, when it is not evident, the technical
  reason it is built that way. Declarative, present tense, short. A comment that
  restates the line below it is deleted.
- No narrative or scenario framing, no editorializing ("simply", "deliberately",
  "actually"), no change history ("previously", "now", "no longer"; that belongs
  in CHANGELOG.md), and no `# --- section ---` banners: a plain `# Section` label
  instead. TOML comments sit on their own line above the key.
- Third person, not first: "the caller drops the tag", not "we drop the tag".
- `anyhow!`/`bail!` messages start lowercase and name the flag or path involved:
  `--summary-json x is a directory`, `opening input reads.fastq`. `main` prints
  them after `Failed after 12ms: ` (or bare, for a parse error), so a capital
  would land mid-sentence. Flags spelled with a leading dash keep their spelling.

## Log messages

- The message starts with a capital letter, at every level. The message is prose;
  the structured fields after it are data. A message opening on a literal that is
  spelled lowercase keeps that spelling, since it is a token and not a word:
  `--adapter-sample is ignored with --adapter-fasta`, `inferred_1 support=0.82`.
- The level is never restated in the text. The formatter already prints `[WARN]`,
  so a message beginning "warning:" says it twice.
- Values go in fields rather than in the sentence, so they can be filtered and
  parsed: `tracing::debug!(reads = n, "Processing finished")`, not
  `tracing::debug!("Processing finished, {n} reads")`.
- Field values are data and stay lowercase: `action="trim 5'"`, `reason="too short"`.
- One wording serves one concept. A rejection reason comes from
  `DropReason::label` so the summary line and the trace event cannot disagree.
