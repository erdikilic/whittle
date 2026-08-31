# Contributing

## Building and testing

whittle needs Rust 1.91 or newer; `rust-toolchain.toml` pins the exact stable the
project builds with.

```bash
cargo build --release
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

There's no external `htslib` to install for a normal build. The `rust-htslib`
dependency is dev-only, used by the decode-equivalence tests described in
[docs/tags.md](docs/tags.md) as an independent `MM`/`ML` oracle.

Some tests are `#[ignore]`d because they need a real uBAM. Point them at one:

```bash
WHITTLE_UBAM=/path/to/real.ubam cargo test --test bam_mods_oracle -- --ignored
```

## Generated artifacts

`man/whittle.1` is rendered from the live clap definition, so regenerate and
commit it after changing any flag or its help text:

```bash
cargo run --example gen-man
```

`clap_mangen` is a dev-dependency and the generator is an example, so neither
reaches the shipped binary. The release workflow re-renders the page per build,
which is what stamps the correct version into it.

## Commit messages

Commit subjects follow [Conventional Commits](https://www.conventionalcommits.org/),
for example `feat(adapter): add a trimming mode`, `fix: handle empty input`, or
`perf: improve processing throughput`. Enable the repository's versioned
validation hook once after cloning:

```bash
git config core.hooksPath .githooks
```

The hook rejects invalid subjects before a commit is created. Git won't enable
repository-provided hooks on its own, for security reasons, so CI validates every
pushed or pull-request commit as the repo-wide backstop. You can bypass the local
hook with `--no-verify`; to make it non-optional, configure branch protection to
require the `commit-message` CI job before merging.

Give manual merge commits a compliant subject instead of Git's default `Merge ...`
message:

```bash
git merge --no-ff feature-branch -m "perf: integrate throughput improvements"
```

## Style

- American English throughout: comments, doc comments, identifiers, and user-facing
  strings.
- No em dashes or en dashes in comments, prose, or commit messages. Use commas,
  periods, semicolons, colons, or parentheses.
- Comments explain why, not what. A comment that restates the line below it should
  be deleted.
- Third person, not first: "the caller drops the tag", not "we drop the tag".

## Log messages

- Start the message with a capital letter, at every level. The message is prose;
  the structured fields after it are data.
- Never restate the level in the text. The formatter already prints `[WARN]`, so
  a message beginning "warning:" says it twice.
- Put values in fields rather than in the sentence, so they can be filtered and
  parsed: `tracing::debug!(reads = n, "Processing finished")`, not
  `tracing::debug!("Processing finished, {n} reads")`.
- Field values are data and stay lowercase: `action="trim 5'"`, `reason="too short"`.
- Reuse one wording for one concept. A rejection reason comes from
  `DropReason::label` so the summary line and the trace event cannot disagree.
