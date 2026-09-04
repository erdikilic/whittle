#!/bin/sh
# Builds a profile-guided release binary for one target: instrument, generate the
# training corpus, train, merge the profile, then rebuild against it. The target
# must run natively on this machine, because training executes the instrumented
# binary; nothing here emulates.
#
# Usage: scripts/pgo-build.sh <target-triple>
#
# Needs the llvm-tools rustup component for llvm-profdata. The binary lands where
# a plain `cargo build --release --target <triple>` would put it. Both the
# release workflow and the CI pgo job call this, so CI covers the release path.

set -eu

if [ "$#" -ne 1 ]; then
    printf >&2 'usage: %s <target-triple>\n' "$0"
    exit 2
fi

target=$1
root=$(CDPATH='' cd -- "$(dirname -- "$0")/.." && pwd)
host=$(rustc -vV | sed -n 's/^host: //p')

target_dir=${CARGO_TARGET_DIR:-$root/target}
mkdir -p "$target_dir"
# rustc resolves -Cprofile-generate and -Cprofile-use against its own working
# directory, which is not the caller's, so both paths are absolute.
target_dir=$(CDPATH='' cd -- "$target_dir" && pwd)
profile_dir=$target_dir/pgo
corpus_dir=$target_dir/training-data
merged=$profile_dir/whittle.profdata

llvm_profdata=$(rustc --print sysroot)/lib/rustlib/$host/bin/llvm-profdata
if [ ! -x "$llvm_profdata" ]; then
    printf >&2 'llvm-profdata not found at %s\n' "$llvm_profdata"
    printf >&2 'install it with: rustup component add llvm-tools\n'
    exit 1
fi

# RUSTFLAGS replaces the rustflags in .cargo/config.toml rather than adding to
# them, so the target-cpu the SIMD gate requires is read back from the config and
# carried into both builds. The config sets it for x86_64 only. Flags the
# environment already carries are kept ahead of it.
inherited=${RUSTFLAGS:-}
cpu_flag=
case "$target" in
    x86_64-*)
        cpu=$(grep -o 'target-cpu=[A-Za-z0-9._-]*' "$root/.cargo/config.toml" | head -n 1)
        if [ -z "$cpu" ]; then
            printf >&2 'no target-cpu in .cargo/config.toml; refusing to build without it\n'
            exit 1
        fi
        cpu_flag="-C $cpu"
        ;;
esac

printf '\n== instrumented build (%s) ==\n' "$target"
rm -rf "$profile_dir"
mkdir -p "$profile_dir"
RUSTFLAGS="$inherited $cpu_flag -Cprofile-generate=$profile_dir" \
    cargo build --release --locked --target "$target"

printf '\n== training corpus ==\n'
# The generator runs on this machine, so it is built for the host: the
# cross-toolchain action exports CARGO_BUILD_TARGET, which would otherwise
# cross-build it. RUSTFLAGS is cleared so .cargo/config.toml supplies target-cpu
# and the generator itself is not instrumented.
env -u CARGO_BUILD_TARGET -u RUSTFLAGS \
    cargo run --locked --example gen-training-data -- "$corpus_dir"

printf '\n== training ==\n'
"$root/scripts/pgo-train.sh" "$target_dir/$target/release/whittle" "$corpus_dir"

printf '\n== merging profile ==\n'
# `set -e` does not see a failing glob expansion, so the raw profiles are counted
# before the merge and an empty profile is rejected after it.
raw_count=$(find "$profile_dir" -name '*.profraw' | wc -l | tr -d '[:space:]')
if [ "$raw_count" -eq 0 ]; then
    printf >&2 'training produced no .profraw files in %s\n' "$profile_dir"
    exit 1
fi
printf 'raw profiles: %s\n' "$raw_count"
find "$profile_dir" -name '*.profraw' -print0 \
    | xargs -0 "$llvm_profdata" merge -o "$merged"
if [ ! -s "$merged" ]; then
    printf >&2 'llvm-profdata produced no profile at %s\n' "$merged"
    exit 1
fi

printf '\n== optimized build (%s) ==\n' "$target"
RUSTFLAGS="$inherited $cpu_flag -Cprofile-use=$merged" \
    cargo build --release --locked --target "$target"

binary=$target_dir/$target/release/whittle
if [ ! -x "$binary" ]; then
    printf >&2 'no binary at %s after the optimized build\n' "$binary"
    exit 1
fi
printf '\n%s\n' "$binary"
