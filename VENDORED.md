# Vendored modules

`jennystats/typed.janet` is vendored from [janet-kit](https://github.com/jennystats/janet-kit)
(AGPLv3) as a regular file - janet-kit 0.1.2 or newer. Vendoring keeps
this repository self-contained: the suite runs with no install step
(the module resolves in-repo from the repo root's `jennystats/` tree),
and janet-kit's own
`defn/typed` docstring-hoisting behaviour (0.1.2+) is what attaches
this library's docstrings.

Refresh this copy when janet-kit's `jennystats/typed.janet` changes
upstream, in the same commit that consumes the change.

## `libs/janet-num.c` and `libs/install-janet-num.sh`

Copied from [janet-num](https://github.com/jennystats/janet-num) at
0.1.2 (AGPLv3). The C source is the pinned kernel library the
`jennystats/num` module loads through core FFI; the installer verifies
its own recorded sha256 before compiling. The suite needs the built
library (`libjn.so`) whenever the `janet-num` path is exercised, and a
clean machine has no way to obtain it without the source, so it is
vendored here rather than fetched by the package manager.

Refresh these two files only alongside an explicit janet-num release
sync: the installer carries the digest, so a partial or hand-edited copy
fails at its checksum instead of building the wrong kernels.
