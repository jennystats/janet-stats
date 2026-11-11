# Vendored C source

`libs/janet-num.c` and `libs/install-janet-num.sh` are copied from
[janet-num](https://github.com/jennystats/janet-num) at 0.1.2
(AGPLv3). The C source is the pinned kernel library the
`jennystats/num` module loads through core FFI; the installer verifies
its own recorded sha256 before compiling. The suite needs the built
library (`libjn.so`) whenever the `janet-num` path is exercised, and a
clean machine has no way to obtain it without the source, so it is
carried here rather than fetched by the package manager.

Refresh these two files only alongside an explicit janet-num release
sync: the installer carries the digest, so a partial or hand-edited copy
fails at its checksum instead of building the wrong kernels.

The Janet modules themselves are not vendored: `typed` and `kit` come
from the declared [janet-kit](https://github.com/jennystats/janet-kit)
dependency and `num` from [janet-num](https://github.com/jennystats/janet-num).
