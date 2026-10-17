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
