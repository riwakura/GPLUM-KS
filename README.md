# GPLUM
GPLUM is an N-body simulation code for planetary system formation using the particle–particle particle–tree (P3T) scheme.

The current version is 2.5.

The user guide in Japanese is available [here](doc/UsersGuide_japanese.pdf).
An English user guide is in preparation.

- 2026/10/9 Ryutaro Iwakura added KS-Hermite scheme for hard binary integration.

## Differences from the original code

RI implemented a KS-regularized Hermite scheme based on Mikkola & Aarseth (1998, New Astronomy, 3, 309; https://ui.adsabs.harvard.edu/abs/1998NewA....3..309M/abstract) for binaries within hard clusters.
New implemented shceme reduces the cost of integrating hard binaries in planetary disks.
The speedup depends on the system being simulated.
RI only tested for planetary disks.

## How to use
Binary integration is enabled by the following setting in `src/Makefile`:

```makefile
use_KS = yes
```

To build GPLUM without KS integration, comment out that line (`#use_KS = yes`).

With `use_KS = yes`, the binary's relative motion includes perturbations from the central star and other members of the hard cluster.
To retain only the central star's perturbation, also uncomment:

```makefile
ks_solar_only = yes
```

This is disabled by default.

KS profiling is disabled by default. To collect integration statistics and write
`binary_profile_rank*.csv`, set `ks_profile = 1` in the parameter file.
Set `ks_profile = 0` to disable profiling; no rebuild is required.
