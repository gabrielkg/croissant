# Croissant

Fast and simple haplotype counting.

## Installation

Only ```sbt``` is required; it handles the rest. Once it is installed:

```
git clone https://github.com/gabrielkg/croissant.git
cd croissant
sbt "runMain Croissant --help"
```

Will download dependencies, compile and run. You should eventually see something like:

```
[info] welcome to sbt 1.5.5 (Oracle Corporation Java 1.8.0_202)
[info] loading settings for project croissant-build from plugins.sbt ...
[info] loading project definition from /group/home/gk0a/temp/croissant/project
[info] loading settings for project croissant from build.sbt ...
[info] set current project to croissant (in build file:/group/home/gk0a/temp/croissant/)
[info] running Croissant --help
Croissant 0.3 (c) 2015-2021 Gabriel Keeble-Gagnere, Agriculture Victoria
  -a, --alignment  <arg>   Sorted BAM alignment file
  -c, --cov-only           Only report coverage (no haplotype information)
  -m, --mp                 Aligned data is mate pair
  -v, --verbose
  -w, --window  <arg>      Number of mismatches around target to consider in
                           haplotype calculation
  -h, --help               Show help message
      --version            Show version of this program

```
