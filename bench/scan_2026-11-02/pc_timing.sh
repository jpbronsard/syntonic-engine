#!/bin/sh
# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
#
# The gain recomputation's time on this PC (g++ -O2, one core), before and after the speed-ups of 0.3.2: builds
# pc_timing.cpp against cpp/SyntonicLoopLaw.h of 0.3.1 (commit 5073c27) and of this tree with each combination of
# prune and quad, and runs each on the logs named (REPS replays each, default 5). Run from the folder that holds the logs:
#   sh bench/scan_2026-11-02/pc_timing.sh replay_loop_bench2.bin replay_loop_bench.bin replay_loop_two.bin
set -eu
LOGS="$*"
HERE=$(cd "$(dirname "$0")" && pwd -P)
REPO=$(cd "$HERE/../.." && pwd -P)
T=$(mktemp -d)
mkdir "$T/old"
git -C "$REPO" show 5073c27:cpp/SyntonicLoopLaw.h > "$T/old/SyntonicLoopLaw.h"
cp "$REPO/cpp/SyntonicEngine.h" "$T/old/"
F="-O2 -std=c++11 -ffp-contract=off"
g++ $F -I"$T/old" "$HERE/pc_timing.cpp" -o "$T/t_0.3.1"
for pq in "0 0" "1 0" "0 1" "1 1"; do
    set -- $pq
    g++ $F -DPRUNE=$1 -DQUAD=$2 -I"$REPO/cpp" "$HERE/pc_timing.cpp" -o "$T/t_p$1q$2"
done
for log in $LOGS; do
    echo "== $log"
    for b in 0.3.1 p0q0 p1q0 p0q1 p1q1; do printf '  %-6s ' "$b"; "$T/t_$b" "$log" "${REPS:-5}"; done
done
rm -rf "$T"
