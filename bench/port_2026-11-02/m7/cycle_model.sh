#!/bin/sh
# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
#
# Cortex-M7 cycle model of the loop law's gain recomputation (0.3.1). Compiles cpp/SyntonicLoopLaw.h for two joints as
# Teensyduino compiles for a Teensy 4.x (arm-none-eabi-g++ -O2 -mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard), takes
# the candidate loop (one frequency point of one candidate, both joints) from the assembly, and runs llvm-mca with
# LLVM's Cortex-M7 scheduling model on it; also the model's throughput of single double-precision operations, and the
# loop of a proposed quadratic form (not in the law). A model, not a measurement: examples/teensy/LoopLawTiming.ino
# measures on the board. Needs arm-none-eabi-g++ and llvm-mca.   sh bench/port_2026-11-02/m7/cycle_model.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd -P)
CPP=$(cd "$HERE/../../../cpp" && pwd -P)
T=$(mktemp -d)
F="-mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb -O2 -std=gnu++14 -fno-exceptions -fno-rtti -ffp-contract=off"
MCA="llvm-mca -mtriple=thumbv7em-none-eabihf -mcpu=cortex-m7"
cat > "$T/law.cpp" <<'CPP'
#include "SyntonicLoopLaw.h"
static syntonic::SyntonicLoopLaw<2> law;
extern "C" bool law_begin(const double* q0, const double* w0) { syntonic::LoopLawConfig c; c.dead = 0.0025; return law.begin(q0, w0, c); }
extern "C" const double* law_step(const double* q, const double* qd, const double* a, const double* b, bool s) { return law.step(q, qd, a, b, s); }
CPP
cat > "$T/quad.cpp" <<'CPP'
extern "C" void candq(const double* __restrict P, const double* __restrict Q, const double* __restrict R,
                      const double* __restrict P0, const double* __restrict P1, double ki, double* J) {
  double J0 = 0, J1 = 0;
  for (int w = 0; w < 1500; w++) { double s2 = 1.0 / (P[w] + ki * (Q[w] + ki * R[w])); J0 += P0[w] * s2; J1 += P1[w] * s2; }
  J[0] = J0; J[1] = J1;
}
CPP
arm-none-eabi-g++ $F -I"$CPP" -S "$T/law.cpp" -o "$T/law.s"
arm-none-eabi-g++ $F -S "$T/quad.cpp" -o "$T/quad.s"
# the innermost loop of work_one_ that holds the division: a label, straight-line code, a branch back to the label
loop_of() {   # loop_of FILE FUNCTION -> loop body on stdout
    awk -v fn="$2" '$0 ~ "^"fn":" {p = 1} p && /\.size/ {exit} p' "$1" | awk '
        /^\.L[0-9]+:/ {lab = substr($1, 1, length($1) - 1); n = 0; next}
        lab != "" && $1 !~ /^[@.]/ && NF {body[n++] = $0; if ($2 == lab && $1 ~ /^b/) {
            s = ""; for (i = 0; i < n; i++) s = s body[i] "\n"; if (s ~ /vdiv\.f64/) {printf "%s", s; exit} }}'
}
loop_of "$T/law.s" "_ZN8syntonic15SyntonicLoopLawILi2ELi16EE9work_one_Ev" > "$T/loop.s"
loop_of "$T/quad.s" "candq" > "$T/quadloop.s"
cyc() { $MCA -iterations=1000 "$1" 2>/dev/null | awk '/Total Cycles/ {printf "%.1f", $3 / 1000}'; }
L=$(cyc "$T/loop.s"); Q=$(cyc "$T/quadloop.s")
echo "loop law, two joints: $(grep -c . "$T/loop.s") instructions and $L cycles per (candidate, frequency) point"
for op in vmul vadd vmla vdiv; do
    : > "$T/op.s"; for r in 3 4 5 6 7 8 9 10 11 12; do printf '\t%s.f64\td%d, d13, d14\n' $op $r >> "$T/op.s"; done
    echo "  model throughput of $op.f64: $(cyc "$T/op.s" | awk '{printf "%.1f", $1 / 10}') cycles"
done
echo "proposed quadratic form (not in the law): $Q cycles per point"
awk -v c="$L" 'BEGIN { for (n = 878; n <= 1246; n += 368) printf "recomputation, %d robust candidates: %.0f Mcycles = %.0f ms at 600 MHz (spectrum not included)\n", n, n * 1500 * c / 1e6, n * 1500 * c / 600e3 }'
rm -rf "$T"
