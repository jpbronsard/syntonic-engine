#!/bin/sh
# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
#
# Cortex-M7 cycle model of the loop law's gain recomputation, before and after the two speed-ups of 0.3.2. Compiles
# for two joints as Teensyduino compiles for a Teensy 4.x (arm-none-eabi-g++ -O2 -mcpu=cortex-m7 -mfpu=fpv5-d16
# -mfloat-abi=hard -ffp-contract=off) and runs llvm-mca with LLVM's Cortex-M7 scheduling model on the innermost loops:
#   before  the candidate loop of 0.3.1's work_one_ (cpp/SyntonicLoopLaw.h at commit 5073c27), as
#           bench/port_2026-11-02/m7/cycle_model.sh measured it
#   L       0.3.2's kernel loop_block_L (|S|^2 from L; the same arithmetic as 0.3.1)
#   quad    0.3.2's kernel loop_block_quad (|S|^2 by the quadratic form)
#   row     the loop that fills a Kp row's P (once per row)
# per (candidate, frequency) point (per frequency for row). A model, not a measurement: examples/teensy/LoopLawTiming.ino
# measures on the board. Writes cycles.json next to this script. Needs git, arm-none-eabi-g++ and llvm-mca.
#   sh bench/scan_2026-11-02/m7/cycle_model.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd -P)
REPO=$(cd "$HERE/../../.." && pwd -P)
T=$(mktemp -d)
F="-mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb -O2 -std=gnu++14 -fno-exceptions -fno-rtti -ffp-contract=off"
MCA="llvm-mca -mtriple=thumbv7em-none-eabihf -mcpu=cortex-m7"
mkdir "$T/old"
git -C "$REPO" show 5073c27:cpp/SyntonicLoopLaw.h > "$T/old/SyntonicLoopLaw.h"
cp "$REPO/cpp/SyntonicEngine.h" "$T/old/"
cat > "$T/old.cpp" <<'CPP'
#include "SyntonicLoopLaw.h"
static syntonic::SyntonicLoopLaw<2> law;
extern "C" bool law_begin(const double* q0, const double* w0) { syntonic::LoopLawConfig c; c.dead = 0.0025; return law.begin(q0, w0, c); }
extern "C" const double* law_step(const double* q, const double* qd, const double* a, const double* b, bool s) { return law.step(q, qd, a, b, s); }
CPP
cat > "$T/new.cpp" <<'CPP'
#include "SyntonicLoopLaw.h"
using namespace syntonic;
extern "C" void kern_L(double kp, double ki, const double* iW, const double* Gr, const double* Gi,
                       const double (*PhiW)[LOOP_NW], double* J) { loop_block_L<2>(kp, ki, iW, Gr, Gi, PhiW, 0, 50, J); }
extern "C" void kern_quad(double ki, const double* P, const double* Q, const double* R, const double (*PhiW)[LOOP_NW],
                          double* J) { loop_block_quad<2>(ki, P, Q, R, PhiW, 0, 50, J); }
extern "C" void kern_row(double kp, const double* Gr, const double* Gi, double* P) {   // the expression of row_()
  for (int w = 0; w < LOOP_NW; w++) P[w] = 1.0 + kp * (2.0 * Gr[w] + kp * (Gr[w] * Gr[w] + Gi[w] * Gi[w]));
}
CPP
arm-none-eabi-g++ $F -I"$T/old" -S "$T/old.cpp" -o "$T/old.s"
arm-none-eabi-g++ $F -I"$REPO/cpp" -S "$T/new.cpp" -o "$T/new.s"
# the innermost loop of FUNCTION holding PATTERN: a label, straight-line code, a branch back to the label
loop_of() {   # loop_of FILE FUNCTION PATTERN -> loop body on stdout
    awk -v fn="$2" '$0 ~ "^"fn":" {p = 1} p && /\.size/ {exit} p' "$1" | awk -v pat="$3" '
        /^\.L[0-9]+:/ {lab = substr($1, 1, length($1) - 1); n = 0; next}
        lab != "" && $1 !~ /^[@.]/ && NF {body[n++] = $0; if ($2 == lab && $1 ~ /^b/) {
            s = ""; for (i = 0; i < n; i++) s = s body[i] "\n"; if (s ~ pat) {printf "%s", s; exit} }}'
}
loop_of "$T/old.s" "_ZN8syntonic15SyntonicLoopLawILi2ELi16EE9work_one_Ev" 'vdiv\.f64' > "$T/before.s"
loop_of "$T/new.s" "kern_L" 'vdiv\.f64' > "$T/L.s"
loop_of "$T/new.s" "kern_quad" 'vdiv\.f64' > "$T/quad.s"
loop_of "$T/new.s" "kern_row" 'vstr|vstmia' > "$T/row.s"
cyc() { $MCA -iterations=1000 "$1" 2>/dev/null | awk '/Total Cycles/ {printf "%.1f", $3 / 1000}'; }
B=$(cyc "$T/before.s"); L=$(cyc "$T/L.s"); Q=$(cyc "$T/quad.s"); R=$(cyc "$T/row.s")
for k in before L quad row; do
    echo "$k: $(grep -c . "$T/$k.s") instructions, $(cyc "$T/$k.s") cycles per iteration"
done
for op in vmul vadd vdiv; do
    : > "$T/op.s"; for r in 3 4 5 6 7 8 9 10 11 12; do printf '\t%s.f64\td%d, d13, d14\n' $op $r >> "$T/op.s"; done
    echo "  model throughput of $op.f64: $(cyc "$T/op.s" | awk '{printf "%.1f", $1 / 10}') cycles"
done
printf '{"before": %s, "L": %s, "quad": %s, "row": %s, "toolchain": "%s", "model": "%s"}\n' "$B" "$L" "$Q" "$R" \
    "$(arm-none-eabi-g++ --version | head -1)" "$(llvm-mca --version | grep -i 'llvm version' | head -1 | sed 's/^ *//')" \
    > "$HERE/cycles.json"
cat "$HERE/cycles.json"
rm -rf "$T"
