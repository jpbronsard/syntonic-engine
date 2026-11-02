#!/bin/sh
# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
#
# The law replays on an emulated Cortex-M7 (0.3.1): builds qemu_replay.cpp with arm-none-eabi-g++ and newlib for the
# Teensy 4.x core (-O2 -mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard), twice: with -ffp-contract=off (the PC's
# arithmetic) and with GCC's default for the gnu++ dialects (fused multiply-adds), and runs both on QEMU's mps2-an500
# board with semihosting. QEMU checks the arithmetic (IEEE double, newlib's libm), not the timing.
# Needs arm-none-eabi-g++ with newlib (Ubuntu: gcc-arm-none-eabi libnewlib-arm-none-eabi libstdc++-arm-none-eabi-newlib)
# and qemu-system-arm. Run from the folder that holds the logs of rust/tools/make_law_replay_log.py:
#   sh bench/port_2026-11-02/m7/run_qemu.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd -P)
CPP=$(cd "$HERE/../../../cpp" && pwd -P)
OUT=$(mktemp -d)
F="-mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb -O2"
arm-none-eabi-gcc $F -c "$HERE/startup.c" -o "$OUT/startup.o"
for fc in off fast; do
    arm-none-eabi-g++ $F -std=gnu++14 -fno-exceptions -fno-rtti -ffp-contract=$fc -Wno-unused-result -I"$CPP" \
        -c "$HERE/qemu_replay.cpp" -o "$OUT/replay_$fc.o"
    arm-none-eabi-g++ $F --specs=rdimon.specs -nostartfiles -T "$HERE/link.ld" "$OUT/startup.o" "$OUT/replay_$fc.o" \
        -o "$OUT/replay_$fc.elf" -lm 2>&1 | grep -v "is not implemented and will always fail" || true
    echo "== Cortex-M7 (QEMU mps2-an500), newlib, -ffp-contract=$fc"
    qemu-system-arm -M mps2-an500 -cpu cortex-m7 -nographic -semihosting-config enable=on,target=native \
        -kernel "$OUT/replay_$fc.elf" || true
done
rm -rf "$OUT"
