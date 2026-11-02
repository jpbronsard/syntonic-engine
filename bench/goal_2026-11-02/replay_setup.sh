#!/bin/sh
# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
#
# Replay helper for bench/goal_2026-11-02. The frozen scripts of this folder name absolute paths of the working session
# in which they ran. Editing them would change the sha256 that the confirmation scripts check (gate G0), so this script
# creates those paths as symbolic links instead:
#
#   /home/claude/repo/syntonic-engine           -> the repository (left alone if it is the repository already)
#   /home/claude/work/goal/confirm/flight_snap  -> DIR/flight_looped/flight_snap
#   /home/claude/work/goal/confirm/flight_dev   -> DIR/flight_looped/flight_dev
#   /home/claude/work/goal/r3                   -> DIR/flight_realistic
#   /home/claude/work/goal/flight               -> DIR/exploration_harness
#   /home/claude/work/goal/servo                -> DIR/servo/bench   (used by dev_pwm.py only)
#
# DIR is this folder by default. The run modes of the scripts append to their result files, so to replay them link a
# copy instead:   cp -r bench/goal_2026-11-02 /tmp/goal && sh bench/goal_2026-11-02/replay_setup.sh /tmp/goal
# It needs write access to /home/claude (as root, or in a container). It replaces only symbolic links, never a file or
# a directory.
#
#   sh replay_setup.sh [DIR]       make the links
#   sh replay_setup.sh --remove    remove them (symbolic links only)
set -eu
HERE=$(cd "$(dirname "$0")" && pwd -P)
REPO=$(cd "$HERE/../.." && pwd -P)
W=/home/claude/work/goal
R=/home/claude/repo/syntonic-engine
LINKS="confirm/flight_snap:flight_looped/flight_snap confirm/flight_dev:flight_looped/flight_dev r3:flight_realistic flight:exploration_harness servo:servo/bench"

if [ "${1:-}" = "--remove" ]; then
    for pair in $LINKS; do
        p="$W/${pair%%:*}"
        if [ -L "$p" ]; then rm "$p"; echo "removed $p"; fi
    done
    if [ -L "$R" ]; then rm "$R"; echo "removed $R"; fi
    exit 0
fi

DIR=$(cd "${1:-$HERE}" && pwd -P)
for d in servo/bench flight_looped/flight_snap flight_looped/flight_dev flight_realistic exploration_harness; do
    if [ ! -d "$DIR/$d" ]; then
        echo "$DIR/$d not found: DIR must be bench/goal_2026-11-02 or a copy of it" >&2; exit 1
    fi
done

FAIL=0
link() {    # link TARGET PATH
    if [ -L "$2" ]; then
        ln -sfn "$1" "$2"
    elif [ -e "$2" ]; then
        echo "$2 exists and is not a link: left as it is" >&2; FAIL=1; return 0
    else
        mkdir -p "$(dirname "$2")"; ln -s "$1" "$2"
    fi
    echo "$2 -> $1"
}

if [ -e "$R" ] && [ "$(cd "$R" && pwd -P)" = "$REPO" ]; then
    echo "$R is the repository"
else
    link "$REPO" "$R"
fi
for pair in $LINKS; do
    link "$DIR/${pair#*:}" "$W/${pair%%:*}"
done
exit $FAIL
