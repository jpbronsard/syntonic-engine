# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Baton-tip position of MAESTRO 1's right arm from its six joint angles (numpy only).

The chain (body offsets, joint axes and anchors, tip site) is exported from the MuJoCo model of the robot and
evaluated the way MuJoCo's mj_kinematics does; checked against MuJoCo to 1e-12 m. The turntable J0 is held at zero:
it turns the whole arm, which does not change the distance between the actual and the planned tip.
"""
import os
import numpy as np

_Z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'right_arm_chain.npz'))
_BODIES, _JOINTS, _SITE = _Z['bodies'], _Z['joints'], _Z['site_pos']


def _qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def _rot(q):
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def tip(q6):
    """q6: (6,) right-arm joint angles RJ1..RJ6 (rad). Returns the baton-tip position (m) in the arm's base frame."""
    xpos = np.zeros(3); xquat = np.array([1.0, 0.0, 0.0, 0.0])
    for i, row in enumerate(_BODIES):
        xpos = xpos + _rot(xquat) @ row[:3]; xquat = _qmul(xquat, row[3:7])
        for jr in _JOINTS[_JOINTS[:, 0] == i]:
            idx, qpos0, axis, jpos = int(jr[1]), jr[2], jr[3:6], jr[6:9]
            ang = (q6[idx] if idx >= 0 else qpos0) - qpos0
            anchor = xpos + _rot(xquat) @ jpos
            n = np.linalg.norm(axis); s = np.sin(ang / 2)
            xquat = _qmul(xquat, np.array([np.cos(ang / 2), *(s * axis / n)]))
            xquat = xquat / np.linalg.norm(xquat)
            xpos = anchor - _rot(xquat) @ jpos
    return xpos + _rot(xquat) @ _SITE
