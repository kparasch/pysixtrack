import pickle
import time

import numpy as np

from pyoptics import madlang, optics
import pysixtrack

# see sps/madx/a001_track_thin.madx
mad = madlang.open("madx/SPS_Q20_thin.seq")
mad.acta_31637.volt = 4.5
mad.acta_31637.lag = 0.5
out, rest = mad.sps.expand_struct()


out, rest = mad.sps.expand_struct(pysixtrack.convert)
elems = list(zip(*out))[1]
sps = pysixtrack.Block(elems)


pickle.dump(sps, open("sps.pickle", "w"))


npart = 10


def check_el(p):
    madt = optics.open("sps/madx/track.obs0002.p0001")
    out1 = [madt[l][0] for l in "x px y py y t pt".split()]
    out2 = [getattr(p, l) for l in "x px y py y tau pt".split()]
    diff = 0
    for a, b in zip(out1, out2):
        diff += (a - b) ** 2
        print("%24.17e %24.17e %24.17e" % (a, b, a - b))
    print(np.sqrt(diff))


def trackn(n):
    p = pysixtrack.Particles(
        x=1e-3,
        px=0.0,
        y=-0.5e-3,
        py=0.0,
        tau=0.74,
        pt=0.0,
        e0=26.01692438e9,
        m0=0.93827205e9,
    )
    for iel, el in enumerate(sps.elems[:n]):
        print(iel, el)
        el.track(p)
        print(
            "%12.9f %12.9f %12.9f %12.9f %12.9f"
            % (p.s, p.x, p.px, p.tau, p.pt)
        )
        if abs(p.x) > 1:
            break
    print((out[n - 1][0]))
    return p


p = trackn(27113)
check_el(p)


def track_turn(n):
    p = pysixtrack.Particles(
        x=1e-3,
        px=np.zeros(5000),
        y=-0.5e-3,
        py=np.zeros(5000),
        tau=0.74,
        pt=np.zeros(5000),
        e0=26.01692438e9,
        m0=0.93827205e9,
    )
    out = []
    before = time.time()
    for i in range(n):
        out.append(p.copy())
        sps.track(p)
        now = time.time()
        print((i, now - before))
        before = now
    return out


tt = track_turn(10)
