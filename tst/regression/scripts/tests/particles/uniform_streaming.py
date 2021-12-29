# Test for uniform streaming between gas and particles.

# Modules
from collections import namedtuple
import math
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data # noqa


def prepare(**kwargs):
    """Configure and make the executable. """

    athena.configure("p", "mpi", "hdf5", "h5double", prob="uniform_streaming", **kwargs)
    athena.make()


def run(**kwargs):
    """Run the executable. """
    import os
    import subprocess

    # Construct a list of arguments to the executable.
    arguments = []

    # Run the executable.
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'], 2,
                  'particles/athinput.uniform_streaming', arguments)


def analyze():
    """Analyze the output and determine if the test passes. """
    import athena_read
    from glob import glob

    # Read the input file.
    athinput = athena_read.athinput("../../inputs/particles//athinput.uniform_streaming")
    input_mesh = athinput["mesh"]
    input_particles = athinput["particles"]
    input_problem = athinput["problem"]

    # Define the base name.
    base = "bin/" + athinput["job"]["problem_id"]

    # Get the initial conditions.
    xlen = float(input_mesh["x1max"]) - float(input_mesh["x1min"])
    ylen = float(input_mesh["x2max"]) - float(input_mesh["x2min"])
    zlen = float(input_mesh["x3max"]) - float(input_mesh["x3min"])

    dtog = float(input_problem["dtog"])
    ux0 = float(input_problem["ux0"]) if "ux0" in input_problem else 0
    uy0 = float(input_problem["uy0"]) if "uy0" in input_problem else 0
    uz0 = float(input_problem["uz0"]) if "uz0" in input_problem else 0
    vpx0 = float(input_problem["vpx0"]) if "vpx0" in input_problem else 0
    vpy0 = float(input_problem["vpy0"]) if "vpy0" in input_problem else 0
    vpz0 = float(input_problem["vpz0"]) if "vpz0" in input_problem else 0

    taus = float(input_particles["taus0"])
    backreaction = False
    if "backreaction" in input_particles:
        if input_particles["backreaction"] == "true":
            backreaction = True

    # Get the initial particle positions.
    t, dp = athena_read.particles(base + ".pout.00000.dat")
    xp0, yp0, zp0 = np.copy(dp.xp), np.copy(dp.yp), np.copy(dp.zp)
    npar = len(xp0)

    # Initialization
    t = []
    dxpavg, dxpmin, dxpmax = [], [], []
    dypavg, dypmin, dypmax = [], [], []
    dzpavg, dzpmin, dzpmax = [], [], []
    vpxavg, vpxmin, vpxmax = [], [], []
    vpyavg, vpymin, vpymax = [], [], []
    vpzavg, vpzmin, vpzmax = [], [], []
    uxavg, uxmin, uxmax = [], [], []
    uyavg, uymin, uymax = [], [], []
    uzavg, uzmin, uzmax = [], [], []
    sx = np.zeros((npar,), dtype=int)
    sy = np.zeros((npar,), dtype=int)
    sz = np.zeros((npar,), dtype=int)
    xpold, ypold, zpold = np.copy(xp0), np.copy(yp0), np.copy(zp0)

    # Collect particle data.
    for fname in sorted(glob(base + ".pout.*.dat")):
        time, dp = athena_read.particles(fname)
        t.append(time)

        # Process particle positions.
        sx[np.where(xpold - dp.xp > 0.5 * xlen)] += 1
        sx[np.where(dp.xp - xpold > 0.5 * xlen)] -= 1
        sy[np.where(ypold - dp.yp > 0.5 * ylen)] += 1
        sy[np.where(dp.yp - ypold > 0.5 * ylen)] -= 1
        sz[np.where(zpold - dp.zp > 0.5 * zlen)] += 1
        sz[np.where(dp.zp - zpold > 0.5 * zlen)] -= 1

        dxp = (dp.xp - xp0) + sx * xlen
        dyp = (dp.yp - yp0) + sy * ylen
        dzp = (dp.zp - zp0) + sz * zlen

        dxpavg.append(dxp.mean())
        dxpmin.append(dxp.min())
        dxpmax.append(dxp.max())

        dypavg.append(dyp.mean())
        dypmin.append(dyp.min())
        dypmax.append(dyp.max())

        dzpavg.append(dzp.mean())
        dzpmin.append(dzp.min())
        dzpmax.append(dzp.max())

        xpold = np.copy(dp.xp)
        ypold = np.copy(dp.yp)
        zpold = np.copy(dp.zp)

        # Process particle velocities.
        vpxavg.append(dp.vpx.mean())
        vpxmin.append(dp.vpx.min())
        vpxmax.append(dp.vpx.max())
        vpyavg.append(dp.vpy.mean())
        vpymin.append(dp.vpy.min())
        vpymax.append(dp.vpy.max())
        vpzavg.append(dp.vpz.mean())
        vpzmin.append(dp.vpz.min())
        vpzmax.append(dp.vpz.max())

    t = np.array(t)
    dxpavg, dxpmin, dxpmax = np.array(dxpavg), np.array(dxpmin), np.array(dxpmax)
    dypavg, dypmin, dypmax = np.array(dypavg), np.array(dypmin), np.array(dypmax)
    dzpavg, dzpmin, dzpmax = np.array(dzpavg), np.array(dzpmin), np.array(dzpmax)
    vpxavg, vpxmin, vpxmax = np.array(vpxavg), np.array(vpxmin), np.array(vpxmax)
    vpyavg, vpymin, vpymax = np.array(vpyavg), np.array(vpymin), np.array(vpymax)
    vpzavg, vpzmin, vpzmax = np.array(vpzavg), np.array(vpzmin), np.array(vpzmax)

    # Collect gas data.
    for fname in sorted(glob(base + ".out1.[0-9][0-9][0-9][0-9][0-9].athdf")):
        dg = athena_read.athdf(fname)
        dg = namedtuple("HDF5", dg.keys())(**dg)

        # Process gas velocities.
        uxavg.append(dg.vel1.mean())
        uxmin.append(dg.vel1.min())
        uxmax.append(dg.vel1.max())
        uyavg.append(dg.vel2.mean())
        uymin.append(dg.vel2.min())
        uymax.append(dg.vel2.max())
        uzavg.append(dg.vel3.mean())
        uzmin.append(dg.vel3.min())
        uzmax.append(dg.vel3.max())

    uxavg, uxmin, uxmax = np.array(uxavg), np.array(uxmin), np.array(uxmax)
    uyavg, uymin, uymax = np.array(uyavg), np.array(uymin), np.array(uymax)
    uzavg, uzmin, uzmax = np.array(uzavg), np.array(uzmin), np.array(uzmax)

    # Find the analytical solution.
    us1 = UniformStreaming(taus, 0, dtog, ux0, vpx0, backreaction=backreaction)
    us2 = UniformStreaming(taus, 0, dtog, uy0, vpy0, backreaction=backreaction)
    us3 = UniformStreaming(taus, 0, dtog, uz0, vpz0, backreaction=backreaction)
    te = np.linspace(t[0], t[-1], 201)
    dxpe = us1.displacement(te)
    dype = us2.displacement(te)
    dzpe = us3.displacement(te)
    uxe, vpxe = us1.velocities(te)
    uye, vpye = us2.velocities(te)
    uze, vpze = us3.velocities(te)

    # Define the norm.
    def norm(u, v):
        s = 0
        for a, b in zip(u, v):
            s += (a - b)**2
        return math.sqrt(s)

    # Find the absolute errors.
    err_drp = norm([dxpavg[-1], dypavg[-1], dzpavg[-1]],
                   [dxpe[-1], dype[-1], dzpe[-1]])
    err_v = norm([vpxavg[-1], vpyavg[-1], vpzavg[-1]],
                 [vpxe[-1], vpye[-1], vpze[-1]])
    err_u = norm([uxavg[-1], uyavg[-1], uzavg[-1]],
                 [uxe[-1], uye[-1], uze[-1]])
    print("\nAbsolute Errors in:\n")
    print("\tParticle displacement = {}".format(err_drp))
    print("\tParticle velocity     = {}".format(err_v))
    print("\tGas velocity          = {}\n".format(err_u))

    # Evaluate the uniformity.
    ddrp = norm([dxpmax[-1], dypmax[-1], dzpmax[-1]],
                [dxpmin[-1], dypmin[-1], dzpmin[-1]])
    dv = norm([vpxmax[-1], vpymax[-1], vpzmax[-1]],
              [vpxmin[-1], vpymin[-1], vpzmin[-1]])
    du = norm([uxmax[-1], uymax[-1], uzmax[-1]],
              [uxmin[-1], uymin[-1], uzmin[-1]])
    print("\nData Range in:\n")
    print("\tParticle displacement = {}".format(ddrp))
    print("\tParticle velocity     = {}".format(dv))
    print("\tGas velocity          = {}\n".format(du))

    # Detect anomalies.
    ok = True
    if err_drp > 2.00E-6:
        print("*** Too much error in particle displacement")
        ok = False
    if err_v > 4.80E-6:
        print("*** Too much error in particle velocity")
        ok = False
    if err_u > 9.60E-7:
        print("*** Too much error in gas velocity")
        ok = False
    if ddrp > 3E-15:
        print("*** Not uniform in particle displacement")
        ok = False
    if dv > 2E-15:
        print("*** Not uniform in particle velocity")
        ok = False
    if du > 4E-15:
        print("*** Not uniform in gas velocity")
        ok = False

    return ok


# ======================================================================
# Class definition for uniform streaming motion between gas and solid
# particles.
# ======================================================================
class UniformStreaming:
    """Object for uniform streaming motion between gas and solid
    particles.
    """
# ----------------------------------------------------------------------
    def __init__(self, ts, g, epsilon, u0, v0, backreaction=True):
        """Initializes the object.

        Positional Arguments:
            ts
                Stopping time between gas and particles.
            g
                External acceleration on the gas.
            epsilon
                Solid-to-gas density ratio.
            u0
                Initial velocity of the gas.
            v0
                Initial velocity of the particles.
        """
        # Record the input parameters.
        self.ts = ts
        self.g = g
        self.epsilon = epsilon
        self.u0 = u0
        self.v0 = v0
        self.br = backreaction
        # Find the center-of-mass velocity.
        self.ucm = (u0 + epsilon * v0) / (1 + epsilon)

# ----------------------------------------------------------------------
    def velocities(self, t):
        """Find the velocities of the gas and the particles as a
        function of time.

        Positional Argument:
            t
                Time, which can be a list or a numpy array.

        Returned Values:
            Velocities of the gas and the particles, respectively, at
            time(s) t.
        """
        t = np.array(t)
        if self.br:
            m = 1 + self.epsilon
            a = self.g * t / m
            q = self.g * self.ts / m**2
            s = np.exp(-m / self.ts * t)
            r = 1 - s
            u = self.u0 * s + a + (self.ucm + q) * r
            v = self.v0 * s + a + (self.ucm - q) * r
        else:
            r = np.exp(-t / self.ts)
            u = self.u0 * np.ones_like(t)
            v = self.v0 * r + self.u0 * (1 - r)

        return u, v

# ----------------------------------------------------------------------
    def displacement(self, t):
        """Find the displacement of a particle function of time.

        Positional Argument:
            t
                Time, which can be a list or a numpy array.

        Returned Value:
            Displacement(s) of a particle at time(s) t.
        """
        t = np.array(t)
        if self.br:
            a = 1 + self.epsilon
            r = 1 - np.exp(-a / self.ts * t)
            a = 1 / a
            w = self.ucm - a**2 * self.g * self.ts
            dxp = a * self.ts * (self.v0 - w) * r + w * t + 0.5 * a * self.g * t**2
        else:
            dxp = self.u0 * t + (self.v0 - self.u0) * self.ts * (1 - np.exp(-t / self.ts))

        return dxp
# ======================================================================
