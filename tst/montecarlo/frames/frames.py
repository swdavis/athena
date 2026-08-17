#! /usr/bin/env python

"""
Frame-consistency tests for the Monte Carlo radiation moments.

These check identities that must hold exactly, so they need no analytic solution, no
convergence study and no reference data -- only that two numbers agree.  They exist
because the moment machinery has two independent frame conventions (the lab basis the
moments are accumulated in, and the comoving basis they are boosted to) which can drift
apart silently: a bug that put a spurious Doppler factor into the general relativistic
comoving moments went unnoticed precisely because no test compared the two.

  identity  With the fluid at rest the boost to the comoving frame is the identity, so
            moments_com must equal moments_lab zone by zone.  Run on both pushers.

  pushers   In flat spacetime the general pusher, geodesic integration and ConstructTetrad
            must reproduce the legacy cartesian answer exactly, boosted or not.  Any
            disagreement is a bug in the general-pusher machinery, not a modelling
            difference.  This is the same claim tst/montecarlo/minkowski makes, checked
            here directly against the cartesian run rather than against a blackbody.

Deliberately uses tab output rather than hdf5 so it runs without an HDF5 build.

Usage:
  python frames.py [--path /path/to/athena] [--nphot N] [--tol 1e-12]
"""

import argparse
import glob
import os
import sys
import numpy as np

KB = 1.380649e-16
TGAS = 1.0e6
DENS = 1.0e-10
PGAS_CODE = 1.0e-6


def write_athinput(path, kind, vel, nphot, iseed):
    """kind is 'cart' (legacy pusher) or 'mink' (general pusher, GR machinery)."""
    general = "true" if kind == "mink" else "false"
    pid = "mciso" if kind == "cart" else "mcmink"
    o = ["<job>", "problem_id = " + pid, "",
         "<output1>", "file_type = tab", "variable = mclab", "id = lab", "dt = 1.0", "",
         "<output2>", "file_type = tab", "variable = mccom", "id = com", "dt = 1.0", "",
         "<time>", "cfl_number = 0.3", "nlim = 1", "tlim = 1.0", "", "<mesh>"]
    for d in (1, 2, 3):
        o += ["nx{0} = 8".format(d),
              "x{0}min = -5.0e10".format(d), "x{0}max = 5.0e10".format(d),
              "ix{0}_bc = periodic".format(d), "ox{0}_bc = periodic".format(d),
              "ix{0}_mc_bc = periodic".format(d), "ox{0}_mc_bc = periodic".format(d), ""]
    o += ["<meshblock>", "nx1 = 4", "nx2 = 4", "nx3 = 4", "",
          "<hydro>", "gamma = 1.666666666666667", "",
          "<montecarlo>",
          "nphot = {0:d}".format(nphot), "iseed = {0:d}".format(iseed),
          "general_pusher = " + general,
          "scattering = none", "emission = freefree", "absorption = freefree",
          "polarized = false", "checkmove = 1000000", "stepsize = 1.0e-2",
          "varystep = true", "abs_method = weight",
          # boosts stay on even at zero velocity: that is what makes the comoving
          # transform run and the identity meaningful rather than vacuous.
          "boosts = true", "", "<problem>"]
    if kind == "cart":
        o += ["temp = {0:e}".format(TGAS), "dens = {0:e}".format(DENS),
              "velocity = {0:e}".format(vel), "constdens = true"]
    else:
        o += ["dens_code = 1.0", "pgas_code = {0:e}".format(PGAS_CODE),
              "rho_cgs = {0:e}".format(DENS),
              "tgas_cgs = {0:e}".format(TGAS / PGAS_CODE),
              "l_cgs = 1.0", "velocity = {0:e}".format(vel)]
    o += ["emin = {0:e}".format(0.001 * KB * TGAS),
          "emax = {0:e}".format(50.0 * KB * TGAS), ""]
    open(path, "w").write("\n".join(o))


def read_moments(rundir, tag):
    """Concatenate the per-block tab files for one output id."""
    rows, hdr = [], None
    for f in sorted(glob.glob(os.path.join(rundir, "*.{0}.*.tab".format(tag)))):
        for line in open(f):
            if line.startswith("#"):
                if "x1v" in line:
                    hdr = line.split()[1:]
                continue
            rows.append([float(x) for x in line.split()])
    if not rows:
        raise RuntimeError("no {0} output found in {1}".format(tag, rundir))
    return np.array(rows), hdr


def moment_columns(arr, hdr):
    """Return the moment columns only, dropping index and coordinate columns."""
    keep = [i for i, c in enumerate(hdr)
            if c.startswith(("Ermc", "Frmc", "Prmc"))]
    return arr[:, keep]


def max_rel(a, b):
    den = np.maximum(np.abs(a), np.abs(b))
    den[den == 0.0] = 1.0
    return float(np.max(np.abs(a - b) / den))


def build_and_run(src, kind, vel, workdir, nphot, iseed):
    import subprocess
    athena_path = os.path.join(src, "bin")
    prob, coord, extra = (("mc_isoth", "cartesian", []) if kind == "cart"
                          else ("mc_isoth_mink", "minkowski", ["-g"]))
    # configure.py rewrites src/defs.hpp, which the Makefile does not track as a
    # dependency, so a stale build would silently link objects from the previous
    # coordinate system.
    subprocess.call(["make", "clean"], cwd=src,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.check_call(["python", "configure.py", "--prob=" + prob,
                           "--coord=" + coord, "-mc"] + extra,
                          cwd=src, stdout=subprocess.DEVNULL)
    subprocess.check_call(["make", "-j8"], cwd=src,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.makedirs(workdir, exist_ok=True)
    inp = os.path.join(workdir, "athinput.frames")
    write_athinput(inp, kind, vel, nphot, iseed)
    with open(os.path.join(workdir, "run.log"), "w") as log:
        subprocess.check_call([os.path.join(athena_path, "athena"), "-i", inp],
                              cwd=workdir, stdout=log, stderr=subprocess.STDOUT)
    return workdir


def main(**kwargs):
    path = kwargs.pop("path")
    nphot = kwargs.pop("nphot")
    tol = kwargs.pop("tol")
    iseed = kwargs.pop("iseed")
    work = kwargs.pop("workdir")
    src = os.path.abspath(path)

    failures = []
    runs = {}
    for kind in ("cart", "mink"):
        for vel, label in ((0.0, "zero"), (0.1, "boost")):
            name = "{0}_{1}".format(kind, label)
            runs[name] = build_and_run(src, kind, vel,
                                       os.path.join(work, name), nphot, iseed)

    print("{0:<34}{1:>14}   {2}".format("check", "max rel diff", "result"))
    print("-" * 66)

    # identity: fluid at rest => comoving moments equal lab moments
    for name in ("cart_zero", "mink_zero"):
        lab, hl = read_moments(runs[name], "lab")
        com, hc = read_moments(runs[name], "com")
        d = max_rel(moment_columns(lab, hl), moment_columns(com, hc))
        ok = d < tol
        failures += [] if ok else ["identity/" + name]
        print("{0:<34}{1:>14.3e}   {2}".format(
            "identity  " + name + "  com==lab", d, "PASS" if ok else "FAIL"))

    # pushers: general pusher in flat spacetime reproduces the legacy cartesian answer
    for label in ("zero", "boost"):
        for tag in ("lab", "com"):
            a, ha = read_moments(runs["cart_" + label], tag)
            b, hb = read_moments(runs["mink_" + label], tag)
            d = max_rel(moment_columns(a, ha), moment_columns(b, hb))
            ok = d < tol
            failures += [] if ok else ["pushers/" + label + "/" + tag]
            print("{0:<34}{1:>14.3e}   {2}".format(
                "pushers   " + label + "  " + tag + "  cart==mink", d,
                "PASS" if ok else "FAIL"))

    print("-" * 66)
    if failures:
        print("FAILED: " + ", ".join(failures))
        return 1
    print("all frame-consistency checks passed")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--path", default="../../..", help="path to athena root")
    p.add_argument("--workdir", default="/tmp/mc_frames", help="scratch run directory")
    p.add_argument("--nphot", type=int, default=200000, help="photons per run")
    p.add_argument("--iseed", type=int, default=17, help="random seed")
    p.add_argument("--tol", type=float, default=1.0e-12, help="relative tolerance")
    sys.exit(main(**vars(p.parse_args())))
