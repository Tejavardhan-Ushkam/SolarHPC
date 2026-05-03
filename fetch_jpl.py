#!/usr/bin/env python3
"""
fetch_jpl.py  --  SolarHPC Project
Query NASA JPL Horizons API for initial conditions of all 10 solar
system bodies at epoch J2000.0 (JD 2451545.0).

Outputs:
  data/initial_conditions.dat   (10 bodies, 8 columns each)
  data/body_masses.dat          (10 masses, one per line)

Run:  python3 fetch_jpl.py
"""

import requests
import os
import sys
import re

# ── Configuration ────────────────────────────────────────────
API_URL  = "https://ssd.jpl.nasa.gov/api/horizons.api"
J2000_JD = "2451545.0"

# Body order MUST match solar_constants.f90 BODY INDEX MAP
BODIES = [
    ("Sun",     "10",  1.989000e+30),
    ("Mercury", "199", 3.285000e+23),
    ("Venus",   "299", 4.867000e+24),
    ("Earth",   "399", 5.972000e+24),
    ("Moon",    "301", 7.342000e+22),
    ("Mars",    "499", 6.417000e+23),
    ("Jupiter", "599", 1.898000e+27),
    ("Saturn",  "699", 5.683000e+26),
    ("Uranus",  "799", 8.681000e+25),
    ("Neptune", "899", 1.024000e+26),
]

PARAMS_BASE = {
    "format":      "text",
    "OBJ_DATA":    "NO",
    "MAKE_EPHEM":  "YES",
    "TABLE_TYPE":  "VECTORS",
    "CENTER":      "500@0",       # solar system barycentre
    "OUT_UNITS":   "AU-D",        # AU and AU/day
    "VEC_TABLE":   "2",
    "VEC_LABELS":  "YES",
    "CSV_FORMAT":  "YES",
    "TLIST":       J2000_JD,
}

def fetch_body(name, naif_id):
    """Fetch state vector for one body. Returns (x,y,z,vx,vy,vz) in AU/AU_day."""
    params = dict(PARAMS_BASE)
    params["COMMAND"] = naif_id

    try:
        r = requests.get(API_URL, params=params, timeout=30)
    except requests.exceptions.RequestException as e:
        print(f"  [ERROR] Network error fetching {name}: {e}")
        sys.exit(1)

    if r.status_code != 200:
        print(f"  [ERROR] HTTP {r.status_code} for {name}")
        sys.exit(1)

    text = r.text

    # Extract the block between $$SOE and $$EOE
    soe = text.find("$$SOE")
    eoe = text.find("$$EOE")
    if soe == -1 or eoe == -1:
        print(f"  [ERROR] $$SOE/$$EOE markers not found for {name}")
        print(f"  Response excerpt:\n{text[:500]}")
        sys.exit(1)

    block = text[soe+5:eoe].strip()
    lines = [l.strip() for l in block.split("\n") if l.strip()]

    # CSV format: line 0 = JD, CALDATE, ...
    #             line 1 = X, Y, Z, VX, VY, VZ, ...
    # Some API versions put X,Y,Z on line 1 and VX,VY,VZ on line 2
    # We extract all floats from lines 1 and 2 and take first 6.
    # Extract floats from ALL lines inside the block
    all_floats = []
    for line in lines:
        nums = re.findall(r'[-+]?\d+\.\d+[Ee]?[-+]?\d*', line)
        all_floats.extend([float(n) for n in nums])

    # Expect: JD + 6 values = at least 7 floats
    if len(all_floats) < 7:
        print(f"  [ERROR] Could not parse 6 state vector components for {name}")
        print(f"  Parsed floats: {all_floats}")
        print(f"  Raw block:\n{block[:300]}")
        sys.exit(1)

    # Skip JD, take next 6 values
    x, y, z, vx, vy, vz = all_floats[1:7]
    return x, y, z, vx, vy, vz


def main():
    os.makedirs("data", exist_ok=True)

    results = []
    for name, naif_id, mass in BODIES:
        print(f"  Fetching {name} (NAIF {naif_id}) ...", end=" ", flush=True)
        x, y, z, vx, vy, vz = fetch_body(name, naif_id)
        results.append((name, mass, x, y, z, vx, vy, vz))
        print(f"x={x:+.4f}  y={y:+.4f}  z={z:+.4f} AU  [OK]")

    # Write initial_conditions.dat
    ic_path = "data/initial_conditions.dat"
    with open(ic_path, "w") as f:
        f.write("# BODY_NAME MASS_KG X_AU Y_AU Z_AU VX_AU_PER_DAY VY_AU_PER_DAY VZ_AU_PER_DAY\n")
        for name, mass, x, y, z, vx, vy, vz in results:
            f.write(f"{name:<10s} {mass:.6e} {x:+.12e} {y:+.12e} {z:+.12e} "
                    f"{vx:+.12e} {vy:+.12e} {vz:+.12e}\n")

    # Write body_masses.dat
    mass_path = "data/body_masses.dat"
    with open(mass_path, "w") as f:
        for name, mass, *_ in results:
            f.write(f"{mass:.6e}   ! {name}\n")

    print()
    print(f"  [OK] {ic_path}   ({len(results)} bodies)")
    print(f"  [OK] {mass_path}")
    print()
    print("  Verify:  wc -l data/initial_conditions.dat   # should print 11")


if __name__ == "__main__":
    main()

