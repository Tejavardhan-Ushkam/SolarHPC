#!/usr/bin/env python3
"""
fetch_reference_ephemeris.py  --  SolarHPC
Fetches JPL Horizons state vectors at 5 future epochs for validation.
Produces data/reference_ephemeris.csv used by ephemeris_compare.cpp
and plot_validation.py to quantify positional accuracy over time.

Run ONCE before the simulation:  python3 fetch_reference_ephemeris.py
"""

import requests, os, sys, re

API_URL  = "https://ssd.jpl.nasa.gov/api/horizons.api"
J2000_JD = 2451545.0

# Validation epochs: 1, 5, 10, 50 years after J2000
EPOCHS = [
    J2000_JD + 365.25 * 1,
    J2000_JD + 365.25 * 5,
    J2000_JD + 365.25 * 10,
    J2000_JD + 365.25 * 50,
]

BODIES = [
    ("Sun",     "10"),
    ("Mercury", "199"),
    ("Venus",   "299"),
    ("Earth",   "399"),
    ("Moon",    "301"),
    ("Mars",    "499"),
    ("Jupiter", "599"),
    ("Saturn",  "699"),
    ("Uranus",  "799"),
    ("Neptune", "899"),
]

def fetch_at_epoch(naif_id, jd):
    params = {
        "format": "text", "COMMAND": naif_id,
        "OBJ_DATA": "NO", "MAKE_EPHEM": "YES",
        "TABLE_TYPE": "VECTORS", "CENTER": "500@0",
        "OUT_UNITS": "AU-D", "VEC_TABLE": "2",
        "VEC_LABELS": "YES", "CSV_FORMAT": "YES",
        "TLIST": f"{jd:.4f}",
    }
    try:
        r = requests.get(API_URL, params=params, timeout=30)
    except requests.exceptions.RequestException as e:
        print(f"    [ERROR] Network: {e}"); sys.exit(1)
    if r.status_code != 200:
        print(f"    [ERROR] HTTP {r.status_code}"); sys.exit(1)
    text = r.text
    soe = text.find("$$SOE"); eoe = text.find("$$EOE")
    if soe == -1: return None
    block = text[soe+5:eoe].strip()
    lines = [l.strip() for l in block.split("\n") if l.strip()]
    floats = []
    for line in lines:
        floats.extend([float(n) for n in
                    re.findall(r'[-+]?\d+\.\d+[Ee]?[-+]?\d*', line)])

    # Expect at least 7 floats: JD + 6 values
    if len(floats) < 7:
        return None

    # Skip JD, take next 6 values
    return floats[1:7]# x,y,z,vx,vy,vz

def main():
    os.makedirs("data", exist_ok=True)
    outpath = "data/reference_ephemeris.csv"
    rows = []

    print(f"Fetching reference ephemeris at {len(EPOCHS)} epochs x {len(BODIES)} bodies...")
    print()

    for ei, jd in enumerate(EPOCHS):
        yr = 2000 + (jd - J2000_JD) / 365.25
        print(f"  Epoch {ei+1}/{len(EPOCHS)}: JD {jd:.1f}  (year {yr:.1f})")
        for name, naif_id in BODIES:
            sv = fetch_at_epoch(naif_id, jd)
            if sv:
                rows.append((name, jd) + tuple(sv))
                print(f"    [OK] {name}")
            else:
                print(f"    [WARN] {name} -- skipped")

    with open(outpath, "w") as f:
        f.write("body_name,epoch_jd,x_au,y_au,z_au,vx_au_day,vy_au_day,vz_au_day\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]:.4f},{r[2]:+.12e},{r[3]:+.12e},{r[4]:+.12e},"
                    f"{r[5]:+.12e},{r[6]:+.12e},{r[7]:+.12e}\n")

    print()
    print(f"  [OK] {outpath}  ({len(rows)} rows)")
    print(f"  Verify: wc -l {outpath}  # should be {len(rows)+1}")

if __name__ == "__main__":
    main()

