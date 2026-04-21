#!/usr/bin/env python3
"""
plot_orbits.py -- SolarHPC
Generates orbital trajectory figures.
Run LOCALLY after downloading results/ from cluster.

Produces:
  plots/output/fig_orbits_2d.pdf/.png
  plots/output/fig_orbits_3d.pdf/.png
  plots/output/fig_spacecraft_trajectory.pdf/.png  (if data exists)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

try:
    plt.style.use('seaborn-v0_8-paper')
except:
    try: plt.style.use('seaborn-paper')
    except: pass

plt.rcParams.update({
    'font.family':'serif','font.size':8,'axes.labelsize':9,
    'figure.dpi':300,'savefig.dpi':300,'savefig.bbox':'tight',
})

OUT = 'plots/output'
os.makedirs(OUT, exist_ok=True)

BODY_COLORS = {
    'Sun':'gold','Mercury':'#8B8B83','Venus':'#FFA500',
    'Earth':'#1E90FF','Moon':'#C0C0C0','Mars':'#CD3700',
    'Jupiter':'#DAA520','Saturn':'#D2B48C','Uranus':'#7FFFD4','Neptune':'#4169E1',
}

def load_orbits():
    path = 'results/orbits/orbit_trajectories.csv'
    if not os.path.exists(path):
        print(f'  [WARN] {path} not found -- generating synthetic orbits for layout')
        # Synthetic circular orbits for testing
        bodies = ['Sun','Mercury','Venus','Earth','Moon','Mars',
                  'Jupiter','Saturn','Uranus','Neptune']
        radii  = [0, 0.387, 0.723, 1.0, 1.00257, 1.524, 5.2, 9.58, 19.2, 30.1]
        rows   = []
        t_vals = np.linspace(0, 365*2, 500)
        for bi, (b, r) in enumerate(zip(bodies, radii)):
            period = r**1.5 if r > 0 else 1
            for t in t_vals:
                ang = 2*np.pi * t / (period*365.25)
                rows.append({'body_name':b,'julian_day':2451545+t,
                              'x_au':r*np.cos(ang),'y_au':r*np.sin(ang),'z_au':0.0})
        return pd.DataFrame(rows)
    return pd.read_csv(path)


def plot_2d(df):
    fig = plt.figure(figsize=(7.0, 3.5))

    # Left panel: full solar system
    ax1 = fig.add_subplot(1, 2, 1)
    # Right panel: inner solar system
    ax2 = fig.add_subplot(1, 2, 2)

    bodies = df['body_name'].unique()

    for body in bodies:
        sub = df[df['body_name'] == body]
        c   = BODY_COLORS.get(body, '#333333')
        # Full view
        ax1.plot(sub['x_au'], sub['y_au'], '-', color=c, lw=0.6, alpha=0.7)
        ax1.plot(sub['x_au'].iloc[-1], sub['y_au'].iloc[-1], 'o', color=c, ms=4)
        # Inner view (only bodies within 2 AU)
        r_max = np.sqrt(sub['x_au']**2 + sub['y_au']**2).max()
        if r_max < 2.1:
            ax2.plot(sub['x_au'], sub['y_au'], '-', color=c, lw=0.8, alpha=0.8)
            ax2.plot(sub['x_au'].iloc[-1], sub['y_au'].iloc[-1], 'o', color=c, ms=5)
            ax2.annotate(body, (sub['x_au'].iloc[-1], sub['y_au'].iloc[-1]),
                         fontsize=6, xytext=(4,4), textcoords='offset points')

    # Sun marker
    ax1.plot(0, 0, '*', color='gold', ms=10, zorder=5)
    ax2.plot(0, 0, '*', color='gold', ms=10, zorder=5)

    ax1.set_xlim(-32, 32); ax1.set_ylim(-32, 32)
    ax1.set_xlabel('x [AU]'); ax1.set_ylabel('y [AU]')
    ax1.set_title('Solar system (ecliptic view)', fontsize=9)
    ax1.set_aspect('equal'); ax1.grid(True, alpha=0.2)

    ax2.set_xlim(-2.2, 2.2); ax2.set_ylim(-2.2, 2.2)
    ax2.set_xlabel('x [AU]'); ax2.set_ylabel('y [AU]')
    ax2.set_title('Inner solar system', fontsize=9)
    ax2.set_aspect('equal'); ax2.grid(True, alpha=0.2)

    fig.tight_layout()
    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_orbits_2d.{ext}')
    plt.close(fig)
    print('  [OK] fig_orbits_2d')


def plot_3d(df):
    fig = plt.figure(figsize=(5, 4.5))
    ax  = fig.add_subplot(111, projection='3d')

    bodies = df['body_name'].unique()
    for body in bodies:
        sub = df[df['body_name'] == body]
        c   = BODY_COLORS.get(body, '#333333')
        ax.plot(sub['x_au'], sub['y_au'], sub.get('z_au', sub['y_au']*0),
                '-', color=c, lw=0.6, alpha=0.7, label=body)
        ax.scatter([sub['x_au'].iloc[-1]], [sub['y_au'].iloc[-1]],
                   [sub.get('z_au', sub['y_au']*0).iloc[-1]],
                   color=c, s=15, zorder=5)

    ax.scatter([0],[0],[0], color='gold', s=80, zorder=10, marker='*')
    ax.set_xlabel('x [AU]', labelpad=3)
    ax.set_ylabel('y [AU]', labelpad=3)
    ax.set_zlabel('z [AU]', labelpad=3)
    ax.set_title('3D orbital view', fontsize=9)
    ax.view_init(elev=30, azim=45)
    ax.legend(fontsize=5, loc='upper left', ncol=2)

    fig.tight_layout()
    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_orbits_3d.{ext}')
    plt.close(fig)
    print('  [OK] fig_orbits_3d')


def plot_spacecraft(df_orbits):
    sc_path  = 'results/missions/spacecraft_path.csv'
    mis_path = 'results/missions/mission_report.csv'

    if not os.path.exists(sc_path):
        print('  [SKIP] fig_spacecraft_trajectory (no spacecraft_path.csv)')
        return

    sc  = pd.read_csv(sc_path)
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    # Inner solar system background
    bodies_inner = ['Earth','Moon','Mars','Venus']
    for body in bodies_inner:
        sub = df_orbits[df_orbits['body_name'] == body]
        if len(sub) == 0: continue
        c = BODY_COLORS.get(body, '#333333')
        ax.plot(sub['x_au'], sub['y_au'], '-', color=c, lw=0.5, alpha=0.4)

    # Sun
    ax.plot(0, 0, '*', color='gold', ms=10, zorder=5)

    # Spacecraft path
    ax.plot(sc['x_au'], sc['y_au'], '--', color='darkorange', lw=1.2,
            label='Spacecraft', zorder=4)
    ax.plot(sc['x_au'].iloc[0],  sc['y_au'].iloc[0],  '^', color='green',
            ms=8, label='Launch', zorder=6)
    ax.plot(sc['x_au'].iloc[-1], sc['y_au'].iloc[-1], '*', color='red',
            ms=10, label='Arrival', zorder=6)

    if os.path.exists(mis_path):
        try:
            mr = pd.read_csv(mis_path).iloc[-1]
            info = (f"$\\Delta v$={mr['dv_total_SI']/1000:.1f} km/s\n"
                    f"Transfer: {mr['transfer_time_days']:.0f} d\n"
                    f"Fuel: {mr['fuel_mass_kg']:.0f} kg")
            ax.text(0.03, 0.97, info, transform=ax.transAxes, fontsize=7,
                    va='top', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
        except Exception:
            pass

    ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5)
    ax.set_xlabel('x [AU]'); ax.set_ylabel('y [AU]')
    ax.set_title('Spacecraft trajectory', fontsize=9)
    ax.set_aspect('equal'); ax.grid(True, alpha=0.2)
    ax.legend(fontsize=7, loc='lower right')
    fig.tight_layout()

    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_spacecraft_trajectory.{ext}')
    plt.close(fig)
    print('  [OK] fig_spacecraft_trajectory')


def main():
    print()
    print('  Generating orbit figures...')
    print()
    df = load_orbits()
    plot_2d(df)
    plot_3d(df)
    plot_spacecraft(df)
    print()
    print(f'  Orbit figures saved to {OUT}/')
    print()

if __name__ == '__main__':
    main()
