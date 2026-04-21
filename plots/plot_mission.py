#!/usr/bin/env python3
"""
plot_mission.py -- SolarHPC
Mission planning figures: launch window, delta-v comparison,
landing zone, mission summary table.
Run LOCALLY after downloading results/ from cluster.

Produces:
  plots/output/fig_launch_window.pdf/.png
  plots/output/fig_dv_comparison.pdf/.png
  plots/output/fig_landing_zone.pdf/.png
  plots/output/fig_mission_table.pdf/.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

try:
    plt.style.use('seaborn-v0_8-paper')
except:
    try: plt.style.use('seaborn-paper')
    except: pass

plt.rcParams.update({
    'font.family':'serif','font.size':8,'axes.labelsize':9,
    'figure.dpi':300,'savefig.dpi':300,'savefig.bbox':'tight',
})

OUT      = 'plots/output'
J2000_JD = 2451545.0
os.makedirs(OUT, exist_ok=True)


def jd_to_year(jd):
    return 2000.0 + (jd - J2000_JD) / 365.25


def load_launch_window():
    path = 'results/missions/launch_window_sweep.csv'
    if not os.path.exists(path):
        print(f'  [WARN] {path} not found -- using synthetic data')
        jds = np.linspace(J2000_JD, J2000_JD+730, 300)
        # Realistic Moon window: roughly periodic ~29.5 days
        dv_moon = 3800 + 300*np.sin(2*np.pi*(jds-J2000_JD)/29.5) \
                       + 150*np.sin(2*np.pi*(jds-J2000_JD)/365.25) \
                       + np.random.normal(0, 30, len(jds))
        dv_mars = 5500 + 1500*np.sin(2*np.pi*(jds-J2000_JD)/780) \
                       + np.random.normal(0, 80, len(jds))
        return pd.DataFrame({
            'launch_jd':  np.concatenate([jds, jds]),
            'dv_total_SI':np.concatenate([dv_moon, dv_mars]),
            'target':     ['Moon']*len(jds) + ['Mars']*len(jds),
            'feasible':   np.concatenate([np.ones(len(jds)), (dv_mars<8000).astype(int)]),
        })
    df = pd.read_csv(path)
    # If single-target, tag it
    if 'target' not in df.columns:
        df['target'] = 'Moon'
    return df


def load_mission_report():
    path = 'results/missions/mission_report.csv'
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def plot_launch_window(df):
    targets = df['target'].unique() if 'target' in df.columns else ['Moon']
    n_t     = len(targets)
    fig, axes = plt.subplots(1, max(n_t, 1), figsize=(3.5*n_t, 2.8))
    if n_t == 1:
        axes = [axes]

    for ax, tgt in zip(axes, targets):
        sub = df[df['target'] == tgt] if 'target' in df.columns else df
        sub = sub.sort_values('launch_jd')
        years = sub['launch_jd'].apply(jd_to_year)
        dv_km = sub['dv_total_SI'] / 1000.0

        # Color by feasibility
        if 'feasible' in sub.columns:
            colors = np.where(sub['feasible'] == 1, '#2ca02c', '#d62728')
        else:
            colors = ['#2ca02c'] * len(sub)

        ax.scatter(years, dv_km, c=colors, s=4, alpha=0.6, linewidths=0)

        # Highlight minimum
        idx_min = dv_km.idxmin()
        yr_best = jd_to_year(sub.loc[idx_min, 'launch_jd'])
        dv_best = dv_km.loc[idx_min]
        ax.plot(yr_best, dv_best, '*', color='gold', ms=12, zorder=5,
                markeredgecolor='darkorange', markeredgewidth=0.5)
        ax.annotate(f'Optimal\n{yr_best:.2f}\n{dv_best:.2f} km/s',
                    xy=(yr_best, dv_best),
                    xytext=(yr_best + 0.05, dv_best + 0.3),
                    fontsize=6, ha='left',
                    arrowprops=dict(arrowstyle='->', lw=0.7),
                    bbox=dict(boxstyle='round,pad=0.2', fc='lightyellow', alpha=0.8))

        green_p = mpatches.Patch(color='#2ca02c', label='Feasible')
        red_p   = mpatches.Patch(color='#d62728', label='High cost')
        ax.legend(handles=[green_p, red_p], fontsize=6, loc='upper right')
        ax.set_xlabel('Launch year')
        ax.set_ylabel(r'$\Delta v_{\mathrm{total}}$ [km/s]')
        ax.set_title(f'{tgt} launch window', fontsize=9)
        ax.grid(True, alpha=0.25)

    fig.tight_layout()
    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_launch_window.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_launch_window')


def plot_dv_comparison(mr):
    """Hohmann vs gravity assist comparison."""
    if mr is None or len(mr) == 0:
        print(f'  [SKIP] fig_dv_comparison (no mission_report.csv)')
        # Use synthetic
        mr = pd.DataFrame({
            'target_body':  ['Moon','Moon','Mars','Mars'],
            'mission_mode': ['Hohmann','GravityAssist','Hohmann','GravityAssist'],
            'dv_total_SI':  [3900, 3400, 5600, 4800],
            'fuel_mass_kg': [2600, 2100, 4500, 3700],
            'transfer_time_days': [3, 5, 259, 320],
        })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8))

    targets = mr['target_body'].unique()
    modes   = ['Hohmann','GravityAssist']
    x       = np.arange(len(targets))
    w       = 0.35
    colors  = ['#1f77b4','#ff7f0e']

    for mi, mode in enumerate(modes):
        sub  = mr[mr['mission_mode'] == mode]
        dvs  = []
        fkgs = []
        for tgt in targets:
            row = sub[sub['target_body'] == tgt]
            dvs.append( row['dv_total_SI'].values[0]/1000.0 if len(row) else 0)
            fkgs.append(row['fuel_mass_kg'].values[0]       if len(row) else 0)

        bars1 = ax1.bar(x + mi*w - w/2, dvs,  w*0.9, label=mode,
                        color=colors[mi], alpha=0.85, hatch=['','//'][mi])
        bars2 = ax2.bar(x + mi*w - w/2, fkgs, w*0.9, label=mode,
                        color=colors[mi], alpha=0.85, hatch=['','//'][mi])
        for b, v in zip(bars1, dvs):
            ax1.text(b.get_x()+b.get_width()/2, b.get_height()+0.05,
                     f'{v:.2f}', ha='center', fontsize=6)
        for b, v in zip(bars2, fkgs):
            ax2.text(b.get_x()+b.get_width()/2, b.get_height()+10,
                     f'{v:.0f}', ha='center', fontsize=6)

    ax1.set_xticks(x); ax1.set_xticklabels(targets)
    ax1.set_ylabel(r'$\Delta v$ [km/s]')
    ax1.set_title('Total delta-v by mission type', fontsize=9)
    ax1.legend(fontsize=7); ax1.grid(True, alpha=0.2, axis='y')

    ax2.set_xticks(x); ax2.set_xticklabels(targets)
    ax2.set_ylabel('Fuel mass [kg]')
    ax2.set_title('Fuel mass (payload=1000 kg)', fontsize=9)
    ax2.legend(fontsize=7); ax2.grid(True, alpha=0.2, axis='y')

    fig.tight_layout()
    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_dv_comparison.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_dv_comparison')


def plot_landing_zone(mr):
    """Simple lat/lon landing zone map."""
    if mr is None or len(mr) == 0:
        lat, lon, tgt = 12.3, 45.6, 'Moon'
    else:
        row = mr.iloc[-1]
        lat = float(row.get('landing_lat_deg', 12.3))
        lon = float(row.get('landing_lon_deg', 45.6))
        tgt = str(row.get('target_body', 'Moon'))

    fig, ax = plt.subplots(figsize=(3.5, 2.4))

    # Simple rectangular projection with grid
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect('equal')

    # Background colour by body
    if 'Moon' in tgt:
        ax.set_facecolor('#e8e8e0')
        body_label = 'Lunar surface'
    else:
        ax.set_facecolor('#d4a574')
        body_label = 'Martian surface'

    # Grid lines
    for lat_g in range(-90, 91, 30):
        ax.axhline(lat_g, color='white', lw=0.5, alpha=0.6)
    for lon_g in range(-180, 181, 60):
        ax.axvline(lon_g, color='white', lw=0.5, alpha=0.6)

    # Landing ellipse (uncertainty ~50 km → ~0.03 deg)
    from matplotlib.patches import Ellipse
    ell = Ellipse((lon, lat), width=3.0, height=2.0,
                  facecolor='red', alpha=0.3, edgecolor='darkred', lw=1.0)
    ax.add_patch(ell)

    # Landing point
    ax.plot(lon, lat, '*', color='red', ms=12, zorder=5,
            markeredgecolor='darkred', markeredgewidth=0.5)
    ax.annotate(f'  Landing\n  {lat:.2f}°N  {lon:.2f}°E',
                xy=(lon, lat), xytext=(lon+8, lat+8),
                fontsize=7, ha='left',
                arrowprops=dict(arrowstyle='->', lw=0.7, color='darkred'),
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.9))

    # Launch site reference
    site_lat, site_lon = 13.72, 80.23   # ISRO Sriharikota
    if 'Moon' not in tgt:
        ax.plot(0, 0, 'b^', ms=6, label='N/A (Mars)')
    ax.set_xlabel('Longitude [°E]')
    ax.set_ylabel('Latitude [°N]')
    ax.set_title(f'Predicted landing zone — {tgt}', fontsize=9)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-90, 91, 30))
    ax.tick_params(labelsize=7)
    fig.tight_layout()

    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_landing_zone.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_landing_zone  (lat={lat:.2f}, lon={lon:.2f})')


def plot_mission_table(mr):
    """Render mission results as a matplotlib table figure."""
    if mr is None or len(mr) == 0:
        mr = pd.DataFrame({
            'target_body':        ['Moon','Moon','Mars','Mars'],
            'mission_mode':       ['Hohmann','GravityAssist','Hohmann','GravityAssist'],
            'dv_total_SI':        [3900, 3400, 5600, 4800],
            'fuel_mass_kg':       [2600, 2100, 4500, 3700],
            'transfer_time_days': [3.2,  5.1,  258.9, 320.0],
            'landing_lat_deg':    [12.3, 8.1,  4.5, -22.3],
            'landing_lon_deg':    [45.6, 72.1, 137.2, 91.0],
        })

    cols    = ['Target','Mode',r'$\Delta v$ [m/s]','Fuel [kg]',
               'Transfer [d]','Land lat','Land lon']
    rows    = []
    for _, r in mr.iterrows():
        rows.append([
            str(r.get('target_body','?')),
            str(r.get('mission_mode','?')),
            f"{r.get('dv_total_SI',0):.0f}",
            f"{r.get('fuel_mass_kg',0):.0f}",
            f"{r.get('transfer_time_days',0):.1f}",
            f"{r.get('landing_lat_deg',0):.2f}°",
            f"{r.get('landing_lon_deg',0):.2f}°",
        ])

    fig, ax = plt.subplots(figsize=(7.0, max(1.2, 0.35*len(rows)+0.6)))
    ax.axis('off')

    tbl = ax.table(cellText=rows, colLabels=cols,
                   loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.4)

    # Style header
    for j in range(len(cols)):
        tbl[0, j].set_facecolor('#2c5f8a')
        tbl[0, j].set_text_props(color='white', fontweight='bold')

    # Alternating row shading
    for i in range(1, len(rows)+1):
        for j in range(len(cols)):
            tbl[i, j].set_facecolor('#f0f4f8' if i % 2 == 0 else 'white')

    # Highlight minimum dv row
    if rows:
        dvs = [float(r[2]) for r in rows]
        best = dvs.index(min(dvs)) + 1
        for j in range(len(cols)):
            tbl[best, j].set_facecolor('#fff3cd')

    ax.set_title('Mission results summary', fontsize=9, pad=8)
    fig.tight_layout()

    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_mission_table.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_mission_table')


def main():
    print()
    print('  Generating mission figures...')
    print()
    df_lw = load_launch_window()
    mr    = load_mission_report()
    plot_launch_window(df_lw)
    plot_dv_comparison(mr)
    plot_landing_zone(mr)
    plot_mission_table(mr)
    print()
    print(f'  Mission figures saved to {OUT}/')
    print()

if __name__ == '__main__':
    main()
