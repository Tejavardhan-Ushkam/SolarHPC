#!/usr/bin/env python3
"""
plot_eclipse_timeline.py -- SolarHPC
Eclipse prediction timeline and Saros cycle figures.
Run LOCALLY after downloading results/ from cluster.

Produces:
  plots/output/fig_eclipse_timeline.pdf/.png
  plots/output/fig_saros_pattern.pdf/.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta

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

J2000_JD       = 2451545.0
SAROS_DAYS     = 6585.3211
TARGET_2024_JD = 2460409.0    # 2024-Apr-08 total solar eclipse

def jd_to_year(jd):
    """Convert Julian Day to decimal year (approximate)."""
    return 2000.0 + (jd - J2000_JD) / 365.25

def load_eclipses():
    path = 'results/eclipses/eclipse_predictions.csv'
    if not os.path.exists(path):
        print(f'  [WARN] {path} not found -- using synthetic data')
        np.random.seed(42)
        rows = []
        # Generate ~30 solar, ~40 lunar eclipses over 2 years
        jd = J2000_JD
        for _ in range(35):
            jd += np.random.uniform(25, 35)
            rows.append({
                'julian_day':    jd,
                'gregorian_date': '2000-01-01',
                'eclipse_type':  0,
                'event_class':   np.random.choice([1,2,3], p=[0.5,0.3,0.2]),
                'umbra_fraction': np.random.uniform(0.3, 1.0),
                'duration_minutes': np.random.uniform(60, 180),
            })
        jd = J2000_JD
        for _ in range(45):
            jd += np.random.uniform(22, 30)
            rows.append({
                'julian_day':    jd,
                'gregorian_date': '2000-01-01',
                'eclipse_type':  1,
                'event_class':   np.random.choice([1,2], p=[0.6,0.4]),
                'umbra_fraction': np.random.uniform(0.2, 1.0),
                'duration_minutes': np.random.uniform(40, 100),
            })
        # Add the 2024 eclipse for validation demo
        rows.append({
            'julian_day': TARGET_2024_JD,
            'gregorian_date': '2024-04-08',
            'eclipse_type': 0, 'event_class': 2,
            'umbra_fraction': 1.0, 'duration_minutes': 268,
        })
        return pd.DataFrame(rows).sort_values('julian_day')
    return pd.read_csv(path)


def plot_timeline(df):
    solar = df[df['eclipse_type'] == 0].copy()
    lunar = df[df['eclipse_type'] == 1].copy()

    solar['year'] = solar['julian_day'].apply(jd_to_year)
    lunar['year'] = lunar['julian_day'].apply(jd_to_year)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(7.0, 4.5),
                                           gridspec_kw={'height_ratios':[3,2]})

    # ── Top panel: full timeline ──────────────────────────────
    # Solar eclipses: above x-axis
    sc_top = ax_top.scatter(solar['year'], solar['umbra_fraction'],
                             c='#d62728', s=40,
                             marker='^', label='Solar eclipse', zorder=4, alpha=0.8)
    # Lunar eclipses: below x-axis
    sc_bot = ax_top.scatter(lunar['year'], -lunar['umbra_fraction']*0.6,
                             c='#1f77b4', s=30,
                             marker='o', label='Lunar eclipse', zorder=4, alpha=0.7)

    # Mark 2024-04-08 validation point if in range
    yr_2024 = jd_to_year(TARGET_2024_JD)
    yr_min  = solar['year'].min() if len(solar) else 2000
    yr_max  = solar['year'].max() if len(solar) else 2002
    if yr_min <= yr_2024 <= yr_max:
        ax_top.scatter([yr_2024], [1.05], c='gold', s=120, marker='*',
                       zorder=6, label='2024-Apr-08 (validation)')
        ax_top.annotate('2024-04-08\nTotal solar\n(validation)',
                        xy=(yr_2024, 1.0), xytext=(yr_2024+0.1, 1.2),
                        fontsize=6, ha='left',
                        arrowprops=dict(arrowstyle='->', lw=0.8),
                        bbox=dict(boxstyle='round,pad=0.2', fc='lightyellow', alpha=0.8))

    ax_top.axhline(0, color='k', lw=0.5, alpha=0.4)
    ax_top.set_xlim(yr_min - 0.1, yr_max + 0.1)
    ax_top.set_ylim(-0.8, 1.4)
    ax_top.set_xlabel('Year')
    ax_top.set_ylabel('Umbra fraction')
    ax_top.set_title(f'Eclipse predictions ({len(solar)} solar, {len(lunar)} lunar)', fontsize=9)
    ax_top.legend(fontsize=7, loc='upper right', ncol=3)
    ax_top.grid(True, alpha=0.2)
    ax_top.text(0.01, 0.97, 'Solar (above) / Lunar (below)',
                transform=ax_top.transAxes, fontsize=7, va='top', style='italic')

    # ── Bottom panel: zoom 5 years ────────────────────────────
    zoom_start = yr_min
    zoom_end   = min(yr_min + 5, yr_max)
    s_zoom = solar[(solar['year'] >= zoom_start) & (solar['year'] <= zoom_end)]
    l_zoom = lunar[(lunar['year'] >= zoom_start) & (lunar['year'] <= zoom_end)]

    ax_bot.vlines(s_zoom['year'], 0, s_zoom['umbra_fraction'],
                  colors='#d62728', lw=1.5, alpha=0.7, label='Solar')
    ax_bot.vlines(l_zoom['year'], 0, -l_zoom['umbra_fraction']*0.6,
                  colors='#1f77b4', lw=1.5, alpha=0.7, label='Lunar')
    ax_bot.axhline(0, color='k', lw=0.5)
    ax_bot.set_xlim(zoom_start - 0.05, zoom_end + 0.05)
    ax_bot.set_xlabel('Year')
    ax_bot.set_ylabel('Umbra fraction')
    ax_bot.set_title(f'First 5-year zoom ({zoom_start:.1f}–{zoom_end:.1f})', fontsize=9)
    ax_bot.legend(fontsize=7, loc='upper right')
    ax_bot.grid(True, alpha=0.2)

    fig.tight_layout()
    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_eclipse_timeline.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_eclipse_timeline  '
          f'({len(solar)} solar, {len(lunar)} lunar eclipses)')


def plot_saros(df):
    """Show Saros cycle repetition: find eclipse pairs ~6585 days apart."""
    jds     = df['julian_day'].values
    etypes  = df['eclipse_type'].values
    n       = len(jds)

    pairs   = []   # (idx_early, idx_late, interval)
    for i in range(n):
        for j in range(i+1, n):
            interval = jds[j] - jds[i]
            if interval > SAROS_DAYS * 1.5:
                break
            if (abs(interval - SAROS_DAYS) < 5.0 and etypes[i] == etypes[j]):
                pairs.append((i, j, interval))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

    # ── Left: Saros pair timeline ─────────────────────────────
    years = [jd_to_year(jd) for jd in jds]
    solar_yr = [years[i] for i in range(n) if etypes[i] == 0]
    lunar_yr = [years[i] for i in range(n) if etypes[i] == 1]

    ax1.scatter(solar_yr, [1]*len(solar_yr), c='#d62728', s=12,
                marker='^', label='Solar', zorder=3)
    ax1.scatter(lunar_yr, [0]*len(lunar_yr), c='#1f77b4', s=12,
                marker='o', label='Lunar', zorder=3)

    # Draw Saros connections
    for (i, j, _) in pairs[:30]:   # limit lines to avoid clutter
        color = '#d62728' if etypes[i] == 0 else '#1f77b4'
        y_val = 1 if etypes[i] == 0 else 0
        ax1.plot([years[i], years[j]], [y_val, y_val],
                 '-', color=color, alpha=0.3, lw=0.8)

    ax1.set_yticks([0, 1])
    ax1.set_yticklabels(['Lunar', 'Solar'], fontsize=8)
    ax1.set_xlabel('Year')
    ax1.set_title(f'Saros pairs found: {len(pairs)}', fontsize=9)
    ax1.legend(fontsize=7)
    ax1.grid(True, alpha=0.2, axis='x')

    # ── Right: histogram of Saros intervals ───────────────────
    if pairs:
        intervals = [p[2] for p in pairs]
        ax2.hist(intervals, bins=20, color='#2ca02c', alpha=0.75, edgecolor='white')
        ax2.axvline(SAROS_DAYS, color='red', lw=1.5, ls='--',
                    label=f'Saros = {SAROS_DAYS:.1f} d')
        mean_int = np.mean(intervals)
        ax2.axvline(mean_int, color='orange', lw=1.2, ls=':',
                    label=f'Mean = {mean_int:.1f} d')
        ax2.set_xlabel('Inter-eclipse interval [days]')
        ax2.set_ylabel('Count')
        ax2.set_title('Saros period distribution', fontsize=9)
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.2)
        err = abs(mean_int - SAROS_DAYS)
        ax2.text(0.03, 0.95, f'Error from nominal: {err:.2f} d',
                 transform=ax2.transAxes, fontsize=7, va='top',
                 bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))
    else:
        ax2.text(0.5, 0.5, 'No Saros pairs found\n(run longer simulation)',
                 ha='center', va='center', transform=ax2.transAxes, fontsize=9)

    fig.tight_layout()
    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_saros_pattern.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_saros_pattern  ({len(pairs)} Saros pairs identified)')


def main():
    print()
    print('  Generating eclipse figures...')
    print()
    df = load_eclipses()
    plot_timeline(df)
    plot_saros(df)
    print()
    print(f'  Eclipse figures saved to {OUT}/')
    print()

if __name__ == '__main__':
    main()
