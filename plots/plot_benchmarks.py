#!/usr/bin/env python3
"""
plot_benchmarks.py -- SolarHPC
Generates publication-quality benchmark figures for the paper.
Run LOCALLY after downloading results/ from cluster.

Produces:
  plots/output/fig_speedup_all_kernels.pdf/.png
  plots/output/fig_efficiency.pdf/.png
  plots/output/fig_absolute_time.pdf/.png
  plots/output/fig_mpi_scaling.pdf/.png
  plots/output/fig_amdahl_analysis.pdf/.png
  plots/output/fig_energy_drift.pdf/.png
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit

# ── Style ─────────────────────────────────────────────────────
try:
    plt.style.use('seaborn-v0_8-paper')
except:
    try:
        plt.style.use('seaborn-paper')
    except:
        pass

plt.rcParams.update({
    'font.family':       'serif',
    'font.size':         9,
    'axes.labelsize':    10,
    'axes.titlesize':    10,
    'legend.fontsize':   8,
    'xtick.labelsize':   8,
    'ytick.labelsize':   8,
    'lines.linewidth':   1.5,
    'lines.markersize':  5,
    'figure.dpi':        300,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
})

OUT  = 'plots/output'
os.makedirs(OUT, exist_ok=True)

KERNEL_COLORS   = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd']
KERNEL_MARKERS  = ['o','s','^','D','v']
KERNEL_STYLES   = ['-','--','-.',':','-']

def amdahl(N, f):
    return 1.0 / (f + (1.0 - f) / np.array(N, dtype=float))


def load_sweep(path='results/benchmarks/thread_sweep_results.csv'):
    if not os.path.exists(path):
        print(f'  [WARN] {path} not found -- using synthetic data for layout test')
        kernels = ['force_computation','leapfrog_step','rk4_step','eclipse_scan','trajectory_prop']
        rows = []
        for k in kernels:
            t1 = np.random.uniform(0.05, 0.5)
            for nt in [1,2,4,8,16]:
                f  = np.random.uniform(0.02, 0.15)
                t  = 1.0/(f + (1-f)/nt) * t1 * np.random.uniform(0.9,1.1)
                rows.append({'kernel_name':k,'n_threads':nt,'mean_sec':t,
                             'stddev_sec':t*0.03,'speedup':t1/t,
                             'efficiency_pct':t1/t/nt*100})
        return pd.DataFrame(rows)
    df = pd.read_csv(path)
    if 'mean_sec' not in df.columns:
        run_cols = [c for c in df.columns if c.startswith('run')]
        df['mean_sec']   = df[run_cols].mean(axis=1)
        df['stddev_sec'] = df[run_cols].std(axis=1)
    return df


def load_amdahl(path='results/benchmarks/amdahl_fit.csv'):
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def load_mpi(path='results/benchmarks/mpi_scaling.csv'):
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


# ── Figure 1: Speedup per kernel ─────────────────────────────
def plot_speedup(df, amdahl_df):
    kernels = df['kernel_name'].unique()
    n_k     = len(kernels)
    ncols   = 3
    nrows   = (n_k + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(7.0, nrows*2.2))
    axes = np.array(axes).flatten()

    ideal_N = np.array([1,2,4,8,16])

    for idx, kernel in enumerate(kernels):
        ax   = axes[idx]
        sub  = df[df['kernel_name'] == kernel].sort_values('n_threads')
        nt   = sub['n_threads'].values
        spd  = sub['speedup'].values
        std  = sub['stddev_sec'].values / sub['mean_sec'].values * spd  # propagated

        ax.plot(nt, nt, 'k--', lw=1.0, label='Ideal', alpha=0.5)
        ax.errorbar(nt, spd, yerr=std, fmt=KERNEL_MARKERS[idx%5],
                    color=KERNEL_COLORS[idx%5], capsize=3,
                    linestyle=KERNEL_STYLES[idx%5], label='Measured')

        # Amdahl fit
        if amdahl_df is not None:
            row = amdahl_df[amdahl_df['kernel_name'] == kernel]
            if not row.empty:
                f  = float(row['serial_fraction'].values[0])
                Nf = np.linspace(1, 16, 100)
                ax.plot(Nf, amdahl(Nf, f), 'r-', lw=1.0,
                        label=f'Amdahl f={f:.3f}', alpha=0.8)
        else:
            try:
                popt, _ = curve_fit(amdahl, nt, spd, p0=[0.1], bounds=(0,1))
                f = popt[0]
                Nf = np.linspace(1,16,100)
                ax.plot(Nf, amdahl(Nf, f), 'r-', lw=1.0,
                        label=f'Amdahl f={f:.3f}', alpha=0.8)
            except Exception:
                pass

        ax.set_xlim(0.5, 18)
        ax.set_xticks([1,2,4,8,16])
        ax.set_xlabel(r'$N_{\mathrm{threads}}$')
        ax.set_ylabel(r'Speedup $S(N)$')
        ax.set_title(kernel.replace('_',' '), fontsize=9)
        ax.legend(fontsize=6, loc='upper left')
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for j in range(len(kernels), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('OpenMP Speedup — All Kernels', fontsize=10, y=1.01)
    fig.tight_layout()

    for ext in ['pdf','png']:
        p = f'{OUT}/fig_speedup_all_kernels.{ext}'
        fig.savefig(p)
    plt.close(fig)
    print(f'  [OK] fig_speedup_all_kernels')


# ── Figure 2: Parallel efficiency ────────────────────────────
def plot_efficiency(df):
    kernels = df['kernel_name'].unique()
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))

    for idx, kernel in enumerate(kernels):
        sub = df[df['kernel_name'] == kernel].sort_values('n_threads')
        ax.plot(sub['n_threads'], sub['efficiency_pct'],
                color=KERNEL_COLORS[idx%5], marker=KERNEL_MARKERS[idx%5],
                linestyle=KERNEL_STYLES[idx%5], label=kernel.replace('_',' '))

    ax.axhline(100, color='k', lw=0.8, ls='--', alpha=0.5, label='100% ideal')
    ax.axhline(50,  color='gray', lw=0.8, ls=':', alpha=0.5, label='50% guideline')
    ax.set_xlim(0.5, 18)
    ax.set_xticks([1,2,4,8,16])
    ax.set_ylim(0, 115)
    ax.set_xlabel(r'$N_{\mathrm{threads}}$')
    ax.set_ylabel(r'Efficiency $E(N)$ [\%]')
    ax.legend(fontsize=6, loc='lower left', ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_efficiency.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_efficiency')


# ── Figure 3: Absolute wall time (bar chart) ─────────────────
def plot_absolute_time(df):
    kernels  = df['kernel_name'].unique()
    threads  = sorted(df['n_threads'].unique())
    n_k      = len(kernels)
    n_t      = len(threads)
    hatches  = ['', '//','xx','..','**']

    fig, ax  = plt.subplots(figsize=(7.0, 2.8))
    x        = np.arange(n_k)
    width    = 0.8 / n_t

    for ti, nt in enumerate(threads):
        sub  = df[df['n_threads'] == nt].set_index('kernel_name')
        vals = [sub.loc[k,'mean_sec']   if k in sub.index else 0.0 for k in kernels]
        errs = [sub.loc[k,'stddev_sec'] if k in sub.index else 0.0 for k in kernels]
        bars = ax.bar(x + (ti - n_t/2 + 0.5)*width, vals, width*0.9,
                      label=f'{nt} threads', hatch=hatches[ti%5],
                      yerr=errs, capsize=2, error_kw={'linewidth':0.8})

    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([k.replace('_','\n') for k in kernels], fontsize=7)
    ax.set_ylabel('Wall time [s]')
    ax.legend(fontsize=7, ncol=n_t)
    ax.grid(True, axis='y', alpha=0.3)
    fig.tight_layout()

    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_absolute_time.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_absolute_time')


# ── Figure 4: MPI scaling ─────────────────────────────────────
def plot_mpi_scaling(mpi_df):
    if mpi_df is None or len(mpi_df) == 0:
        print('  [SKIP] fig_mpi_scaling (no data)')
        return

    fig, ax1 = plt.subplots(figsize=(3.5, 2.8))
    ax2      = ax1.twinx()

    nr  = mpi_df['n_ranks'].values
    wt  = mpi_df['wall_time_sec'].values
    spd = mpi_df['speedup'].values
    eff = mpi_df['efficiency_pct'].values

    ax1.plot(nr, wt,  'b-o', lw=1.5, ms=5, label='Wall time')
    ax1.set_xlabel('MPI ranks')
    ax1.set_ylabel('Wall time [s]', color='b')
    ax1.tick_params(axis='y', colors='b')

    ax2.plot(nr, spd, 'r--s', lw=1.5, ms=5, label='Speedup')
    ax2.plot(nr, nr,  'k:',   lw=0.8, alpha=0.5, label='Ideal')
    ax2.set_ylabel(r'Speedup $S(P)$', color='r')
    ax2.tick_params(axis='y', colors='r')

    lines1, labs1 = ax1.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1+lines2, labs1+labs2, fontsize=7, loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(nr)
    fig.tight_layout()

    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_mpi_scaling.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_mpi_scaling')


# ── Figure 5: Amdahl analysis bar chart ──────────────────────
def plot_amdahl_analysis(df):
    kernels = df['kernel_name'].unique()
    f_vals  = []

    for kernel in kernels:
        sub = df[df['kernel_name'] == kernel].sort_values('n_threads')
        nt  = sub['n_threads'].values
        spd = sub['speedup'].values
        try:
            popt, _ = curve_fit(amdahl, nt, spd, p0=[0.05], bounds=(1e-4, 0.999))
            f_vals.append(popt[0])
        except Exception:
            f_vals.append(0.1)

    fig, ax1 = plt.subplots(figsize=(3.5, 2.8))
    ax2      = ax1.twinx()

    x     = np.arange(len(kernels))
    bars  = ax1.bar(x, f_vals, 0.6, color=KERNEL_COLORS[:len(kernels)], alpha=0.8)
    peak  = [1.0/max(f, 1e-6) for f in f_vals]
    ax2.plot(x, peak, 'ro--', ms=5, lw=1.2, label='Peak speedup 1/f')

    ax1.axhline(0.05, color='gray', ls=':', lw=0.8, label='5% guideline')
    ax1.set_xticks(x)
    ax1.set_xticklabels([k.replace('_','\n') for k in kernels], fontsize=7)
    ax1.set_ylabel('Serial fraction $f$')
    ax2.set_ylabel('Theoretical peak speedup')
    ax1.set_ylim(0, max(f_vals)*1.4)
    ax2.set_ylim(0, min(max(peak)*1.4 if peak else 20, 200))

    for xi, (fv, pk) in enumerate(zip(f_vals, peak)):
        ax1.text(xi, fv+0.002, f'{fv:.3f}', ha='center', fontsize=6)
        ax2.text(xi, pk+0.5,   f'{pk:.1f}x', ha='center', fontsize=6, color='r')

    lines2, labs2 = ax2.get_legend_handles_labels()
    lines1, labs1 = ax1.get_legend_handles_labels()
    ax1.legend(lines1+lines2, labs1+labs2, fontsize=7)
    ax1.grid(True, alpha=0.3, axis='y')
    try:
        fig.tight_layout()
    except Exception:
        fig.subplots_adjust(bottom=0.15, top=0.88)

    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_amdahl_analysis.{ext}', dpi=150)
    plt.close(fig)
    print(f'  [OK] fig_amdahl_analysis')


# ── Figure 6: Integrator energy drift ────────────────────────
def plot_energy_drift(path='results/benchmarks/integrator_comparison.csv'):
    if not os.path.exists(path):
        print(f'  [SKIP] fig_energy_drift (run test_energy first)')
        return

    df  = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(3.5, 2.8))

    ax.semilogy(df['sim_year'], df['leapfrog_energy_drift'],
                'b-o', ms=4, lw=1.5, label='Leapfrog (symplectic)')
    ax.semilogy(df['sim_year'], df['rk4_energy_drift'],
                'r--s', ms=4, lw=1.5, label='RK4 (non-symplectic)')

    ax.set_xlabel('Simulation time [years]')
    ax.set_ylabel(r'Energy drift $|E-E_0|/|E_0|$')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    ratio_10yr = df['rk4_energy_drift'].iloc[-1] / max(df['leapfrog_energy_drift'].iloc[-1], 1e-20)
    ax.text(0.05, 0.95, f'RK4/Leapfrog at yr{df["sim_year"].iloc[-1]}: {ratio_10yr:.1f}x worse',
            transform=ax.transAxes, fontsize=7, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.tight_layout()
    for ext in ['pdf','png']:
        fig.savefig(f'{OUT}/fig_energy_drift.{ext}')
    plt.close(fig)
    print(f'  [OK] fig_energy_drift')


# ── Main ──────────────────────────────────────────────────────
def main():
    print()
    print('  Generating benchmark figures...')
    print()

    df       = load_sweep()
    amdahl_d = load_amdahl()
    mpi_df   = load_mpi()

    plot_speedup(df, amdahl_d)
    plot_efficiency(df)
    plot_absolute_time(df)
    plot_mpi_scaling(mpi_df)
    plot_amdahl_analysis(df)
    plot_energy_drift()

    print()
    print(f'  All benchmark figures saved to {OUT}/')
    print()

if __name__ == '__main__':
    main()
