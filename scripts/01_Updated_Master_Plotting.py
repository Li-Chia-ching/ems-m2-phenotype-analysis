#!/usr/bin/env python3
# =============================================================================
# EMS Mutagenesis — Updated Master Plotting Script (Python) — Publication-Ready
# Based on updated data with Root_StrongSystem REMOVED
#
# Figures generated:
#   Fig. 1: EMS Dose-Response (A: LL.4 survival | B: Quadratic mutation rate)
#   Fig. 3: Dual-Perspective Phenotype Frequency (A: Lollipop | B: Bubble matrix)
#   Fig. 4: Trait Correlation & Response (A: Pearson heatmap | B: Line plots)
#
# Aesthetics: seaborn paper context, custom academic palette, mathtext math
#
# Data Sources:
#   - Updated_EMS_M2_Cleaned_1425plants.csv
#   - Updated_Supplementary_Table_S1.csv
#   - Updated_Supplementary_Table_S2.csv
#   - EMS_M2_MP2603V2_Sheet1.csv (dose-response summary)
# =============================================================================

import os, sys, warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mticker
import seaborn as sns

warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════════════
# 0. GLOBAL CONFIGURATION — Publication-Ready Aesthetics
# ══════════════════════════════════════════════════════════════════════════════

PROJECT_ROOT = r"C:\Users\lijia\Documents\R Workplace\Lei_Thesis\Lei_EMS"
UPDATED_DATA_DIR = os.path.join(PROJECT_ROOT, "03_Outputs_and_Figures", "Updated_Data")
SHEET1_PATH = os.path.join(PROJECT_ROOT, "01_Raw_Data", "Archived", "EMS_M2_MP2603V2_Sheet1.csv")
OUT_DIR = os.path.join(PROJECT_ROOT, "03_Outputs_and_Figures")
os.makedirs(OUT_DIR, exist_ok=True)

# --- Seaborn theme: clean, paper-quality ---
sns.set_theme(style="ticks", context="paper", font_scale=1.15)

# --- Custom Academic Palette ---
# Tech Blue -> curves & main data; Accent Orange -> data points; Green -> annotations
CLR = {
    'blue':     '#2b7bba',
    'orange':   '#d95f02',
    'green':    '#1b9e77',
    'dark':     '#2d3436',
    'grey':     '#636e72',
    'light_bg': '#F8F9FA',
    'border':   '#DEE2E6',
    'faint':    '#dfe6e9',
}

# --- Matplotlib math rendering configuration ---
rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Arial', 'DejaVu Sans', 'Helvetica'],
    'mathtext.fontset':  'stix',
    'axes.unicode_minus': False,
    'figure.dpi':        150,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'savefig.pad_inches': 0.08,
    'pdf.fonttype':      42,
    'ps.fonttype':       42,
    'xtick.direction':   'in',
    'ytick.direction':   'in',
    'xtick.major.size':  3.5,
    'ytick.major.size':  3.5,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
})

N_TOTAL = 1425

# Shared annotation box style: elegant rounded box with light background
ANNO_BOX = dict(
    boxstyle='round,pad=0.5',
    facecolor=CLR['light_bg'],
    edgecolor=CLR['border'],
    alpha=0.92,
    linewidth=0.6,
)


# ══════════════════════════════════════════════════════════════════════════════
# 1. DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════
def load_data():
    """Load all input datasets required for plotting."""
    print("=" * 60)
    print("Loading data ...")

    m2 = pd.read_csv(os.path.join(UPDATED_DATA_DIR, "Updated_EMS_M2_Cleaned_1425plants.csv"))
    print(f"  M2_Cleaned: {m2.shape[0]} plants × {m2.shape[1]} cols")

    s1 = pd.read_csv(os.path.join(UPDATED_DATA_DIR, "Updated_Supplementary_Table_S1.csv"))
    print(f"  Table_S1:   {s1.shape[0]} families × {s1.shape[1]} cols")

    s2 = pd.read_csv(os.path.join(UPDATED_DATA_DIR, "Updated_Supplementary_Table_S2.csv"))
    print(f"  Table_S2:   {s2.shape[0]} rows × {s2.shape[1]} cols")

    dose = pd.read_csv(SHEET1_PATH)
    dose.columns = [c.strip() for c in dose.columns]
    print(f"  Sheet1:     {dose.shape[0]} concentrations")

    # Verify absence of 'Root_StrongSystem'
    pheno_cols = ['Phenotype_day120','Phenotype_day90','Phenotype_day60','Phenotype_day30']
    hits = sum(m2[c].astype(str).str.contains('Root', na=False).sum() for c in pheno_cols)
    assert hits == 0, f"ERROR: {hits} Root_StrongSystem entries found!"
    print(f"  ✓ 0 'Root_StrongSystem' entries confirmed\n")
    
    return {'m2': m2, 's1': s1, 's2': s2, 'dose': dose}


# ══════════════════════════════════════════════════════════════════════════════
# 2. PHENOTYPE EXTRACTION HELPER
# ══════════════════════════════════════════════════════════════════════════════
def extract_phenotype_categories(m2_df, col='Phenotype_day90_AllCategories'):
    """Split ';'-delimited AllCategories into (Family, category) record entries."""
    recs = []
    for _, row in m2_df.iterrows():
        s = str(row[col])
        if not s or s == 'nan':
            continue
        for cat in s.split(';'):
            cat = cat.strip()
            if cat and cat != 'nan':
                recs.append({'Family': row['Family'], 'category': cat})
    return pd.DataFrame(recs)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — EMS Dose-Response Analysis
# ══════════════════════════════════════════════════════════════════════════════
def generate_fig1(dose_df, out_dir):
    """
    Fig 1A: LL.4 Log-Logistic Dose-Response (M1 seedling survival)
    Fig 1B: Quadratic Mutation Rate Optimization (M2 mutation rate vs EMS conc.)
    
    ASSUMPTIONS:
    - Survival data from M1 dose-response (unaffected by Root_StrongSystem removal).
    - EMS >= 0.6% have no viable M2 offspring (mutation rate = nd).
    - Quadratic fitted to exactly 3 points (0.2, 0.3, 0.4%); df = 0, so no CI band.
    """
    print("=" * 60)
    print("Fig 1: EMS Dose-Response Analysis")

    # ── 1A. LL.4 model ───────────────────────────────────────────────────────
    surv = dose_df[['EMS_Concentration','Seedling Survival Rate_Mean']].copy()
    surv = surv[surv['EMS_Concentration'] > 0].dropna()
    surv.columns = ['conc','survival']
    surv['conc_pct'] = surv['conc'] * 100

    xd, yd = surv['conc'].values, surv['survival'].values

    # 4-parameter log-logistic function
    def ll4(x, c, d, e, b):
        return c + (d - c) / (1. + np.exp(b * (np.log(np.maximum(x,1e-10)) - np.log(e))))

    # Inverse LL.4 function for LD calculation
    def inv_ll4(p, c, d, e, b):
        target = d - p * (d - c)
        ratio  = (d - c) / (target - c) - 1
        if ratio <= 0:
            return np.nan
        return e * ratio ** (1. / b)

    try:
        # Fit LL.4 model
        popt, _ = curve_fit(ll4, xd, yd, p0=[0.075, 0.937, 0.006, 5.0],
                            bounds=([0.,0.85,0.002,1.],[0.2,1.,0.015,8.]),
                            maxfev=50000, ftol=1e-12)
        cV, dV, eV, bV = popt
        ld50V  = eV
        ld10V  = inv_ll4(0.10, cV, dV, eV, bV)
        ld90V  = inv_ll4(0.90, cV, dV, eV, bV)

        # Bootstrap CIs for LD50
        np.random.seed(2026)
        b50 = np.full(1000, np.nan)
        for i in range(1000):
            idx = np.random.choice(len(xd), len(xd), replace=True)
            try:
                pb, _ = curve_fit(ll4, xd[idx], yd[idx], p0=popt,
                                  bounds=([0.,0.85,0.002,1.],[0.2,1.,0.015,8.]),
                                  maxfev=5000)
                b50[i] = pb[2]
            except Exception:
                pass
        ci = (np.nanpercentile(b50,2.5), np.nanpercentile(b50,97.5))
        xp = np.linspace(xd.min(), xd.max(), 200)
        yp = ll4(xp, *popt)
        ll4_ok = True
        print(f"  LL.4: LD50={ld50V*100:.2f}%  (95% CI {ci[0]*100:.2f}–{ci[1]*100:.2f}%)")
        print(f"         LD10={ld10V*100:.2f}%   LD90={ld90V*100:.2f}%")
    except Exception as e:
        print(f"  WARNING — LL.4 fit failed: {e}")
        ll4_ok = False

    # ── 1B. Quadratic mutation-rate optimisation ──────────────────────────────
    mut = dose_df[['EMS_Concentration','M2 Mutation Rate_Mean']].dropna(subset=['M2 Mutation Rate_Mean'])
    mut.columns = ['conc','mut_rate']
    mut['conc_pct'] = mut['conc'] * 100
    mf = mut[mut['conc'] > 0]
    xm, ym = mf['conc'].values, mf['mut_rate'].values

    # Fit quadratic curve (degree 2)
    coeffs = np.polyfit(xm, ym, 2)
    aC, bC, cC = coeffs
    if aC < 0:
        optC = -bC / (2*aC)
    else:
        optC = np.nan

    # NOTE: With only 3 data points and a 3-parameter quadratic, the fit is
    # exact (df = 0). Neither residual-based nor bootstrap-based confidence
    # bands are statistically valid under these conditions — we intentionally
    # omit the 95 % CI ribbon from Fig 1B.
    xpm = np.linspace(xm.min(), xm.max(), 100)
    ypm = np.polyval(coeffs, xpm)

    print(f"  Quadratic: a={aC:.4f}, b={bC:.4f}, c={cC:.4f}")
    print(f"  C_opt = {optC*100:.2f}%  (df = 0, CI not estimable)")

    # ── Plotting ──────────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 5.0))

    # ---- Panel A -------------------------------------------------------------
    if ll4_ok:
        ax1.plot(xp*100, yp, color=CLR['blue'], lw=1.6, zorder=3)
        ax1.scatter(xd*100, yd, s=52, color=CLR['orange'],
                    edgecolors='white', linewidth=1.5, zorder=4, alpha=0.92)
        ax1.axhline(0.5, color=CLR['grey'], ls='--', lw=0.7, alpha=0.5, zorder=2)
        ax1.scatter([ld50V*100], [0.5], marker='D', s=64, color=CLR['green'],
                    edgecolors='white', linewidth=1.2, zorder=5)

        # LD50 floating annotation box (bottom-right of panel)
        l50, cLo, cHi = ld50V*100, ci[0]*100, ci[1]*100
        ld_str = (r'$\mathbf{LD_{50}}$ = ' + f'{l50:.2f}%'
                  + '\n95% CI: [{:.2f}%, {:.2f}%]'.format(cLo, cHi))
        ax1.text(0.97, 0.05, ld_str, transform=ax1.transAxes,
                 fontsize=7.5, ha='right', va='bottom', color=CLR['dark'],
                 bbox=ANNO_BOX, zorder=10)
    else:
        ax1.scatter(xd*100, yd, s=52, color=CLR['orange'],
                    edgecolors='white', linewidth=1.5, zorder=4)

    ax1.set_title('A   LL.4 Log-Logistic Dose-Response', fontsize=10.5,
                  fontweight='bold', loc='left', color=CLR['dark'])
    ax1.set_xlabel('EMS Concentration  (% v/v)', fontsize=9, color=CLR['dark'])
    ax1.set_ylabel(r'M$_1$ Seedling Survival Rate', fontsize=9, color=CLR['dark'])
    ax1.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax1.set_ylim(0, 1.06)
    sns.despine(ax=ax1)

    # ---- Panel B -------------------------------------------------------------
    ax2.plot(xpm*100, ypm, color=CLR['blue'], lw=1.6, zorder=3)
    ax2.scatter(xm*100, ym, s=58, color=CLR['orange'],
                edgecolors='white', linewidth=1.5, zorder=4, alpha=0.92)

    if not np.isnan(optC):
        ax2.axvline(optC*100, color=CLR['green'], ls='--', lw=0.9, zorder=2)

        # C_opt + regression formula in a floating box at BOTTOM-RIGHT
        info_lines = [
            r'$y = ' + f'{aC:.4f}' + r'x^2 + ' + f'{bC:.4f}' + r'x + ' + f'{cC:.4f}' + r'$',
            r'$\mathbf{C_{opt}}$ = ' + f'{optC*100:.2f}%',
        ]
        ax2.text(0.97, 0.05, '\n'.join(info_lines), transform=ax2.transAxes,
                 fontsize=7.5, ha='right', va='bottom', color=CLR['dark'],
                 bbox=ANNO_BOX, zorder=10)

    ax2.set_title('B   Quadratic Mutation Rate Optimization', fontsize=10.5,
                  fontweight='bold', loc='left', color=CLR['dark'])
    ax2.set_xlabel('EMS Concentration  (% v/v)', fontsize=9, color=CLR['dark'])
    ax2.set_ylabel(r'M$_2$ Mutation Rate', fontsize=9, color=CLR['dark'])
    ax2.yaxis.set_major_formatter(mticker.PercentFormatter(1.0, decimals=1))
    sns.despine(ax=ax2)

    # Shared: light y-grid only
    for ax in (ax1, ax2):
        ax.grid(True, axis='y', alpha=0.25, lw=0.4, color=CLR['faint'])
        ax.tick_params(labelsize=8, colors=CLR['dark'])

    fig.suptitle('EMS Dose-Response Analysis in Medicago polymorpha',
                 fontsize=12, fontweight='bold', y=1.01, color=CLR['dark'])
    fig.tight_layout()

    for fmt in ('pdf','svg','png'):
        fig.savefig(os.path.join(out_dir, f'Fig1_DoseResponse_Updated.{fmt}'),
                    facecolor='white')
    print(f"  -> Saved Fig1_DoseResponse_Updated.[pdf|svg|png]")
    plt.close(fig)

    # Export data for replotting if needed
    surv.copy().to_csv(os.path.join(out_dir, 'Fig1A_PlotData_Updated.csv'), index=False)
    mf.copy().to_csv(os.path.join(out_dir, 'Fig1B_PlotData_Updated.csv'), index=False)

    return {'ld50': ld50V*100 if ll4_ok else None, 'opt_conc': optC*100}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Dual-Perspective Phenotype Frequency
# ══════════════════════════════════════════════════════════════════════════════
def generate_fig3(m2_df, s1_df, out_dir):
    print("=" * 60)
    print("Fig 3: Dual-Perspective Phenotype Frequency")

    # Phenotype extraction (Day 90)
    pdf = extract_phenotype_categories(m2_df, 'Phenotype_day90_AllCategories')

    # ── Global frequencies GF_j ───────────────────────────────────────────────
    gf = pdf.groupby('category').size().reset_index(name='n')
    gf['GF']     = gf['n'] / N_TOTAL
    gf['GF_pct'] = gf['GF'] * 100
    gf = gf.sort_values('GF', ascending=True)

    # Dictionary to clean label names for plot readability
    CLEAN = {
        'Dwarfism_StuntedGrowth':'Dwarfism','Flowering_Early':'Early Flowering',
        'Multifoliate':'Multifoliate','Chlorosis_Pallor':'Chlorosis',
        'Branching_Clustered':'Clustered Branch','Leaf_Morphology_Small':'Small Leaves',
        'Stem_Elongated':'Elongated Stem','Stem_Thin':'Thin Stem',
        'Pigmentation_Anthocyanin':'Anthocyanin','Leaf_Morphology_Narrow':'Narrow Leaves',
        'Chimerism':'Chimerism','Leaf_Morphology_Folded':'Folded Leaves',
        'Senescence_Early':'Early Senescence','Lethal':'Lethal',
        'GrowthHabit_Prostrate':'Prostrate','GrowthHabit_Ascending':'Ascending',
        'Chlorosis_Albinism':'Albinism','GrowthHabit_Compact':'Compact Growth',
        'Branching_Reduced':'Reduced Branch','Disease_Susceptible':'Disease Susc.',
        'Pigmentation_GreenPetiole':'Green Petiole','Leaf_Morphology_DeepVein':'Deep Vein',
        'Leaf_Morphology_ShallowVein':'Shallow Vein','Leaf_Morphology_Round':'Round Leaves',
        'Leaf_Morphology_Large':'Large Leaves','Leaf_Morphology_Thin':'Thin Leaves',
        'Leaf_Morphology_NecroticMargin':'Necrotic Margin','Leaf_Morphology_Reduced':'Reduced Leaves',
        'Flowering_Sterile':'Flower Sterile','GrowthHabit_Erect':'Erect Growth',
        'Unclassified':'Unclassified',
    }
    gf['label'] = gf['category'].map(CLEAN).fillna(gf['category'])

    total = gf['n'].sum()
    print(f"  Categories: {len(gf)}  |  Total obs: {total}")
    for _, r in gf.nlargest(5,'GF').iterrows():
        print(f"    {r['label']:22s}  {r['n']:3d} / {N_TOTAL} = {r['GF_pct']:.2f}%")

    # ── Family frequencies FF_ij — Top 10 ─────────────────────────────────────
    fsz = m2_df.groupby('Family').size().reset_index(name='Nf')
    fbd = pdf.groupby('Family').size().reset_index(name='tot')
    fbd = fbd.merge(fsz, on='Family').sort_values('tot', ascending=False)
    
    # Extract the top 10 mutant-dense lineages
    top10 = fbd.head(10)['Family'].tolist()

    print(f"\n  Top 10 families:")
    for _, r in fbd.head(10).iterrows():
        print(f"    Fam {r['Family']:3d}: {r['tot']:2d} mutants / {r['Nf']:2d} plants")

    ff = pdf[pdf['Family'].isin(top10)]
    ff = ff.groupby(['Family','category']).size().reset_index(name='n')
    ff = ff.merge(fsz, on='Family')
    
    # Calculate family-specific frequency
    ff['FF']     = ff['n'] / ff['Nf']
    ff['FF_pct'] = ff['FF'] * 100
    ff['label']  = ff['category'].map(CLEAN).fillna(ff['category'])
    ff['Flab']   = ff['Family'].apply(
        lambda x: f'Fam {x}\n(n={fsz.loc[fsz.Family==x,"Nf"].values[0]})')

    valid = gf['category'].tolist()
    ff = ff[ff['category'].isin(valid)]
    ff['label'] = pd.Categorical(ff['label'], categories=gf['label'].tolist(), ordered=True)

    # Pivot table for bubble matrix
    piv = ff.pivot_table(index='label', columns='Flab', values='FF_pct',
                          aggfunc='sum', fill_value=0)
    fam_order = [
        f'Fam {f}\n(n={fsz.loc[fsz.Family==f,"Nf"].values[0]})' for f in top10
    ]
    fam_order = [f for f in fam_order if f in piv.columns]
    piv = piv.reindex(columns=fam_order)

    # ── Plotting ──────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14.5, 6.5))

    # Panel A — Lollipop chart
    axA = fig.add_axes([0.055, 0.10, 0.33, 0.84])
    yy = range(len(gf))
    axA.hlines(yy, 0, gf['GF_pct'], color=CLR['faint'], lw=1.0, zorder=2)
    axA.scatter(gf['GF_pct'], yy, s=55, color=CLR['blue'],
                edgecolors='white', linewidth=1.2, zorder=3, alpha=0.88)
    for i, (_, r) in enumerate(gf.iterrows()):
        axA.text(r['GF_pct'] + 0.35, i, f'{r["GF_pct"]:.2f}%',
                 fontsize=6.2, va='center', color=CLR['grey'])
    axA.set_yticks(yy)
    axA.set_yticklabels(gf['label'], fontsize=7, color=CLR['dark'])
    axA.set_xlabel('Global Frequency (%)', fontsize=9, color=CLR['dark'])
    axA.set_title('A   Global Phenotypic Mutation Spectrum\n'
                  r'$\mathbf{GF_j = n_j\ /\ 1425}$',
                  fontsize=10.5, fontweight='bold', loc='left', color=CLR['dark'])
    axA.set_xlim(0, gf['GF_pct'].max() * 1.22)
    axA.invert_yaxis()
    axA.grid(True, axis='x', alpha=0.25, lw=0.4, color=CLR['faint'])
    sns.despine(ax=axA)
    axA.tick_params(labelsize=7)

    # Panel B — Bubble matrix for top 10 families
    axB = fig.add_axes([0.43, 0.10, 0.55, 0.84])

    xl = piv.columns.tolist()
    yl = piv.index.tolist()

    bdata = []
    for ci, col in enumerate(xl):
        for ri, idx in enumerate(yl):
            v = piv.loc[idx, col] if (idx in piv.index and col in piv.columns) else 0
            if v > 0:
                bdata.append({'x':ci,'y':ri,'v':v})
    if bdata:
        bdf = pd.DataFrame(bdata)
        # Using rocket_r reversed gradient -> dark = high frequency
        sc = axB.scatter(bdf['x'], bdf['y'], s=bdf['v']*4.2,
                         c=bdf['v'], cmap='rocket_r', alpha=0.86,
                         edgecolors='white', linewidth=0.4, zorder=3)
        cbar = plt.colorbar(sc, ax=axB, shrink=0.75, aspect=22, pad=0.02)
        cbar.set_label('Family Freq. (%)', fontsize=8, color=CLR['dark'])
        cbar.ax.tick_params(labelsize=7, colors=CLR['dark'])
        cbar.outline.set_visible(False)

    axB.set_xticks(range(len(xl)))
    axB.set_xticklabels(xl, rotation=0, fontsize=6.8, color=CLR['dark'])
    axB.set_yticks(range(len(yl)))
    axB.set_yticklabels(yl, fontsize=6.8, color=CLR['dark'])
    axB.set_title('B   Family-Level Mutation Frequencies\n'
                  r'$\mathbf{FF_{ij} = n_{ij}\ /\ N_i}$  ·  Top 10 families',
                  fontsize=10.5, fontweight='bold', loc='left', color=CLR['dark'])
    axB.set_xlim(-0.5, len(xl)-0.5)
    axB.set_ylim(-0.5, len(yl)-0.5)
    axB.invert_yaxis()
    axB.grid(True, alpha=0.15, lw=0.3, color=CLR['faint'])
    sns.despine(ax=axB)
    axB.tick_params(labelsize=7)

    fig.suptitle('Phenotype Segregation and Trait Correlation Analysis\n'
                 f'M2 population:  N = {N_TOTAL} plants,  {m2_df["Family"].nunique()} families',
                 fontsize=12, fontweight='bold', y=1.01, color=CLR['dark'])

    for fmt in ('pdf','svg','png'):
        fig.savefig(os.path.join(out_dir, f'Fig3_DualFrequency_Updated.{fmt}'),
                    facecolor='white')
    print(f"  -> Saved Fig3_DualFrequency_Updated.[pdf|svg|png]")
    plt.close(fig)

    gf.to_csv(os.path.join(out_dir, 'Fig3A_GlobalFrequency_Updated.csv'), index=False)
    piv.to_csv(os.path.join(out_dir, 'Fig3B_FamilyFrequencyMatrix_Updated.csv'))
    return gf, piv


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Trait Correlation & Response
# ══════════════════════════════════════════════════════════════════════════════
def generate_fig4(dose_df, out_dir):
    print("=" * 60)
    print("Fig 4: Trait Correlation & Phenotypic Responses")

    mean_cols = [c for c in dose_df.columns if c.endswith('_Mean')]
    tdf = dose_df[['EMS_Concentration'] + mean_cols].copy()

    # Rename variables for aesthetic plotting
    RENAME = {
        'Germination Rate_Mean':'Germ. Rate',
        'Seedling Survival Rate_Mean':'Survival Rate',
        'Flowering Rate_Mean':'Flower Rate',
        'Pod Set Rate_Mean':'Pod Set',
        'Fertile Mutation Rate_Mean':'Fertile Mut.',
        'Harvested Pods / Plant_Mean':'Pods/Plant',
        'Seeds / Pod_Mean':'Seeds/Pod',
        'Harvested Seeds / Plant (Calc)_Mean':'Seeds/Plant',
        'Dead Embryo Ratio_Mean':'Dead Embryo',
        'M2 Germination Rate_Mean':'M2 Germ. Rate',
        'M2 Mutation Rate_Mean':'M2 Mut. Rate',
    }
    tdf = tdf.rename(columns=RENAME)
    tdf = tdf.rename(columns={'EMS_Concentration':'EMS Conc.'})

    # Keep columns with >= 3 non-NaN values
    keep = []
    for c in tdf.columns:
        if c == 'EMS Conc.' or tdf[c].notna().sum() >= 3:
            keep.append(c)
    cdat = tdf[keep]

    # Calculate Pearson correlation matrix
    cormat = cdat.corr(method='pearson')

    # ── Plotting ──────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14.5, 6.5))

    # --- Panel A: Correlation Heatmap ---
    axA = fig.add_axes([0.05, 0.10, 0.40, 0.84])
    # Mask the upper triangle
    mask = np.triu(np.ones_like(cormat, dtype=bool), k=1)

    sns.heatmap(cormat, mask=mask, annot=True, fmt='.2f',
                cmap='RdBu_r', center=0, vmin=-1, vmax=1, square=True,
                linewidths=0.6, linecolor='white',
                cbar_kws={'shrink':0.68, 'label':'Pearson r', 'pad':0.02},
                annot_kws={'fontsize':7, 'color':CLR['dark']},
                ax=axA)
    axA.set_title('A   Trait Correlation Matrix', fontsize=10.5,
                  fontweight='bold', loc='left', color=CLR['dark'])
    axA.tick_params(labelsize=7, colors=CLR['dark'])
    
    # Rotate y labels horizontally for readability
    axA.set_yticklabels(axA.get_yticklabels(), rotation=0, fontsize=7)
    axA.set_xticklabels(axA.get_xticklabels(), rotation=40, ha='right', fontsize=7)

    # --- Panel B: Multi-trait line plots (responses to EMS) ---
    recs = []
    for _, row in dose_df.iterrows():
        conc_pct = row['EMS_Concentration'] * 100
        for mc in mean_cols:
            mv = row[mc]
            sv = row.get(mc.replace('_Mean','_SD'), np.nan)
            if pd.notna(mv):
                recs.append({
                    'EMS_pct': conc_pct,
                    'trait': RENAME.get(mc, mc.replace('_Mean','')),
                    'mean': mv,
                    'sd': sv if pd.notna(sv) else 0,
                })
    pdf = pd.DataFrame(recs)
    utraits = pdf['trait'].unique()
    ncols = 3
    nrows = int(np.ceil(len(utraits) / ncols))

    gs = fig.add_gridspec(nrows, ncols, left=0.51, right=0.98, top=0.90,
                           bottom=0.10, hspace=0.60, wspace=0.50)

    for i, trait in enumerate(utraits):
        r, c = divmod(i, ncols)
        ax = fig.add_subplot(gs[r, c])
        tr = pdf[pdf['trait'] == trait].sort_values('EMS_pct')

        # Line plot with standard deviation error bars
        ax.plot(tr['EMS_pct'], tr['mean'], '-o', color=CLR['blue'],
                lw=1.0, ms=3.5, mfc=CLR['orange'], mec='white', mew=1.0, zorder=3)
        for _, trow in tr.iterrows():
            if trow['sd'] > 0:
                ax.errorbar(trow['EMS_pct'], trow['mean'], yerr=trow['sd'],
                            fmt='none', ecolor=CLR['grey'], capsize=2.5,
                            lw=0.6, alpha=0.55, zorder=2)
        ax.set_title(trait, fontsize=7.8, fontweight='bold', color=CLR['dark'])
        ax.tick_params(labelsize=6.2, colors=CLR['grey'])
        ax.grid(True, alpha=0.2, lw=0.3, color=CLR['faint'])
        sns.despine(ax=ax)
        if r == nrows - 1:
            ax.set_xlabel('EMS (% v/v)', fontsize=7, color=CLR['dark'])
        if c == 0:
            ax.set_ylabel('Value', fontsize=7, color=CLR['dark'])

    # Hide unused subplots
    for i in range(len(utraits), nrows*ncols):
        r, c = divmod(i, ncols)
        fig.add_subplot(gs[r, c]).set_visible(False)

    fig.suptitle('Trait Correlation and Phenotypic Trait Responses',
                 fontsize=12, fontweight='bold', y=1.01, color=CLR['dark'])
    fig.text(0.51, 0.94, 'B   Phenotypic Trait Responses to EMS Treatment',
             fontsize=10.5, fontweight='bold', ha='left', color=CLR['dark'])

    for fmt in ('pdf','svg','png'):
        fig.savefig(os.path.join(out_dir, f'Fig4_Correlation_Trends_Updated.{fmt}'),
                    facecolor='white')
    print(f"  -> Saved Fig4_Correlation_Trends_Updated.[pdf|svg|png]")
    plt.close(fig)

    cormat.to_csv(os.path.join(out_dir, 'Fig4A_CorrelationMatrix_Updated.csv'))
    pdf.to_csv(os.path.join(out_dir, 'Fig4B_TraitTrends_Updated.csv'), index=False)
    return cormat


# ══════════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 60)
    print("EMS M2 — Publication-Ready Master Plotting  (Root_StrongSystem excluded)")
    print("=" * 60)

    data = load_data()
    r1 = generate_fig1(data['dose'], OUT_DIR)
    r3 = generate_fig3(data['m2'], data['s1'], OUT_DIR)
    r4 = generate_fig4(data['dose'], OUT_DIR)

    print("\n" + "=" * 60)
    print("ALL FIGURES EXPORTED SUCCESSFULLY")
    print(f"  Export Directory: {OUT_DIR}")
    print(f"  LD50 = {r1['ld50']:.2f}%  |  C_opt = {r1['opt_conc']:.2f}%")
    print(f"  Categories = {len(r3[0])}  |  Total observations = {r3[0]['n'].sum()}")
    print("=" * 60)


if __name__ == '__main__':
    main()
