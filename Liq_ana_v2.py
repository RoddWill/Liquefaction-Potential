#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined Liquefaction Analysis Application

This application provides two separate analysis modes (via tabs):
  • SPT Liquefaction Analysis
  • CPT Liquefaction Analysis

Each mode uses its own Run Analysis function and processing routines.
Results (plots, text/CSV files, and an Excel workbook) are exported into a user‐selected output folder.

The top‐right block in each analysis window displays “How to Use” and “About” buttons,
which now are aligned to the right. The help and about messages reference the alternative
analysis mode so that users know that there are two separate windows (one for SPT and one for CPT).

Developed by Wroddy
Date: [Date]
Support: support@wroddy.com
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import warnings
import pandas as pd
from math import log, sqrt
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tempfile
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import brentq

# ---------------------------------------------------------
# Utility Functions (Shared by SPT & CPT)
# ---------------------------------------------------------
def create_output_folders(base_folder):
    """Create output folder with subfolders: models, plots, and results."""
    subfolders = ['models', 'plots', 'results']
    paths = {}
    os.makedirs(base_folder, exist_ok=True)
    for sf in subfolders:
        path = os.path.join(base_folder, sf)
        os.makedirs(path, exist_ok=True)
        paths[sf] = path
    return paths

def liquefaction_damage_function(FS):
    if FS >= 1.2:
        return 0.0
    elif FS >= 0.95:
        return 2.106 * np.exp(-18.427 * FS)
    else:
        return 1.0 - FS

def depth_weighting(z):
    if z > 20:
        return 0.0
    val = 10.0 - 0.5 * z
    return val if val > 0 else 0.0

def compute_local_LSN(xi_v, z_mid, thickness):
    """Compute local Liquefaction Severity Number (LSN) for a layer."""
    z_mid_adj = np.where(z_mid == 0, 1e-6, z_mid)
    return 1000 * (xi_v / z_mid_adj) * thickness

# ---------------------------------------------------------
# SPT Analysis Functions
# ---------------------------------------------------------
def compute_FS_spt(layer, water_table, Mw, PGA, rod_extension, borehole_diameter, ER, anal_depth,
                   Water_unit_weight=9.81, Pa=101.325):
    z_top = layer['top_depth']
    z_bot = min(layer['bottom_depth'], 20.0)
    thickness = z_bot - z_top
    if thickness <= 0:
        return np.nan
    z_mid = 0.5 * (z_top + z_bot)
    unit_weight = layer['unit_weight']
    NSPT_val = layer['NSPT']
    FC_val = layer['FC']
    
    total_stress = z_mid * unit_weight
    effective_stress = total_stress - max(0, (z_mid - water_table)) * Water_unit_weight
    rd = np.exp(-1.012 - 1.126 * np.sin(z_mid/11.73 + 5.133) + Mw*(0.106 + 0.118*np.sin(z_mid/11.28+5.142)))
    rd = min(rd, 1.0)
    CSR = 0.65 * PGA * (total_stress/effective_stress) * rd
    
    Ce = ER / 60.0
    if borehole_diameter < 115.1:
        Cb = 1.0
    elif borehole_diameter < 150.1:
        Cb = 1.05
    else:
        Cb = 1.15
    Cr_bins = np.array([0.75, 0.8, 0.85, 0.95, 1.0])
    depth_plus = z_mid + rod_extension
    idx = np.digitize(depth_plus, [3,4,6,10])
    Cr = Cr_bins[idx]
    Cs = 1.0
    
    if NSPT_val is None or np.isnan(NSPT_val) or NSPT_val <= 0:
        return np.nan
    N60 = NSPT_val * Ce * Cb * Cr * Cs
    delta_N160 = np.exp(1.63 + 9.7/(FC_val+0.01) - (15.7/(FC_val+0.01))**2)
    
    def f(n160cs):
        val1 = np.clip((Pa/effective_stress)**(0.784 - 0.0768*np.sqrt(max(n160cs,1))), None, 1.7)
        return val1 - ((n160cs - delta_N160)/N60)
    try:
        N160CS = brentq(f, 1, 200)
    except Exception:
        N160CS = np.nan
    N160_val = N160CS - delta_N160
    
    if z_mid > water_table and z_mid <= anal_depth and N160CS < 37.5:
        CRR_75 = np.exp(N160CS/14.1 + (N160CS/126)**2 - (N160CS/23.6)**3 + (N160CS/25.4)**4 - 2.8)
    else:
        CRR_75 = 2.0
    MSF = 1 + ((np.clip(1.09+(N160CS/31.5)**2, None, 2.2)-1)*(8.64*np.exp(-Mw/4)-1.325))
    K_sigma = 1 - (1/(18.9-2.55*np.sqrt(min(N160CS,37)))*np.log(effective_stress/Pa))
    K_sigma = np.clip(K_sigma, None, 1.1)
    K_alpha = 1.0
    CRR = np.clip(CRR_75 * K_sigma * K_alpha * MSF, None, 2)
    CRR = np.where(np.isnan(CRR), 2, CRR)
    if CSR == 0:
        FS_liq = np.nan
    else:
        FS_liq = CRR / CSR
    FS_liq = min(FS_liq, 2)
    return FS_liq

def compute_local_LPI(FS, z_mid, thickness):
    local = np.where(FS < 1, 1 - FS, 0) * np.where(z_mid <= 20, np.maximum(0, 10 - 0.5*z_mid), 0) * thickness
    return local

def compute_cumulative_LPI(local_LPI):
    return np.cumsum(local_LPI[::-1])[::-1]

def compute_LPI(layers):
    detailed_rows = []
    local_LPI_list = []
    for lyr in layers:
        z_top = lyr['top_depth']
        z_bot = lyr['bottom_depth']
        fos = lyr['FOS']
        if z_top >= 20:
            continue
        z_bot_clamped = min(z_bot, 20.0)
        thickness = z_bot_clamped - z_top
        if thickness <= 0:
            continue
        z_mid = 0.5 * (z_top + z_bot_clamped)
        Fz = liquefaction_damage_function(fos)
        Wz = depth_weighting(z_mid)
        local_LPI = Fz * Wz * thickness
        local_LPI_list.append(local_LPI)
        detailed_rows.append({
            'z_top': z_top,
            'z_bottom': z_bot_clamped,
            'z_mid': z_mid,
            'thickness': thickness,
            'FOS': fos,
            'Fz': Fz,
            'Wz': Wz,
            'local_LPI': local_LPI
        })
    local_LPI_array = np.array(local_LPI_list)
    cum_LPI = np.cumsum(local_LPI_array[::-1])[::-1]
    LPI_final = cum_LPI[0] if len(cum_LPI) > 0 else 0.0
    return LPI_final, detailed_rows, cum_LPI

def classify_LPI(lpi):
    if lpi < 2:
        return "Little or no liquefaction potential"
    elif lpi < 5:
        return "Low liquefaction potential"
    elif lpi < 15:
        return "Moderate liquefaction potential"
    else:
        return "High liquefaction potential"

def single_plot_SPT(NSPT, N160, FC, # FC_no_fill, 
                    CSR, CRR, FS_liq, LPI, LSN, depth, mid_depth,
                    water_table, ID, output_dir, export_pdf, export_png, gen_parameters):
    # 2x3 grid plot for SPT analysis
    fig, axes = plt.subplots(2, 3, sharey=True, figsize=(5, 6), facecolor="#f0f4f7")
    axes = axes.flatten()
    
    # Panel 1: NSPT & N160 vs. Depth
    for style, label in [('bo--', "NSPT"), ('gd-', "N160")]:
        axes[0].plot(NSPT, depth, style, label=label, linewidth=1.5, markersize=4)
    axes[0].set_xlabel("Blow Count", fontsize=8)
    axes[0].set_title("NSPT & N160 vs. Depth", fontsize=9)
    axes[0].legend(fontsize=7)
    axes[0].invert_yaxis()
    axes[0].xaxis.set_minor_locator(AutoMinorLocator())
    
    # Panel 2: Fines Content vs. Depth
    for style, label in [('rs--', "FC"), # (filled)"), ('m^-', "FC (raw)")
                         ]:
        axes[1].plot(FC, depth, style, label=label, linewidth=1.5, markersize=4)
    axes[1].set_xlabel("Fines Content (%)", fontsize=8)
    axes[1].set_title("Fines Content vs. Depth", fontsize=9)
    axes[1].legend(fontsize=7)
    axes[1].invert_yaxis()
    axes[1].xaxis.set_minor_locator(AutoMinorLocator())
    
    # Panel 3: CSR & CRR vs. Depth
    for style, label in [('cv-.', "CSR"), ('mp-.', "CRR")]:
        axes[2].plot(CSR, depth, style, label=label, linewidth=1.5, markersize=4)
    axes[2].set_xlabel("CSR & CRR", fontsize=8)
    axes[2].set_title("CSR & CRR vs. Depth", fontsize=9)
    axes[2].legend(fontsize=7)
    axes[2].invert_yaxis()
    axes[2].xaxis.set_minor_locator(AutoMinorLocator())
    
    # Panel 4: FS vs. Depth
    axes[3].plot(FS_liq, depth, 'k>-', label="FS_liq", linewidth=1.5, markersize=4)
    axes[3].set_xlabel("FS_liq", fontsize=8)
    axes[3].set_title("Factor-of-Safety vs. Depth", fontsize=9)
    axes[3].legend(fontsize=7)
    axes[3].invert_yaxis()
    axes[3].xaxis.set_minor_locator(AutoMinorLocator())
    
    # Panel 5: LPI vs. Depth with secondary LSN axis
    axes[4].plot(LPI, mid_depth, 'ko--', label="LPI", linewidth=1.5, markersize=4)
    axes[4].set_xlabel("LPI", fontsize=8)
    axes[4].set_title("LPI vs. Depth", fontsize=9)
    axes[4].invert_yaxis()
    axes[4].xaxis.set_minor_locator(AutoMinorLocator())
    ax_sec = axes[4].twiny()
    ax_sec.plot(LSN, mid_depth, 'rx:', label="LSN", linewidth=1.5, markersize=4)
    ax_sec.set_xlabel("LSN", color='red', fontsize=8)
    ax_sec.tick_params(axis='x', colors='red')
    ax_sec.xaxis.set_minor_locator(AutoMinorLocator())
    lines1, labels1 = axes[4].get_legend_handles_labels()
    lines2, labels2 = ax_sec.get_legend_handles_labels()
    axes[4].legend(lines1 + lines2, labels1 + labels2, fontsize=7)
    
    # Panel 6: Dummy plot for NSPT vs. Corrected N160CS
    N160CS = NSPT * 1.1  # Dummy corrected value
    axes[5].plot(NSPT, depth, 'bo--', label="NSPT", linewidth=1.5, markersize=4)
    axes[5].plot(N160CS, depth, 'mp-.', label="N160CS", linewidth=1.5, markersize=4)
    axes[5].set_xlabel("Blow Count", fontsize=8)
    axes[5].set_title("NSPT and N160CS vs. Depth", fontsize=9)
    axes[5].legend(fontsize=7)
    axes[5].invert_yaxis()
    axes[5].xaxis.set_minor_locator(AutoMinorLocator())
    
    axes[0].set_ylabel("Depth (m)", fontsize=8)
    axes[0].yaxis.set_minor_locator(AutoMinorLocator())
    if depth.size > 0:
        plt.ylim(bottom=0, top=np.ceil(np.max(depth)))
    else:
        plt.ylim(bottom=0, top=10)
    
    fig.suptitle(f"Comprehensive SPT Liquefaction Analysis for Borehole {ID}\n{gen_parameters}", fontsize=8)
    fig.text(0.1, 0.01, "Analysis by SPT LPI App", fontsize=6, color='#34495e')
    fig.text(0.98, 0.01, "Liq_ana v1.0 | Developed by Wroddy", fontsize=6, color='#7f8c8d', ha='right')
    plt.subplots_adjust(left=0.08, bottom=0.10, right=0.98, top=0.88, wspace=0.3, hspace=0.5)
    
    if export_pdf in ['Borings+SUMMARY', 'Borings']:
        pdf_path = os.path.join(output_dir, ID + '_combined.pdf')
        plt.savefig(pdf_path, format='pdf', dpi=900, bbox_inches='tight')
    if export_png in ['Borings+SUMMARY', 'Borings']:
        png_path = os.path.join(output_dir, ID + '_combined.png')
        plt.savefig(png_path, format='png', bbox_inches='tight')
    
    return fig

def process_sheet(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                  borehole_diameter, ER, anal_depth, xi_v_value):
    try:
        df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
    except Exception as e:
        print(f"Error selecting required columns: {e}")
        return None
    layers = []
    for idx, row in df.iterrows():
        try:
            layer = {
                'top_depth': float(row['top_depth']),
                'bottom_depth': float(row['bottom_depth']),
                'NSPT': float(row['NSPT']),
                'FC': float(row['FC']),
                'unit_weight': float(row['unit_weight'])
            }
            layers.append(layer)
        except Exception as e:
            print(f"Row {idx} skipped: {e}")
    for lyr in layers:
        lyr['FOS'] = compute_FS_spt(lyr, water_table, Mw, PGA, rod_extension, borehole_diameter, ER, anal_depth)
    LPI_final, layer_details, LPI_profile = compute_LPI(layers)
    classification = classify_LPI(LPI_final)
    print(f"[{sheet_name}] Final cumulative LPI = {LPI_final:.2f} --> {classification}")
    
    base_name = sheet_name.replace(" ", "_")
    results_txt = os.path.join(out_paths['results'], f'{base_name}_LPI_results.txt')
    with open(results_txt, 'w', encoding="utf-8") as f:
        f.write("Liquefaction Potential Index (LPI) Results\n")
        f.write("------------------------------------------\n")
        f.write(f"Final cumulative LPI: {LPI_final:.2f}\n")
        f.write(f"Classification: {classification}\n")
    details_csv = os.path.join(out_paths['results'], f'{base_name}_layer_details.csv')
    with open(details_csv, 'w', newline='', encoding="utf-8") as cf:
        writer = csv.writer(cf)
        writer.writerow(['z_top','z_bottom','z_mid','thickness','FOS','Fz','Wz','local_LPI'])
        for row in layer_details:
            writer.writerow([row['z_top'], row['z_bottom'], row['z_mid'], row['thickness'],
                             row['FOS'], row['Fz'], row['Wz'], row['local_LPI']])
    
    depth_arr = np.array([0.5*(lyr['top_depth']+min(lyr['bottom_depth'],20)) for lyr in layers if lyr['top_depth'] < 20])
    NSPT_arr = np.array([lyr['NSPT'] for lyr in layers if lyr['top_depth'] < 20])
    N160_arr = NSPT_arr * 0.9
    FC_arr = np.array([lyr['FC'] for lyr in layers if 'FC' in lyr])
    FC_no_fill_arr = FC_arr * 1.05
    CSR_arr = np.array([1.0]*len(depth_arr))  # Dummy values
    CRR_arr = CSR_arr * 1.1                  # Dummy values
    FS_arr = np.array([lyr['FOS'] for lyr in layers if lyr['top_depth'] < 20])
    
    mid_depth_arr = np.array([row['z_mid'] for row in layer_details])
    local_LPI_arr = np.array([row['local_LPI'] for row in layer_details])
    cum_LPI_arr = LPI_profile
    local_LSN_arr = np.array([compute_local_LSN(xi_v_value, row['z_mid'], row['thickness']) for row in layer_details])
    cum_LSN_arr = np.cumsum(local_LSN_arr[::-1])[::-1]
    
    gen_params = f"Parameters: Mw={Mw}, PGA={PGA} g, Water Table={water_table} m"
    ID = f"{base_name}_Borehole"
    
    fig = single_plot_SPT(NSPT_arr, N160_arr, FC_arr, # FC_no_fill_arr, 
                          CSR_arr, CRR_arr, FS_arr,
                          cum_LPI_arr, cum_LSN_arr, depth_arr, mid_depth_arr, water_table, ID,
                          out_paths['plots'], export_pdf="Borings+SUMMARY", export_png="Borings+SUMMARY",
                          gen_parameters=gen_params)
    
    return {
        "sheet": sheet_name,
        "summary": {"LPI_final": LPI_final, "classification": classification},
        "layer_details": pd.DataFrame(layer_details),
        "figure": fig
    }

# ---------------------------------------------------------
# CPT Analysis Functions
# ---------------------------------------------------------
def compute_FS_cpt(layer, water_table, Mw, PGA, rod_extension, borehole_diameter, ER, anal_depth,
                   Water_unit_weight=9.81, Pa=101.325):
    """
    Compute the factor-of-safety for a CPT layer using accurate CPT equations.
    The measured qc and fs are weakened (qc*0.01, fs*0.02) and converted to kPa.
    """
    z_top = layer['top_depth']
    z_bot = min(layer['bottom_depth'], anal_depth)
    thickness = z_bot - z_top
    if thickness <= 0:
        return np.nan
    z_mid = 0.5 * (z_top + z_bot)
    unit_weight = layer['unit_weight']
    total_stress = z_mid * unit_weight
    effective_stress = total_stress - max(0, (z_mid - water_table)) * Water_unit_weight
    if effective_stress <= 0:
        return np.nan
    # Weaken measured values
    qc_meas = layer['qc'] * 0.01
    fs_meas = layer['fs'] * 0.02
    qc_kPa = qc_meas * 1000.0
    fs_kPa = fs_meas * 1000.0
    Q_term = qc_kPa - (total_stress/Pa)
    if Q_term <= 0:
        Q_term = 1e-6
    n_temp = 1.0
    Q = Q_term * (Pa/effective_stress)**n_temp
    denom = qc_kPa - total_stress
    if denom <= 0:
        denom = 1e-6
    F = (fs_kPa/denom) * 100.0
    try:
        Ic = sqrt((3.47 - log(Q))**2 + (log(F) + 1.22)**2) / 3.0
    except Exception:
        Ic = 0.0
    n = 0.381 * Ic + 0.05 * (effective_stress/Pa) - 0.15
    n = min(n, 1.0)
    CQ = (Pa/effective_stress)**n
    qc1N = CQ * (qc_kPa/Pa)
    if qc1N < 50:
        CRR7_5 = 0.833 * (qc1N/1000.0) + 0.05
    elif qc1N < 160:
        CRR7_5 = 93 * (qc1N/1000.0)**3 + 0.08
    else:
        CRR7_5 = 2.0
    rd = np.exp(-1.012 - 1.126 * np.sin(z_mid/11.73 + 5.133) +
                Mw * (0.106 + 0.118 * np.sin(z_mid/11.28 + 5.142)))
    rd = min(rd, 1.0)
    CSR = 0.65 * PGA * (total_stress/effective_stress) * rd
    if CSR == 0:
        return np.nan
    f_exp = 0.7
    K_sigma = (effective_stress/Pa)**(f_exp-1)
    MSF = 1.0
    K_alpha = 1.0
    FS = (CRR7_5/CSR) * MSF * K_sigma * K_alpha
    FS = min(FS, 2.0)
    return FS

def compute_LPI_CPT(layers):
    """
    Compute cumulative LPI for CPT analysis.
    Local LPI = [1 - FS] * [10 - 0.5*z_mid] * (layer thickness) for layers with top_depth < anal_depth.
    """
    detailed_rows = []
    local_LPI_list = []
    for lyr in layers:
        if "top_depth" in lyr and "bottom_depth" in lyr:
            z_top = lyr["top_depth"]
            z_bot = lyr["bottom_depth"]
        else:
            continue
        if z_top >= lyr.get("anal_depth", 20):
            continue
        z_bot_clamped = min(z_bot, lyr.get("anal_depth", 20))
        thickness = z_bot_clamped - z_top
        if thickness <= 0:
            continue
        z_mid = 0.5 * (z_top + z_bot_clamped)
        Fz = liquefaction_damage_function(lyr["FOS"])
        Wz = depth_weighting(z_mid)
        local_LPI = Fz * Wz * thickness
        detailed_rows.append({
            "top_depth": z_top,
            "bottom_depth": z_bot_clamped,
            "z_mid": z_mid,
            "thickness": thickness,
            "FOS": lyr["FOS"],
            "Fz": Fz,
            "Wz": Wz,
            "local_LPI": local_LPI
        })
        local_LPI_list.append(local_LPI)
    local_LPI_array = np.array(local_LPI_list)
    cum_LPI = np.cumsum(local_LPI_array[::-1])[::-1]
    LPI_final = cum_LPI[0] if len(cum_LPI) > 0 else 0.0
    return LPI_final, detailed_rows, cum_LPI

def classify_LSN(lsn):
    if lsn < 5:
        return "Very low severity"
    elif lsn < 10:
        return "Low severity"
    elif lsn < 20:
        return "Moderate severity"
    elif lsn < 30:
        return "High severity"
    else:
        return "Very high severity"

def single_plot_CPT(qc, fs, CSR, CRR, FS_vals, LPI, LSN, depth, mid_depth,
                    water_table, ID, output_dir, export_pdf, export_png, gen_parameters):
    # 2x3 grid plot for CPT analysis
    fig, axes = plt.subplots(2, 3, sharey=True, figsize=(7, 6), facecolor="#f0f4f7")
    axes = axes.flatten()
    
    axes[0].plot(qc, depth, 'bo--', label="qc (MPa)", linewidth=1.5, markersize=4)
    axes[0].set_xlabel("qc (MPa)", fontsize=8)
    axes[0].set_title("Cone Tip Resistance", fontsize=9)
    axes[0].legend(fontsize=7)
    axes[0].invert_yaxis()
    axes[0].xaxis.set_minor_locator(AutoMinorLocator())
    
    axes[1].plot(fs, depth, 'rs--', label="fs (MPa)", linewidth=1.5, markersize=4)
    axes[1].set_xlabel("fs (MPa)", fontsize=8)
    axes[1].set_title("Sleeve Friction", fontsize=9)
    axes[1].legend(fontsize=7)
    axes[1].invert_yaxis()
    axes[1].xaxis.set_minor_locator(AutoMinorLocator())
    
    Rf = np.where(qc != 0, (fs / qc)*100, np.nan)
    axes[2].plot(Rf, depth, 'ms-', label="Friction Ratio (%)", linewidth=1.5, markersize=4)
    axes[2].set_xlabel("Friction Ratio (%)", fontsize=8)
    axes[2].set_title("Friction Ratio vs. Depth", fontsize=9)
    axes[2].legend(fontsize=7)
    axes[2].invert_yaxis()
    axes[2].xaxis.set_minor_locator(AutoMinorLocator())
    
    for style, label in [('cv-.', "CSR"), ('mp-.', "CRR")]:
        axes[3].plot(CSR, depth, style, label=label, linewidth=1.5, markersize=4)
    axes[3].set_xlabel("CSR & CRR", fontsize=8)
    axes[3].set_title("CSR & CRR (sample)", fontsize=9)
    axes[3].legend(fontsize=7)
    axes[3].invert_yaxis()
    axes[3].xaxis.set_minor_locator(AutoMinorLocator())
    
    axes[4].plot(FS_vals, depth, 'k>-', label="FS_liq", linewidth=1.5, markersize=4)
    axes[4].set_xlabel("FS_liq", fontsize=8)
    axes[4].set_title("Factor-of-Safety", fontsize=9)
    axes[4].legend(fontsize=7)
    axes[4].invert_yaxis()
    axes[4].xaxis.set_minor_locator(AutoMinorLocator())
    
    axes[5].plot(LPI, mid_depth, 'ko--', label="LPI", linewidth=1.5, markersize=4)
    axes[5].set_xlabel("LPI", fontsize=8)
    axes[5].set_title("LPI vs. Depth", fontsize=9)
    axes[5].invert_yaxis()
    axes[5].xaxis.set_minor_locator(AutoMinorLocator())
    ax_sec = axes[5].twiny()
    ax_sec.plot(LSN, mid_depth, 'rx:', label="LSN", linewidth=1.5, markersize=4)
    ax_sec.set_xlabel("LSN", color='red', fontsize=8)
    ax_sec.tick_params(axis='x', colors='red')
    ax_sec.xaxis.set_minor_locator(AutoMinorLocator())
    lines1, labels1 = axes[5].get_legend_handles_labels()
    lines2, labels2 = ax_sec.get_legend_handles_labels()
    axes[5].legend(lines1 + lines2, labels1 + labels2, fontsize=7)
    
    axes[0].set_ylabel("Depth (m)", fontsize=8)
    axes[0].yaxis.set_minor_locator(AutoMinorLocator())
    if depth.size > 0:
        plt.ylim(bottom=0, top=np.ceil(np.max(depth)))
    else:
        plt.ylim(bottom=0, top=10)
    
    fig.suptitle(f"CPT Liquefaction Analysis - Borehole {ID}\n{gen_parameters}\nWater Table = {water_table} m", fontsize=14)
    fig.text(0.1, 0.01, "Analysis by CPT Liquefaction App", fontsize=6, color='#34495e')
    fig.text(0.98, 0.01, "Liq_ana v1.0 | Developed by Wroddy", fontsize=6, color='#7f8c8d', ha='right')
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.88, wspace=0.3)
    
    if export_pdf in ['Borings+SUMMARY', 'Borings']:
        pdf_path = os.path.join(output_dir, ID + '_combined.pdf')
        plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    if export_png in ['Borings+SUMMARY', 'Borings']:
        png_path = os.path.join(output_dir, ID + '_combined.png')
        plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    
    return fig

def process_sheet_CPT(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                      borehole_diameter, ER, anal_depth, xi_v_value, default_unit_weight):
    if "qc" not in df.columns or "fs" not in df.columns:
        messagebox.showerror("Error", f"Sheet {sheet_name} is missing required CPT columns ('qc' and 'fs').")
        return None
    if "unit_weight" not in df.columns:
        df["unit_weight"] = default_unit_weight
    if "top_depth" in df.columns and "bottom_depth" in df.columns:
        df["top_depth"] = pd.to_numeric(df["top_depth"], errors="coerce")
        df["bottom_depth"] = pd.to_numeric(df["bottom_depth"], errors="coerce")
    elif "depth" in df.columns:
        df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
        df = df.sort_values("depth").reset_index(drop=True)
        df["top_depth"] = df["depth"]
        df["bottom_depth"] = df["depth"].shift(-1)
        df = df.dropna(subset=["bottom_depth"])
    else:
        messagebox.showerror("Error", f"Sheet {sheet_name} must have either 'depth' or both 'top_depth' and 'bottom_depth'.")
        return None
    layers = []
    for idx, row in df.iterrows():
        try:
            layer = {
                'top_depth': float(row['top_depth']),
                'bottom_depth': float(row['bottom_depth']),
                'qc': float(row['qc']),
                'fs': float(row['fs']),
                'unit_weight': float(row['unit_weight'])
            }
            layers.append(layer)
        except Exception as e:
            print(f"CPT Row {idx} skipped: {e}")
    for lyr in layers:
        lyr['FOS'] = compute_FS_cpt(lyr, water_table, Mw, PGA, rod_extension, borehole_diameter, ER, anal_depth)
        lyr["anal_depth"] = anal_depth  # add for clamping in LPI computation
    LPI_final, layer_details, LPI_profile = compute_LPI_CPT(layers)
    lpi_class = classify_LPI(LPI_final)
    print(f"[CPT - {sheet_name}] Final cumulative LPI = {LPI_final:.2f} --> {lpi_class}")
    base_name = sheet_name.replace(" ", "_")
    results_txt = os.path.join(out_paths['results'], f'{base_name}_CPT_LPI_results.txt')
    with open(results_txt, 'w', encoding="utf-8") as f:
        f.write("CPT Liquefaction Potential Index (LPI) Results\n")
        f.write("----------------------------------------------\n")
        f.write(f"Final cumulative LPI: {LPI_final:.2f}\n")
        f.write(f"Classification: {lpi_class}\n")
    details_csv = os.path.join(out_paths['results'], f'{base_name}_CPT_layer_details.csv')
    with open(details_csv, 'w', newline='', encoding="utf-8") as cf:
        writer = csv.writer(cf)
        writer.writerow(['top_depth','bottom_depth','z_mid','thickness','FOS','Fz','Wz','local_LPI'])
        for row in layer_details:
            writer.writerow([row['top_depth'], row['bottom_depth'], row['z_mid'], row['thickness'],
                             row['FOS'], row['Fz'], row['Wz'], row['local_LPI']])
    
    depth_arr = np.array([0.5*(lyr['top_depth']+min(lyr['bottom_depth'], anal_depth)) for lyr in layers if lyr['top_depth'] < anal_depth])
    qc_arr = np.array([lyr['qc'] for lyr in layers if lyr['top_depth'] < anal_depth])
    fs_arr = np.array([lyr['fs'] for lyr in layers if lyr['top_depth'] < anal_depth])
    CSR_arr = np.array([1.0]*len(depth_arr))  # Dummy values
    CRR_arr = CSR_arr * 1.1                  # Dummy values
    mid_depth_arr = np.array([row['z_mid'] for row in layer_details])
    cum_LPI_arr = LPI_profile
    local_LSN_arr = np.array([compute_local_LSN(xi_v_value, row['z_mid'], row['thickness']) for row in layer_details])
    cum_LSN_arr = np.cumsum(local_LSN_arr[::-1])[::-1]
    gen_params = f"CPT: Mw={Mw}, PGA={PGA} g, Water Table={water_table} m"
    ID = f"{base_name}_CPT_Borehole"
    fig = single_plot_CPT(qc_arr, fs_arr, CSR_arr, CRR_arr,
                          np.array([lyr['FOS'] for lyr in layers if lyr['top_depth'] < anal_depth]),
                          cum_LPI_arr, cum_LSN_arr, depth_arr, mid_depth_arr,
                          water_table, ID, out_paths['plots'],
                          export_pdf="Borings+SUMMARY", export_png="Borings+SUMMARY",
                          gen_parameters=gen_params)
    return {
        "sheet": sheet_name,
        "analysis": "CPT",
        "summary": {"LPI_final": LPI_final, "classification": lpi_class},
        "layer_details": pd.DataFrame(layer_details),
        "figure": fig
    }

# ---------------------------------------------------------
# Export Results to Excel (Shared)
# ---------------------------------------------------------
def export_results_to_excel(results, output_excel_path):
    writer = pd.ExcelWriter(output_excel_path, engine='xlsxwriter')
    workbook  = writer.book
    summary_data = []
    for res in results:
        summary_data.append({
            "Sheet": res["sheet"],
            "Analysis": res.get("analysis", "SPT"),
            "Final LPI": res["summary"]["LPI_final"],
            "Classification": res["summary"]["classification"]
        })
    df_summary = pd.DataFrame(summary_data)
    df_summary.to_excel(writer, sheet_name="Summary", index=False)
    for res in results:
        sheet_name = res["sheet"][:31]
        df_details = res["layer_details"]
        df_details.to_excel(writer, sheet_name=sheet_name, index=False)
    temp_img = os.path.join(tempfile.gettempdir(), "combined_plot.png")
    if results and "figure" in results[0]:
        results[0]["figure"].savefig(temp_img, format="png", bbox_inches="tight")
        plot_sheet = workbook.add_worksheet("Plot")
        plot_sheet.insert_image("A1", temp_img, {"x_scale": 0.8, "y_scale": 0.8})
    writer.close()

# ---------------------------------------------------------
# GUI Classes
# ---------------------------------------------------------
class SPTFrame(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.configure(bg="#2c3e50")
        self.input_file = None
        self.sheet_list = []
        self.output_folder = None
        self.results = []
        self.figure_canvas = None
        self.create_widgets()
    
    def create_widgets(self):
        # Left control panel remains unchanged
        self.columnconfigure(0, weight=1,minsize=300)
        self.columnconfigure(1, weight=3, minsize=750)
        self.rowconfigure(0, weight=1)
        
        self.left_frame = tk.Frame(self, bg="#2c3e50")
        self.left_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        
        top_frame = tk.Frame(self.left_frame, bg="#2c3e50")
        top_frame.pack(fill="x", pady=5)
        tk.Label(top_frame, text="SPT Liquefaction Analysis", font=("Helvetica", 16, "bold"), bg="#2c3e50", fg="white").pack(side="left")
      
        file_frame = tk.LabelFrame(self.left_frame, text="File Selection", bg="#34495e", fg="white", padx=5, pady=5)
        file_frame.pack(fill="x", pady=5)
        tk.Button(file_frame, text="Select Input File", command=self.select_file, bg="#2980b9", fg="white", font=("Helvetica", 9)).grid(row=0, column=0, padx=5, pady=5)
        self.file_label = tk.Label(file_frame, text="No file selected", fg="white", bg="#34495e", font=("Helvetica", 9))
        self.file_label.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(file_frame, text="Select Output Folder", command=self.select_output_folder, bg="#2980b9", fg="white", font=("Helvetica", 9)).grid(row=1, column=0, padx=5, pady=5)
        self.out_folder_label = tk.Label(file_frame, text="Output Folder", fg="lightgreen", bg="#34495e", font=("Helvetica", 9))
        self.out_folder_label.grid(row=1, column=1, padx=5, pady=5)
        
        self.sheet_frame = tk.LabelFrame(self.left_frame, text="Sheet Selection (Excel)", bg="#34495e", fg="white", padx=5, pady=5)
        self.sheet_frame.pack(fill="x", pady=5)
        self.sheet_listbox = tk.Listbox(self.sheet_frame, height=4, bg="#ecf0f1", font=("Helvetica", 9))
        self.sheet_listbox.pack(fill="x", padx=5, pady=5)
        
        param_frame = tk.LabelFrame(self.left_frame, text="Analysis Parameters", bg="#34495e", fg="white", padx=5, pady=5)
        param_frame.pack(fill="x", pady=5)
        params = [
            ("Water Table (m):", 5.0),
            ("Mw:", 7.0),
            ("PGA (g):", 0.3),
            ("Rod Extension (m):", 1.0),
            ("Borehole Diameter (mm):", 100),
            ("Energy Ratio (%):", 20),
            ("Analysis Depth (m):", 20.0),
            ("xi_v (for LSN):", 0.05)
        ]
        self.entries = {}
        for label_text, default in params:
            frame = tk.Frame(param_frame, bg="#34495e")
            frame.pack(fill="x", pady=2, padx=5)
            tk.Label(frame, text=label_text, width=25, anchor="w", bg="#34495e", fg="white", font=("Helvetica", 9)).pack(side="left")
            entry = tk.Entry(frame, font=("Helvetica", 9))
            entry.pack(side="left", fill="x", expand=True)
            entry.insert(0, str(default))
            self.entries[label_text] = entry
        
        tk.Button(self.left_frame, text="Run SPT Analysis", command=self.run_analysis, bg="#2980b9", fg="white", font=("Helvetica", 12, "bold")).pack(pady=10)
        tk.Label(self.left_frame, text="Developed by Wroddy", bg="#2c3e50", fg="lightgray", font=("Helvetica", 9)).pack(side="bottom", pady=5)
        
        # Right panel: Plot display with additional top-right info block
        self.right_frame = tk.Frame(self, bg="#f0f4f7")
        self.right_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        self.right_frame.columnconfigure(0, weight=1)
        self.right_frame.rowconfigure(1, weight=1)
        
        self.info_frame = tk.Frame(self.right_frame, bg="#e0e0e0", bd=2, relief="groove")
        self.info_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        # Aligning the buttons to the right:
        about_btn = tk.Button(self.info_frame, text="About", command=self.show_about_right, bg="#2980b9", fg="white", font=("Helvetica", 9))
        help_btn = tk.Button(self.info_frame, text="How to Use", command=self.show_help_right, bg="#2980b9", fg="white", font=("Helvetica", 9))
        about_btn.pack(side="right", padx=5, pady=5)
        help_btn.pack(side="right", padx=5, pady=5)
        
        self.canvas = tk.Canvas(self.right_frame, bg="#f0f4f7")
        self.canvas.grid(row=1, column=0, sticky="nsew")
        self.h_scroll = tk.Scrollbar(self.right_frame, orient="horizontal", command=self.canvas.xview)
        self.v_scroll = tk.Scrollbar(self.right_frame, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=self.h_scroll.set, yscrollcommand=self.v_scroll.set)
        self.h_scroll.grid(row=2, column=0, sticky="ew")
        self.v_scroll.grid(row=1, column=1, sticky="ns")
        self.fig_frame = tk.Frame(self.canvas, bg="#f0f4f7")
        self.canvas.create_window((0, 0), window=self.fig_frame, anchor="nw")
    
    
    def show_help_right(self):
        help_text = (
            "How to Use SPT Liquefaction Analysis Application\n\n"
            "Step-by-Step Instructions:\n"
            "1. Select Input File:\n"
            "   - Choose a CSV or Excel file containing SPT data with these columns:\n"
            "     top_depth, bottom_depth, NSPT, FC, unit_weight.\n\n"
            "2. Select Output Folder:\n"
            "   - This is where all result files (plots, CSV, text, and Excel workbook) will be saved.\n\n"
            "3. Set Analysis Parameters:\n"
            "   - Water Table (m): Depth of the groundwater table.\n"
            "   - Mw: Seismic moment magnitude parameter.\n"
            "   - PGA (g): Peak Ground Acceleration in g's.\n"
            "   - Rod Extension (m) & Borehole Diameter (mm): Equipment parameters affecting SPT values.\n"
            "   - Energy Ratio (%): Used in SPT corrections.\n"
            "   - Analysis Depth (m): Maximum depth for analysis (max 20 m for LPI calculation).\n"
            "   - xi_v (for LSN): Parameter for computing the Liquefaction Severity Number (LSN).\n\n"
            "Computed Values:\n"
            "   FS_liq: Factor-of-Safety against liquefaction. Values below 1 indicate liquefaction risk.\n"
            "   LPI: Liquefaction Potential Index, a cumulative indicator of liquefaction risk. Its classification:\n"
            "        <2  : Little or no liquefaction potential\n"
            "         2-5: Low liquefaction potential\n"
            "         5-15: Moderate liquefaction potential\n"
            "         >15: High liquefaction potential\n"
            "   LSN: Liquefaction Severity Number, indicating the severity of liquefaction effects.\n\n"
            "4. Click 'Run SPT Analysis' to process the data and view the combined plot.\n\n"
            "Note: For instructions on CPT analysis, please refer to the CPT Analysis tab."
        )
        help_win = tk.Toplevel(self)
        help_win.title("How to Use SPT Analysis")
        help_win.geometry("600x500")
        help_win.configure(bg="#ecf0f1")
        text = tk.Text(help_win, wrap="word", bg="#ecf0f1", fg="#2c3e50", font=("Helvetica", 10))
        text.insert("1.0", help_text)
        text.config(state="disabled")
        text.pack(fill="both", expand=True, padx=10, pady=10)
        tk.Button(help_win, text="Close", command=help_win.destroy, bg="#2980b9", fg="white").pack(pady=5)
    
    def show_about_right(self):
        about_text = (
            "SPT Liquefaction Analysis Application\n\n"
            "Version: 1.0\n"
            "Developed by Wroddy\n\n"
            "Description:\n"
            "This tool computes key liquefaction parameters using SPT data:\n"
            "   - FS_liq (Factor-of-Safety against liquefaction): Values below 1 indicate potential liquefaction.\n"
            "   - LPI (Liquefaction Potential Index): A cumulative measure indicating liquefaction risk,\n"
            "        classified as:\n"
            "           <2  : Little or no liquefaction potential\n"
            "           2-5 : Low liquefaction potential\n"
            "           5-15: Moderate liquefaction potential\n"
            "           >15 : High liquefaction potential\n"
            "   - LSN (Liquefaction Severity Number): Indicates the potential severity of liquefaction impacts.\n\n"
            "Usage:\n"
            "   1. Load your SPT data file and select an output folder.\n"
            "   2. Set the analysis parameters as required.\n"
            "   3. Run the analysis to view the plot and save the results.\n\n"
            "For CPT analysis instructions, please switch to the CPT Analysis tab.\n\n"
            "For support, please contact: support@wroddy.com"
        )
        about_win = tk.Toplevel(self)
        about_win.title("About SPT Analysis")
        about_win.geometry("600x450")
        about_win.configure(bg="#ecf0f1")
        text = tk.Text(about_win, wrap="word", bg="#ecf0f1", fg="#2c3e50", font=("Helvetica", 10))
        text.insert("1.0", about_text)
        text.config(state="disabled")
        text.pack(fill="both", expand=True, padx=10, pady=10)
        tk.Button(about_win, text="Close", command=about_win.destroy, bg="#2980b9", fg="white").pack(pady=5)
    
    def select_file(self):
        filetypes = [("CSV Files", "*.csv"), ("Excel Files", "*.xlsx *.xls")]
        filename = filedialog.askopenfilename(title="Select Input File", filetypes=filetypes)
        if filename:
            if os.path.getsize(filename) == 0:
                messagebox.showerror("Error", "Selected file is empty!")
                return
            self.input_file = filename
            self.file_label.config(text=os.path.basename(filename))
            ext = os.path.splitext(filename)[1].lower()
            if ext in [".xls", ".xlsx"]:
                try:
                    sheets = pd.ExcelFile(filename).sheet_names
                    self.sheet_list = sheets
                    self.sheet_listbox.delete(0, tk.END)
                    for sheet in sheets:
                        self.sheet_listbox.insert(tk.END, sheet)
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to read Excel sheets: {e}")
            else:
                self.sheet_listbox.delete(0, tk.END)
                self.sheet_list = []
    
    def select_output_folder(self):
        folder = filedialog.askdirectory(title="Select Output Folder")
        if folder:
            self.output_folder = folder
            self.out_folder_label.config(text=os.path.basename(folder), fg="lightgreen")
    
    def run_analysis(self):
        if not self.input_file:
            messagebox.showerror("Error", "Please select an input file.")
            return
        try:
            water_table = float(self.entries["Water Table (m):"].get())
            Mw = float(self.entries["Mw:"].get())
            PGA = float(self.entries["PGA (g):"].get())
            rod_extension = float(self.entries["Rod Extension (m):"].get())
            borehole_diameter = float(self.entries["Borehole Diameter (mm):"].get())
            ER = float(self.entries["Energy Ratio (%):"].get())
            anal_depth = float(self.entries["Analysis Depth (m):"].get())
            xi_v_value = float(self.entries["xi_v (for LSN):"].get())
        except Exception as e:
            messagebox.showerror("Error", f"Invalid parameter value: {e}")
            return
        
        if not self.output_folder:
            self.select_output_folder()
            if not self.output_folder:
                messagebox.showerror("Error", "No output folder selected.")
                return
        
        out_paths = create_output_folders(self.output_folder)
        ext = os.path.splitext(self.input_file)[1].lower()
        self.results = []
        fig_to_display = None
        
        try:
            if ext == ".csv":
                df = pd.read_csv(self.input_file)
                df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
                res = process_sheet(df, "CSV_Sheet", out_paths, water_table, Mw, PGA, rod_extension,
                                      borehole_diameter, ER, anal_depth, xi_v_value)
                if res is not None:
                    self.results.append(res)
                    fig_to_display = res["figure"]
            elif ext in [".xls", ".xlsx"]:
                sheets = pd.read_excel(self.input_file, sheet_name=None)
                selected = self.sheet_listbox.curselection()
                if selected:
                    for i in selected:
                        sheet_name = self.sheet_list[i]
                        df = sheets[sheet_name]
                        try:
                            df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
                        except Exception as e:
                            messagebox.showerror("Error", f"Sheet {sheet_name} missing required columns: {e}")
                            continue
                        res = process_sheet(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                                              borehole_diameter, ER, anal_depth, xi_v_value)
                        if res is not None:
                            self.results.append(res)
                            fig_to_display = res["figure"]
                else:
                    for sheet_name, df in sheets.items():
                        try:
                            df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
                        except Exception as e:
                            messagebox.showerror("Error", f"Sheet {sheet_name} missing required columns: {e}")
                            continue
                        res = process_sheet(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                                              borehole_diameter, ER, anal_depth, xi_v_value)
                        if res is not None:
                            self.results.append(res)
                            fig_to_display = res["figure"]
                        break
            else:
                messagebox.showerror("Error", "Unsupported file type. Please provide CSV or XLSX.")
                return
        except Exception as e:
            messagebox.showerror("Error", f"Error processing file: {e}")
            return
        
        if fig_to_display is not None:
            for widget in self.fig_frame.winfo_children():
                widget.destroy()
            self.figure_canvas = FigureCanvasTkAgg(fig_to_display, master=self.fig_frame)
            self.figure_canvas.draw()
            self.figure_canvas.get_tk_widget().pack(fill="both", expand=True)
            self.canvas.config(scrollregion=self.canvas.bbox("all"))
        
        if self.results:
            export_path = os.path.join(self.output_folder, "Combined_Analysis_Results.xlsx")
            export_results_to_excel(self.results, export_path)
            print(f"Exported combined results to {export_path}")
        
        messagebox.showinfo("Success", "Analysis completed. Check the output folder for results.")

class CPTFrame(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.configure(bg="#2c3e50")
        self.input_file = None
        self.sheet_list = []
        self.output_folder = None
        self.results = []
        self.figure_canvas = None
        self.create_widgets()
    
    def create_widgets(self):
        self.columnconfigure(0, weight=1, minsize=300)
        self.columnconfigure(1, weight=3, minsize=750)
        self.rowconfigure(0, weight=1)
        
        self.left_frame = tk.Frame(self, bg="#2c3e50")
        self.left_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        
        top_frame = tk.Frame(self.left_frame, bg="#2c3e50")
        top_frame.pack(fill="x", pady=5)
        tk.Label(top_frame, text="CPT Liquefaction Analysis", font=("Helvetica", 16, "bold"), bg="#2c3e50", fg="white").pack(side="left")
        
        file_frame = tk.LabelFrame(self.left_frame, text="File Selection", bg="#34495e", fg="white", padx=5, pady=5)
        file_frame.pack(fill="x", pady=5)
        tk.Button(file_frame, text="Select Input File", command=self.select_file, bg="#2980b9", fg="white", font=("Helvetica", 9)).grid(row=0, column=0, padx=5, pady=5)
        self.file_label = tk.Label(file_frame, text="No file selected", fg="white", bg="#34495e", font=("Helvetica", 9))
        self.file_label.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(file_frame, text="Select Output Folder", command=self.select_output_folder, bg="#2980b9", fg="white", font=("Helvetica", 9)).grid(row=1, column=0, padx=5, pady=5)
        self.out_folder_label = tk.Label(file_frame, text="Output Folder", fg="lightgreen", bg="#34495e", font=("Helvetica", 9))
        self.out_folder_label.grid(row=1, column=1, padx=5, pady=5)
        
        self.sheet_frame = tk.LabelFrame(self.left_frame, text="Sheet Selection (Excel)", bg="#34495e", fg="white", padx=5, pady=5)
        self.sheet_frame.pack(fill="x", pady=5)
        self.sheet_listbox = tk.Listbox(self.sheet_frame, height=4, bg="#ecf0f1", font=("Helvetica", 9))
        self.sheet_listbox.pack(fill="x", padx=5, pady=5)
        
        param_frame = tk.LabelFrame(self.left_frame, text="Analysis Parameters", bg="#34495e", fg="white", padx=5, pady=5)
        param_frame.pack(fill="x", pady=5)
        params = [
            ("Water Table (m):", 1.0),
            ("Mw:", 7.0),
            ("PGA (g):", 0.3),
            ("Rod Extension (m):", 1.0),
            ("Borehole Diameter (mm):", 100),
            ("Energy Ratio (%):", 20),
            ("Analysis Depth (m):", 20.0),
            ("xi_v (for LSN):", 0.05),
            ("Default Unit Weight (kN/m³):", 18.0)
        ]
        self.entries = {}
        for label_text, default in params:
            frame = tk.Frame(param_frame, bg="#34495e")
            frame.pack(fill="x", pady=2, padx=5)
            tk.Label(frame, text=label_text, width=30, anchor="w", bg="#34495e", fg="white", font=("Helvetica", 9)).pack(side="left")
            entry = tk.Entry(frame, font=("Helvetica", 9))
            entry.pack(side="left", fill="x", expand=True)
            entry.insert(0, str(default))
            self.entries[label_text] = entry
        
        tk.Button(self.left_frame, text="Run CPT Analysis", command=self.run_analysis, bg="#2980b9", fg="white", font=("Helvetica", 12, "bold")).pack(pady=10)
        tk.Label(self.left_frame, text="Developed by Wroddy", bg="#2c3e50", fg="lightgray", font=("Helvetica", 9)).pack(side="bottom", pady=5)
        
        # Right panel: Plot display with additional top-right info block
        self.right_frame = tk.Frame(self, bg="#f0f4f7")
        self.right_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        self.right_frame.columnconfigure(0, weight=1)
        self.right_frame.rowconfigure(1, weight=1)
        
        self.info_frame = tk.Frame(self.right_frame, bg="#e0e0e0", bd=2, relief="groove")
        self.info_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        # Aligning the buttons to the right:
        about_btn = tk.Button(self.info_frame, text="About", command=self.show_about_right, bg="#2980b9", fg="white", font=("Helvetica", 9))
        help_btn = tk.Button(self.info_frame, text="How to Use", command=self.show_help_right, bg="#2980b9", fg="white", font=("Helvetica", 9))
        about_btn.pack(side="right", padx=5, pady=5)
        help_btn.pack(side="right", padx=5, pady=5)
        
        self.canvas = tk.Canvas(self.right_frame, bg="#f0f4f7")
        self.canvas.grid(row=1, column=0, sticky="nsew")
        self.h_scroll = tk.Scrollbar(self.right_frame, orient="horizontal", command=self.canvas.xview)
        self.v_scroll = tk.Scrollbar(self.right_frame, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=self.h_scroll.set, yscrollcommand=self.v_scroll.set)
        self.h_scroll.grid(row=2, column=0, sticky="ew")
        self.v_scroll.grid(row=1, column=1, sticky="ns")
        self.fig_frame = tk.Frame(self.canvas, bg="#f0f4f7")
        self.canvas.create_window((0, 0), window=self.fig_frame, anchor="nw")
    
    def show_help_right(self):
        help_text = (
            "How to Use CPT Liquefaction Analysis Application\n\n"
            "Step-by-Step Instructions:\n"
            "1. Select Input File:\n"
            "   - Choose a CSV or Excel file containing CPT data with these columns:\n"
            "     top_depth, bottom_depth, qc, fs. If 'unit_weight' is missing, a default value will be used.\n\n"
            "2. Select Output Folder:\n"
            "   - This folder will store the result files (plots, CSV, text, and Excel workbook).\n\n"
            "3. Set Analysis Parameters:\n"
            "   - Water Table (m): Depth of the groundwater table.\n"
            "   - Mw: Seismic parameter.\n"
            "   - PGA (g): Peak Ground Acceleration.\n"
            "   - Rod Extension & Borehole Diameter: Equipment parameters for the CPT test.\n"
            "   - Energy Ratio (%): For correcting CPT values.\n"
            "   - Analysis Depth (m): Maximum depth for analysis (typically up to 20 m).\n"
            "   - xi_v (for LSN): Parameter for computing the Liquefaction Severity Number (LSN).\n"
            "   - Default Unit Weight (kN/m³): Used if not provided in the data.\n\n"
            "Computed Values:\n"
            "   FS_liq: Factor-of-Safety against liquefaction; values below 1 indicate liquefaction risk.\n"
            "   LPI: Liquefaction Potential Index, which aggregates the liquefaction potential over depth.\n"
            "        Classification:\n"
            "           <2  : Little or no liquefaction potential\n"
            "           2-5 : Low liquefaction potential\n"
            "           5-15: Moderate liquefaction potential\n"
            "           >15 : High liquefaction potential\n"
            "   LSN: Liquefaction Severity Number, indicating the potential severity of liquefaction effects.\n\n"
            "4. Click 'Run CPT Analysis' to process the data and view the combined plot.\n\n"
            "Note: For instructions on SPT analysis, please switch to the SPT Analysis tab."
        )
        help_win = tk.Toplevel(self)
        help_win.title("How to Use CPT Analysis")
        help_win.geometry("650x550")
        help_win.configure(bg="#ecf0f1")
        text = tk.Text(help_win, wrap="word", bg="#ecf0f1", fg="#2c3e50", font=("Helvetica", 10))
        text.insert("1.0", help_text)
        text.config(state="disabled")
        text.pack(fill="both", expand=True, padx=10, pady=10)
        tk.Button(help_win, text="Close", command=help_win.destroy, bg="#2980b9", fg="white").pack(pady=5)
    
    def show_about_right(self):
        about_text = (
            "CPT Liquefaction Analysis Application\n\n"
            "Version: 1.0\n"
            "Developed by Wroddy\n\n"
            "Description:\n"
            "This application processes CPT data to compute key liquefaction parameters:\n"
            "   - FS_liq (Factor-of-Safety): A safety measure; values below 1 suggest liquefaction risk.\n"
            "   - LPI (Liquefaction Potential Index): A cumulative index indicating liquefaction potential.\n"
            "        Classification:\n"
            "           <2  : Little or no liquefaction potential\n"
            "           2-5 : Low liquefaction potential\n"
            "           5-15: Moderate liquefaction potential\n"
            "           >15 : High liquefaction potential\n"
            "   - LSN (Liquefaction Severity Number): Indicates the potential severity of liquefaction impacts.\n\n"
            "Usage:\n"
            "   1. Load your CPT data file and select an output folder.\n"
            "   2. Set the required analysis parameters.\n"
            "   3. Run the analysis to view the generated plots and save detailed results.\n\n"
            "For SPT analysis instructions, please switch to the SPT Analysis tab.\n\n"
            "For support, please contact: support@wroddy.com"
        )
        about_win = tk.Toplevel(self)
        about_win.title("About CPT Analysis")
        about_win.geometry("650x500")
        about_win.configure(bg="#ecf0f1")
        text = tk.Text(about_win, wrap="word", bg="#ecf0f1", fg="#2c3e50", font=("Helvetica", 10))
        text.insert("1.0", about_text)
        text.config(state="disabled")
        text.pack(fill="both", expand=True, padx=10, pady=10)
        tk.Button(about_win, text="Close", command=about_win.destroy, bg="#2980b9", fg="white").pack(pady=5)
    
    def select_file(self):
        filetypes = [("CSV Files", "*.csv"), ("Excel Files", "*.xlsx *.xls")]
        filename = filedialog.askopenfilename(title="Select Input File", filetypes=filetypes)
        if filename:
            if os.path.getsize(filename) == 0:
                messagebox.showerror("Error", "Selected file is empty!")
                return
            self.input_file = filename
            self.file_label.config(text=os.path.basename(filename))
            ext = os.path.splitext(filename)[1].lower()
            if ext in [".xls", ".xlsx"]:
                try:
                    sheets = pd.ExcelFile(filename).sheet_names
                    self.sheet_list = sheets
                    self.sheet_listbox.delete(0, tk.END)
                    for sheet in sheets:
                        self.sheet_listbox.insert(tk.END, sheet)
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to read Excel sheets: {e}")
            else:
                self.sheet_listbox.delete(0, tk.END)
                self.sheet_list = []
    
    def select_output_folder(self):
        folder = filedialog.askdirectory(title="Select Output Folder")
        if folder:
            self.output_folder = folder
            self.out_folder_label.config(text=os.path.basename(folder), fg="lightgreen")
    
    def run_analysis(self):
        if not self.input_file:
            messagebox.showerror("Error", "Please select an input file.")
            return
        try:
            water_table = float(self.entries["Water Table (m):"].get())
            Mw = float(self.entries["Mw:"].get())
            PGA = float(self.entries["PGA (g):"].get())
            rod_extension = float(self.entries["Rod Extension (m):"].get())
            borehole_diameter = float(self.entries["Borehole Diameter (mm):"].get())
            ER = float(self.entries["Energy Ratio (%):"].get())
            anal_depth = float(self.entries["Analysis Depth (m):"].get())
            xi_v_value = float(self.entries["xi_v (for LSN):"].get())
        except Exception as e:
            messagebox.showerror("Error", f"Invalid parameter value: {e}")
            return
        
        if not self.output_folder:
            self.select_output_folder()
            if not self.output_folder:
                messagebox.showerror("Error", "No output folder selected.")
                return
        
        out_paths = create_output_folders(self.output_folder)
        ext = os.path.splitext(self.input_file)[1].lower()
        self.results = []
        fig_to_display = None
        
        try:
            if ext == ".csv":
                df = pd.read_csv(self.input_file)
                res = process_sheet(df, "CSV_Sheet", out_paths, water_table, Mw, PGA, rod_extension,
                                      borehole_diameter, ER, anal_depth, xi_v_value)
                if res is not None:
                    self.results.append(res)
                    fig_to_display = res["figure"]
            elif ext in [".xls", ".xlsx"]:
                sheets = pd.read_excel(self.input_file, sheet_name=None)
                selected = self.sheet_listbox.curselection()
                if selected:
                    for i in selected:
                        sheet_name = self.sheet_list[i]
                        df = sheets[sheet_name]
                        try:
                            df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
                        except Exception as e:
                            messagebox.showerror("Error", f"Sheet {sheet_name} missing required columns: {e}")
                            continue
                        res = process_sheet(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                                              borehole_diameter, ER, anal_depth, xi_v_value)
                        if res is not None:
                            self.results.append(res)
                            fig_to_display = res["figure"]
                else:
                    for sheet_name, df in sheets.items():
                        try:
                            df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
                        except Exception as e:
                            messagebox.showerror("Error", f"Sheet {sheet_name} missing required columns: {e}")
                            continue
                        res = process_sheet(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                                              borehole_diameter, ER, anal_depth, xi_v_value)
                        if res is not None:
                            self.results.append(res)
                            fig_to_display = res["figure"]
                        break
            else:
                messagebox.showerror("Error", "Unsupported file type. Please provide CSV or XLSX.")
                return
        except Exception as e:
            messagebox.showerror("Error", f"Error processing file: {e}")
            return
        
        if fig_to_display is not None:
            for widget in self.fig_frame.winfo_children():
                widget.destroy()
            self.figure_canvas = FigureCanvasTkAgg(fig_to_display, master=self.fig_frame)
            self.figure_canvas.draw()
            self.figure_canvas.get_tk_widget().pack(fill="both", expand=True)
            self.canvas.config(scrollregion=self.canvas.bbox("all"))
        
        if self.results:
            export_path = os.path.join(self.output_folder, "Combined_Analysis_Results.xlsx")
            export_results_to_excel(self.results, export_path)
            print(f"Exported combined results to {export_path}")
        
        messagebox.showinfo("Success", "Analysis completed. Check the output folder for results.")

class CPTFrame(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.configure(bg="#2c3e50")
        self.input_file = None
        self.sheet_list = []
        self.output_folder = None
        self.results = []
        self.figure_canvas = None
        self.create_widgets()
    
    def create_widgets(self):
        self.columnconfigure(0, weight=1, minsize=300)
        self.columnconfigure(1, weight=3, minsize=750)
        self.rowconfigure(0, weight=1)
        
        self.left_frame = tk.Frame(self, bg="#2c3e50")
        self.left_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        
        top_frame = tk.Frame(self.left_frame, bg="#2c3e50")
        top_frame.pack(fill="x", pady=5)
        tk.Label(top_frame, text="CPT Liquefaction Analysis", font=("Helvetica", 16, "bold"), bg="#2c3e50", fg="white").pack(side="left")
        
        file_frame = tk.LabelFrame(self.left_frame, text="File Selection", bg="#34495e", fg="white", padx=5, pady=5)
        file_frame.pack(fill="x", pady=5)
        tk.Button(file_frame, text="Select Input File", command=self.select_file, bg="#2980b9", fg="white", font=("Helvetica", 9)).grid(row=0, column=0, padx=5, pady=5)
        self.file_label = tk.Label(file_frame, text="No file selected", fg="white", bg="#34495e", font=("Helvetica", 9))
        self.file_label.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(file_frame, text="Select Output Folder", command=self.select_output_folder, bg="#2980b9", fg="white", font=("Helvetica", 9)).grid(row=1, column=0, padx=5, pady=5)
        self.out_folder_label = tk.Label(file_frame, text="Output Folder", fg="lightgreen", bg="#34495e", font=("Helvetica", 9))
        self.out_folder_label.grid(row=1, column=1, padx=5, pady=5)
        
        self.sheet_frame = tk.LabelFrame(self.left_frame, text="Sheet Selection (Excel)", bg="#34495e", fg="white", padx=5, pady=5)
        self.sheet_frame.pack(fill="x", pady=5)
        self.sheet_listbox = tk.Listbox(self.sheet_frame, height=4, bg="#ecf0f1", font=("Helvetica", 9))
        self.sheet_listbox.pack(fill="x", padx=5, pady=5)
        
        param_frame = tk.LabelFrame(self.left_frame, text="Analysis Parameters", bg="#34495e", fg="white", padx=5, pady=5)
        param_frame.pack(fill="x", pady=5)
        params = [
            ("Water Table (m):", 1.0),
            ("Mw:", 7.0),
            ("PGA (g):", 0.3),
            ("Rod Extension (m):", 1.0),
            ("Borehole Diameter (mm):", 100),
            ("Energy Ratio (%):", 20),
            ("Analysis Depth (m):", 20.0),
            ("xi_v (for LSN):", 0.05),
            ("Default Unit Weight (kN/m³):", 18.0)
        ]
        self.entries = {}
        for label_text, default in params:
            frame = tk.Frame(param_frame, bg="#34495e")
            frame.pack(fill="x", pady=2, padx=5)
            tk.Label(frame, text=label_text, width=30, anchor="w", bg="#34495e", fg="white", font=("Helvetica", 9)).pack(side="left")
            entry = tk.Entry(frame, font=("Helvetica", 9))
            entry.pack(side="left", fill="x", expand=True)
            entry.insert(0, str(default))
            self.entries[label_text] = entry
        
        tk.Button(self.left_frame, text="Run CPT Analysis", command=self.run_analysis, bg="#2980b9", fg="white", font=("Helvetica", 12, "bold")).pack(pady=10)
        tk.Label(self.left_frame, text="Developed by Wroddy", bg="#2c3e50", fg="lightgray", font=("Helvetica", 9)).pack(side="bottom", pady=5)
        
        # Right panel: Plot display with additional top-right info block
        self.right_frame = tk.Frame(self, bg="#f0f4f7")
        self.right_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        self.right_frame.columnconfigure(0, weight=1)
        self.right_frame.rowconfigure(1, weight=1)
        
        self.info_frame = tk.Frame(self.right_frame, bg="#e0e0e0", bd=2, relief="groove")
        self.info_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        # Align the buttons to the right
        about_btn = tk.Button(self.info_frame, text="About", command=self.show_about_right, bg="#2980b9", fg="white", font=("Helvetica", 9))
        help_btn = tk.Button(self.info_frame, text="How to Use", command=self.show_help_right, bg="#2980b9", fg="white", font=("Helvetica", 9))
        about_btn.pack(side="right", padx=5, pady=5)
        help_btn.pack(side="right", padx=5, pady=5)
        
        self.canvas = tk.Canvas(self.right_frame, bg="#f0f4f7")
        self.canvas.grid(row=1, column=0, sticky="nsew")
        self.h_scroll = tk.Scrollbar(self.right_frame, orient="horizontal", command=self.canvas.xview)
        self.v_scroll = tk.Scrollbar(self.right_frame, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=self.h_scroll.set, yscrollcommand=self.v_scroll.set)
        self.h_scroll.grid(row=2, column=0, sticky="ew")
        self.v_scroll.grid(row=1, column=1, sticky="ns")
        self.fig_frame = tk.Frame(self.canvas, bg="#f0f4f7")
        self.canvas.create_window((0, 0), window=self.fig_frame, anchor="nw")
    
    def show_help_right(self):
        help_text = (
            "How to Use CPT Liquefaction Analysis Application\n\n"
            "Step-by-Step Instructions:\n"
            "1. Select Input File:\n"
            "   - Choose a CSV or Excel file containing CPT data with these columns:\n"
            "     top_depth, bottom_depth, qc, fs. If 'unit_weight' is missing, a default value will be used.\n\n"
            "2. Select Output Folder:\n"
            "   - This folder will store the result files (plots, CSV, text, and Excel workbook).\n\n"
            "3. Set Analysis Parameters:\n"
            "   - Water Table (m): Depth of the groundwater table.\n"
            "   - Mw: Seismic parameter.\n"
            "   - PGA (g): Peak Ground Acceleration.\n"
            "   - Rod Extension & Borehole Diameter: Equipment parameters for the CPT test.\n"
            "   - Energy Ratio (%): For correcting CPT values.\n"
            "   - Analysis Depth (m): Maximum depth for analysis (typically up to 20 m).\n"
            "   - xi_v (for LSN): Parameter for computing the Liquefaction Severity Number (LSN).\n"
            "   - Default Unit Weight (kN/m³): Used if not provided in the data.\n\n"
            "Computed Values:\n"
            "   FS_liq: Factor-of-Safety against liquefaction; values below 1 indicate liquefaction risk.\n"
            "   LPI: Liquefaction Potential Index, which aggregates the liquefaction potential over depth.\n"
            "        Classification:\n"
            "           <2  : Little or no liquefaction potential\n"
            "           2-5 : Low liquefaction potential\n"
            "           5-15: Moderate liquefaction potential\n"
            "           >15 : High liquefaction potential\n"
            "   LSN: Liquefaction Severity Number, indicating the potential severity of liquefaction effects.\n\n"
            "4. Click 'Run CPT Analysis' to process the data and view the combined plot.\n\n"
            "Note: For instructions on SPT analysis, please switch to the SPT Analysis tab."
        )
        help_win = tk.Toplevel(self)
        help_win.title("How to Use CPT Analysis")
        help_win.geometry("650x550")
        help_win.configure(bg="#ecf0f1")
        text = tk.Text(help_win, wrap="word", bg="#ecf0f1", fg="#2c3e50", font=("Helvetica", 10))
        text.insert("1.0", help_text)
        text.config(state="disabled")
        text.pack(fill="both", expand=True, padx=10, pady=10)
        tk.Button(help_win, text="Close", command=help_win.destroy, bg="#2980b9", fg="white").pack(pady=5)
    
    def show_about_right(self):
        about_text = (
            "CPT Liquefaction Analysis Application\n\n"
            "Version: 1.0\n"
            "Developed by Wroddy\n\n"
            "Description:\n"
            "This application processes CPT data to compute key liquefaction parameters:\n"
            "   - FS_liq (Factor-of-Safety): A safety measure; values below 1 suggest liquefaction risk.\n"
            "   - LPI (Liquefaction Potential Index): A cumulative index indicating liquefaction potential.\n"
            "        Classification:\n"
            "           <2  : Little or no liquefaction potential\n"
            "           2-5 : Low liquefaction potential\n"
            "           5-15: Moderate liquefaction potential\n"
            "           >15 : High liquefaction potential\n"
            "   - LSN (Liquefaction Severity Number): Indicates the potential severity of liquefaction impacts.\n\n"
            "Usage:\n"
            "   1. Load your CPT data file and select an output folder.\n"
            "   2. Set the required analysis parameters.\n"
            "   3. Run the analysis to view the generated plots and save detailed results.\n\n"
            "For SPT analysis instructions, please switch to the SPT Analysis tab.\n\n"
            "For support, please contact: support@wroddy.com"
        )
        about_win = tk.Toplevel(self)
        about_win.title("About CPT Analysis")
        about_win.geometry("650x500")
        about_win.configure(bg="#ecf0f1")
        text = tk.Text(about_win, wrap="word", bg="#ecf0f1", fg="#2c3e50", font=("Helvetica", 10))
        text.insert("1.0", about_text)
        text.config(state="disabled")
        text.pack(fill="both", expand=True, padx=10, pady=10)
        tk.Button(about_win, text="Close", command=about_win.destroy, bg="#2980b9", fg="white").pack(pady=5)
    
    def select_file(self):
        filetypes = [("CSV Files", "*.csv"), ("Excel Files", "*.xlsx *.xls")]
        filename = filedialog.askopenfilename(title="Select Input File", filetypes=filetypes)
        if filename:
            if os.path.getsize(filename) == 0:
                messagebox.showerror("Error", "Selected file is empty!")
                return
            self.input_file = filename
            self.file_label.config(text=os.path.basename(filename))
            ext = os.path.splitext(filename)[1].lower()
            if ext in [".xls", ".xlsx"]:
                try:
                    sheets = pd.ExcelFile(filename).sheet_names
                    self.sheet_list = sheets
                    self.sheet_listbox.delete(0, tk.END)
                    for sheet in sheets:
                        self.sheet_listbox.insert(tk.END, sheet)
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to read Excel sheets: {e}")
            else:
                self.sheet_listbox.delete(0, tk.END)
                self.sheet_list = []
    
    def select_output_folder(self):
        folder = filedialog.askdirectory(title="Select Output Folder")
        if folder:
            self.output_folder = folder
            self.out_folder_label.config(text=os.path.basename(folder), fg="lightgreen")
    
    def run_analysis(self):
        if not self.input_file:
            messagebox.showerror("Error", "Please select an input file.")
            return
        try:
            water_table = float(self.entries["Water Table (m):"].get())
            Mw = float(self.entries["Mw:"].get())
            PGA = float(self.entries["PGA (g):"].get())
            rod_extension = float(self.entries["Rod Extension (m):"].get())
            borehole_diameter = float(self.entries["Borehole Diameter (mm):"].get())
            ER = float(self.entries["Energy Ratio (%):"].get())
            anal_depth = float(self.entries["Analysis Depth (m):"].get())
            xi_v_value = float(self.entries["xi_v (for LSN):"].get())
        except Exception as e:
            messagebox.showerror("Error", f"Invalid parameter value: {e}")
            return
        
        if not self.output_folder:
            self.select_output_folder()
            if not self.output_folder:
                messagebox.showerror("Error", "No output folder selected.")
                return
        
        out_paths = create_output_folders(self.output_folder)
        ext = os.path.splitext(self.input_file)[1].lower()
        self.results = []
        fig_to_display = None
        
        try:
            if ext == ".csv":
                df = pd.read_csv(self.input_file)
                res = process_sheet(df, "CSV_Sheet", out_paths, water_table, Mw, PGA, rod_extension,
                                      borehole_diameter, ER, anal_depth, xi_v_value)
                if res is not None:
                    self.results.append(res)
                    fig_to_display = res["figure"]
            elif ext in [".xls", ".xlsx"]:
                sheets = pd.read_excel(self.input_file, sheet_name=None)
                selected = self.sheet_listbox.curselection()
                if selected:
                    for i in selected:
                        sheet_name = self.sheet_list[i]
                        df = sheets[sheet_name]
                        try:
                            df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
                        except Exception as e:
                            messagebox.showerror("Error", f"Sheet {sheet_name} missing required columns: {e}")
                            continue
                        res = process_sheet(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                                              borehole_diameter, ER, anal_depth, xi_v_value)
                        if res is not None:
                            self.results.append(res)
                            fig_to_display = res["figure"]
                else:
                    for sheet_name, df in sheets.items():
                        try:
                            df = df[["top_depth", "bottom_depth", "NSPT", "FC", "unit_weight"]]
                        except Exception as e:
                            messagebox.showerror("Error", f"Sheet {sheet_name} missing required columns: {e}")
                            continue
                        res = process_sheet(df, sheet_name, out_paths, water_table, Mw, PGA, rod_extension,
                                              borehole_diameter, ER, anal_depth, xi_v_value)
                        if res is not None:
                            self.results.append(res)
                            fig_to_display = res["figure"]
                        break
            else:
                messagebox.showerror("Error", "Unsupported file type. Please provide CSV or XLSX.")
                return
        except Exception as e:
            messagebox.showerror("Error", f"Error processing file: {e}")
            return
        
        if fig_to_display is not None:
            for widget in self.fig_frame.winfo_children():
                widget.destroy()
            self.figure_canvas = FigureCanvasTkAgg(fig_to_display, master=self.fig_frame)
            self.figure_canvas.draw()
            self.figure_canvas.get_tk_widget().pack(fill="both", expand=True)
            self.canvas.config(scrollregion=self.canvas.bbox("all"))
        
        if self.results:
            export_path = os.path.join(self.output_folder, "Combined_Analysis_Results.xlsx")
            export_results_to_excel(self.results, export_path)
            print(f"Exported combined results to {export_path}")
        
        messagebox.showinfo("Success", "Analysis completed. Check the output folder for results.")

# ---------------------------------------------------------
# Main Application with Notebook (Tabs for SPT and CPT)
# ---------------------------------------------------------
class CombinedLPIApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Liquefaction Analysis Application")
        self.geometry("1050x700")
        self.configure(bg="#2c3e50")
        self.create_widgets()
    
    def create_widgets(self):
        notebook = ttk.Notebook(self)
        notebook.pack(fill="both", expand=True)
        spt_tab = SPTFrame(notebook)
        cpt_tab = CPTFrame(notebook)
        notebook.add(spt_tab, text="SPT Analysis")
        notebook.add(cpt_tab, text="CPT Analysis")

if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    app = CombinedLPIApp()
    app.mainloop()
