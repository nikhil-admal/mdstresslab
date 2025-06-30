##
# @file visualize.py
# @brief Generates 2D contour plots of stress components from .stress files.
#
# This script reads atomistic simulation output and visualizes a selected
# stress component (e.g., σ<sub>yy</sub>) using a filled contour plot.
#
# ## Input Format
# The input file should be structured as follows:
#
# ```
# 125
# Properties=pos:R:3:stress:R:6
# x y z σ_xx σ_yy σ_zz σ_xy σ_xz σ_yz
# ...
# ```
#
# Each data line must contain 9 floating-point values: the atomic position `(x, y, z)`
# followed by the six unique components of the symmetric Cauchy stress tensor.
#
# ## Command-Line Usage
# ```
# python visualize.py <filename.stress> <component: xx|yy|zz|xy|xz|yz>
# ```
#
# Converts the selected stress component from internal units to GPa
# (by multiplying by 160.0), and saves a PDF contour plot.
#
# ## Related Pages
# - @ref visualization_intro "Visualization Utilities Overview"
#

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

##
# @brief Mapping from stress component names to (column index, label).
STRESS_COMPONENTS = {
    'xx': (3, r'$\sigma_{xx}$'),
    'yy': (4, r'$\sigma_{yy}$'),
    'zz': (5, r'$\sigma_{zz}$'),
    'xy': (6, r'$\sigma_{xy}$'),
    'xz': (7, r'$\sigma_{xz}$'),
    'yz': (8, r'$\sigma_{yz}$'),
}

##
# @brief Reads a .stress file containing atomistic stress data.
# @param filename The path to the .stress file.
# @return A list of rows with [x, y, z, σ_xx, σ_yy, σ_zz, σ_xy, σ_xz, σ_yz].
#
# The file is expected to have:
# - Line 1: number of atoms (ignored)
# - Line 2: metadata (e.g., Properties=...)
# - Line 3 onward: whitespace-delimited numeric data (one line per atom)
#
# Malformed or non-numeric lines are skipped.
def read_file(filename):
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[2:]:  # Skip count and metadata header
            line = line.strip()
            if not line:
                continue
            try:
                parts = list(map(float, line.split()))
                if len(parts) >= 9:
                    data.append(parts[:9])
            except ValueError:
                continue
    return data

##
# @brief Plots a 2D contour of the selected stress component.
# @param data List of rows with [x, y, z, σ_xx, ..., σ_yz].
# @param output_fname Output PDF filename.
# @param component_index Index of the stress component to plot.
# @param component_label Label for the plot title.
#
# Uses SciPy's `griddata` to interpolate the scattered data, and matplotlib
# to create a filled contour plot with 100 levels using the RdBu_r colormap.
def plot_contour(data, output_fname, component_index, component_label):
    x = [row[0] for row in data]
    y = [row[1] for row in data]  # using x-y plane for plotting
    sig = [row[component_index] for row in data]

    xi = np.linspace(min(x), max(x), 300)
    yi = np.linspace(min(y), max(y), 300)
    xi, yi = np.meshgrid(xi, yi)

    sig_grid = griddata((x, y), sig, (xi, yi), method='linear')

    plt.figure(figsize=(6, 5))
    contour = plt.contourf(xi, yi, sig_grid, levels=100, cmap='RdBu_r')
    plt.colorbar(contour, label='Stress (GPa)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(component_label)
    plt.tight_layout()
    plt.savefig(output_fname)
    plt.close()

# Main script execution
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python visualize.py <stress_file> <component: xx|yy|zz|xy|xz|yz>")
        sys.exit(1)

    filename = sys.argv[1]
    component_key = sys.argv[2].lower()

    if component_key not in STRESS_COMPONENTS:
        print("Component must be one of: xx, yy, zz, xy, xz, yz")
        sys.exit(1)

    component_index, component_label = STRESS_COMPONENTS[component_key]

    base, _ = os.path.splitext(filename)
    plt_fname = f"{base}_{component_key}.pdf"

    data = read_file(filename)

    # Convert stress to GPa
    for row in data:
        row[component_index] *= 160.0

    plot_contour(data, plt_fname, component_index, component_label)

