import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#conf.

RESOLUTION = 500
A_RANGE = np.linspace(0, 1.5, RESOLUTION)
B_RANGE = np.linspace(0, 1.0, RESOLUTION)
A_GRID, B_GRID = np.meshgrid(A_RANGE, B_RANGE)

# Shared x values (sin^2(alpha/2))
x_vals = np.array([
    0.06698729810778929, 0.25, 0.5, 0.75, 0.9330127018922625
], dtype=float)

# helper functions

def compute_analytical_wls(x, y, weights):
    """Computes closed-form a, b for weighted least squares."""
    # Weighted means
    sum_w = np.sum(weights)
    x_w_mean = np.sum(weights * x) / sum_w
    y_w_mean = np.sum(weights * y) / sum_w
    
    # Numerator and Denominator for slope 'a'
    num = np.sum(weights * (x - x_w_mean) * (y - y_w_mean))
    den = np.sum(weights * (x - x_w_mean)**2)
    
    a_opt = num / den
    b_opt = y_w_mean - a_opt * x_w_mean
    return a_opt, b_opt

def compute_pivot_rule_bounds(x, y, sigma, weights):
    """
    Calculates a_min and a_max using the Pivot Rule (End Points Method).
    1. Finds weighted centroid.
    2. Uses ONLY the leftmost and rightmost points to determine slope limits.
    """
    # 1. Calculate Weighted Centroid
    sum_w = np.sum(weights)
    x_cent = np.sum(weights * x) / sum_w
    y_cent = np.sum(weights * y) / sum_w
    
    # 2. Identify Extreme Points
    idx_min = np.argmin(x)
    idx_max = np.argmax(x)
    
    # Leftmost Point (L)
    x_L = x[idx_min]
    y_L = y[idx_min]
    sig_L = sigma[idx_min]
    
    # Rightmost Point (R)
    x_R = x[idx_max]
    y_R = y[idx_max]
    sig_R = sigma[idx_max]
    
    dx_L = x_L - x_cent
    dx_R = x_R - x_cent

    # 3. Calculate Limits based on End Points
    # Max Slope: Centroid to Top of Right OR Bottom of Left
    # We take the tighter constraint (minimum of the calculated max-slopes)
    # Slope to Top Right: (y_R + sig_R - y_c) / dx_R
    # Slope to Bottom Left: (y_L - sig_L - y_c) / dx_L  (Note: dx_L is negative)
    
    slope_max_R = (y_R + sig_R - y_cent) / dx_R
    slope_max_L = (y_L - sig_L - y_cent) / dx_L 
    # Since dx_L is negative, (y_low - y_c)/neg is a positive slope upper bound
    
    a_pivot_max = min(slope_max_R, slope_max_L)
    
    # Min Slope: Centroid to Bottom of Right OR Top of Left
    slope_min_R = (y_R - sig_R - y_cent) / dx_R
    slope_min_L = (y_L + sig_L - y_cent) / dx_L
    
    a_pivot_min = max(slope_min_R, slope_min_L)
    
    return a_pivot_min, a_pivot_max

def compute_r2_grid(x, y, weights):
    """Computes R2 surface over the global A_GRID, B_GRID."""
    y_w_mean = np.sum(weights * y) / np.sum(weights)
    ss_tot = np.sum(weights * (y - y_w_mean)**2)
    preds = A_GRID[:, :, np.newaxis] * x + B_GRID[:, :, np.newaxis]
    residuals = y - preds
    ss_res = np.sum(weights * (residuals**2), axis=2)
    return 1 - (ss_res / ss_tot)

def analyze_and_plot(name, x, y, weights, sigma, r2_grid, output_prefix, dat_filename, do_pivot_bounds=True):
    """Finds optima, extracts contours, plots heatmap, saves data."""
    print(f"\n--- ANALYSIS: {name} ---")
    
    # 1. Grid Optimum
    max_idx = np.unravel_index(np.argmax(r2_grid), r2_grid.shape)
    grid_a = A_RANGE[max_idx[1]]
    grid_b = B_RANGE[max_idx[0]]
    grid_r2 = r2_grid[max_idx]
    
    # 2. Analytical Optimum
    ana_a, ana_b = compute_analytical_wls(x, y, weights)
    
    print(f"Grid Max: a={grid_a:.4f}, b={grid_b:.4f}, R2={grid_r2:.5f}")
    print(f"Analytical: a={ana_a:.4f}, b={ana_b:.4f}")

    # 3. Pivot Rule Calculation
    if do_pivot_bounds:
     pivot_min, pivot_max = compute_pivot_rule_bounds(x, y, sigma, weights)
     print(f"Pivot Rule Bounds (End Points): a_min={pivot_min:.4f}, a_max={pivot_max:.4f}")
    else:
     print("Pivot Rule Bounds: skipped by configuration.")
    
    # 4. Plotting
    plt.rcParams['font.family'] = 'serif'
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Heatmap
    im = ax.imshow(r2_grid, extent=[0, 1.5, 0, 1], origin='lower',
                   aspect='auto', cmap='viridis', vmin=0, vmax=1)
    
    cbar_label = r'Weighted $R^2$' if "Experimental" in name else r'Ordinary $R^2$'
    plt.colorbar(im, label=cbar_label)
    
    # Contours
    levels = [0.90, 0.95, 0.98, 0.99]
    contours = ax.contour(A_GRID, B_GRID, r2_grid, levels=levels, 
                          colors='white', linestyles='dashed', linewidths=1)
    ax.clabel(contours, inline=True, fontsize=10, fmt='%1.2f')
    
    # Analytical Optimum (Red Cross)
    ax.plot(ana_a, ana_b, 'rx', markersize=10, markeredgewidth=2, 
            label='Max $R^2$ (Analytical)')
    
        
    ax.set_xlabel('Slope Parameter $a$')
    ax.set_ylabel('Intercept Parameter $b$')
    ax.set_title(f'{name} Parameter Space Heatmap')
    ax.legend(loc='lower right')
    
    img_name = f"{output_prefix}_param_space.png"
    plt.tight_layout()
    plt.savefig(img_name, dpi=300)
    print(f"Saved plot: {img_name}")
    
    # 5. Save Data to .dat
    df = pd.DataFrame({
        'a': A_GRID.flatten(),
        'b': B_GRID.flatten(),
        'R2': r2_grid.flatten()
    })
    df.to_csv(dat_filename, index=False, float_format='%.6f')
    print(f"Saved data: {dat_filename}")


# 1.exp analysis (WLS)

y_exp = np.array([
    0.527576297466166, 0.6419620607599249, 0.862154655100411, 
    0.9622421979824501, 1.146605134114744
], dtype=float)

sigma_exp = np.array([
    0.13379016000002, 0.14838547909538, 0.17356698767590, 
    0.18439605617689, 0.20269902398342
], dtype=float)

# Weights = 1/sigma^2 for WLS
w_exp = 1.0 / (sigma_exp ** 2)
r2_exp = compute_r2_grid(x_vals, y_exp, w_exp)

analyze_and_plot(
    name="Experimental (WLS)", 
    x=x_vals, 
    y=y_exp, 
    weights=w_exp, 
    sigma=sigma_exp,
    r2_grid=r2_exp, 
    output_prefix="exp_a_b", 
    dat_filename="exp_a(0_1.5)_b(0_1)_param_space_heatmap_data.dat",
    do_pivot_bounds=True
)

# 2. cfd analysis (WLS with weight=1)
y_cfd = np.array([
    0.383698983, 0.64161545, 0.857249858, 
    1.039745268, 1.16709476
], dtype=float)

# OLS implies no uncertainty weighting
sigma_cfd = np.zeros_like(y_cfd)
w_cfd = np.ones_like(y_cfd)
r2_cfd = compute_r2_grid(x_vals, y_cfd, w_cfd)

analyze_and_plot(
    name="CFD (OLS)", 
    x=x_vals, 
    y=y_cfd, 
    weights=w_cfd,
    sigma=sigma_cfd, 
    r2_grid=r2_cfd, 
    output_prefix="cfd_a_b", 
    dat_filename="cfd_a(0_1.5)_b(0_1)_param_space_heatmap_data.dat",
    do_pivot_bounds=False
)