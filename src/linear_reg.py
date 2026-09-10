import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Raw Data Constants
X_DATA = np.array([0.06698729810778929, 0.25, 0.5, 0.75, 0.9330127018922625])
EXP_Y_DATA = np.array([0.527576297466166, 0.6419620607599249, 0.862154655100411, 0.9622421979824501, 1.146605134114744])
EXP_SIGMA = np.array([0.13379016000002, 0.14838547909538, 0.17356698767590, 0.18439605617689, 0.20269902398342])
CFD_Y_DATA = np.array([0.383698983, 0.64161545, 0.857249858, 1.039745268, 1.16709476])

# Parameter Space Generator
def generate_parameter_space(a_opt, b_opt, a_span=0.5, b_span=0.5, resolution=500):
    """Generate dynamic parameter grid for slope (a) and intercept (b)."""
    a_min = a_opt - a_span
    a_max = a_opt + a_span
    b_min = b_opt - b_span
    b_max = b_opt + b_span
    
    a_values = np.linspace(a_min, a_max, resolution)
    b_values = np.linspace(b_min, b_max, resolution)
    A_GRID, B_GRID = np.meshgrid(a_values, b_values)
    return A_GRID, B_GRID

# Analytical Solvers 
def analytical_wls(x, y, weights):
    """
    Calculate exact WLS parameters using closed-form formulas.
    Returns optimal a, b, and the true maximum R^2 value.
    """
    w_sum = np.sum(weights)
    
    # Weighted means
    x_bar = np.sum(weights * x) / w_sum
    y_bar = np.sum(weights * y) / w_sum
    
    # Calculate slope (a) using wls
    numerator = np.sum(weights * (x - x_bar) * (y - y_bar))
    denominator = np.sum(weights * (x - x_bar)**2)
    
    a_opt = numerator / denominator
    b_opt = y_bar - a_opt * x_bar
    
    # Calculate exact true maximum R^2
    ss_tot_weighted = np.sum(weights * (y - y_bar)**2)
    y_pred = a_opt * x + b_opt
    ss_res_weighted = np.sum(weights * (y - y_pred)**2)
    
    r2_max = 1.0 - (ss_res_weighted / ss_tot_weighted) if ss_tot_weighted > 0 else 0.0
    
    return a_opt, b_opt, r2_max

# Analytical Solvers
def pivot_rule_bounds(x, y, sigma, weights):
    """
    Calculate theoretical slope bounds based on error bars bounding the centroid.
    Uses only first (leftmost) and last (rightmost) points.
    Returns a_min, a_max, b_min, b_max for the smallest valid interval.
    """
    w_sum = np.sum(weights)
    
    # Weighted centroid
    x_centroid = np.sum(weights * x) / w_sum
    y_centroid = np.sum(weights * y) / w_sum
    
    # First point (leftmost) and last point (rightmost)
    x_left, y_left, sigma_left = x[0], y[0], sigma[0]
    x_right, y_right, sigma_right = x[-1], y[-1], sigma[-1]
    
    # Calculate dx from centroid to each extreme point
    dx_L = x_left - x_centroid
    dx_R = x_right - x_centroid
    
    # Calculate slopes from centroid to error boundaries of left point
    slope_L_min = (y_left - sigma_left - y_centroid) / dx_L if dx_L != 0 else 0.0
    slope_L_max = (y_left + sigma_left - y_centroid) / dx_L if dx_L != 0 else 0.0
    
    # Calculate slopes from centroid to error boundaries of right point
    slope_R_min = (y_right - sigma_right - y_centroid) / dx_R if dx_R != 0 else 0.0
    slope_R_max = (y_right + sigma_right - y_centroid) / dx_R if dx_R != 0 else 0.0
    
    # Return the smallest interval across all bounds
    a_min = max(slope_L_min, slope_R_min)
    a_max = min(slope_L_max, slope_R_max)
    
    # Calculate corresponding intercepts 
    b_min = y_centroid - a_max * x_centroid  
    b_max = y_centroid - a_min * x_centroid 
    
    return a_min, a_max, b_min, b_max


# R^2 Grid Calculator
def compute_r2_surface(x, y, weights, a_grid, b_grid):
    """Calculate R^2 surface across the entire parameter grid."""
    w_sum = np.sum(weights)
    y_bar = np.sum(weights * y) / w_sum
    
    # Total Sum of Squares (weighted)
    ss_tot = np.sum(weights * (y - y_bar)**2)
    
    # Calculate predictions for entire grid using broadcasting
    preds = a_grid[:, :, None] * x[None, None, :] + b_grid[:, :, None]
    # Residual Sum of Squares (weighted)
    ss_res = np.sum(weights[None, None, :] * (y[None, None, :] - preds)**2, axis=2)
    
    # Calculate R^2
    r2_array = 1.0 - (ss_res / ss_tot)
    
    # Clip the negative values
    r2_array = np.clip(r2_array, a_min=0, a_max=1)
    
    return r2_array

# Visualization

def analyze_and_plot(a_grid, b_grid, r2_surface, a_opt, b_opt, r2_opt, 
                     pivot_bounds=None, data_label="Data", output_dir="."):

    plt.rcParams['font.family'] = 'serif'
    fig, ax = plt.subplots(figsize=(8, 6))

    # Dynamically extract extent from grid min/max values
    extent = [a_grid.min(), a_grid.max(), b_grid.min(), b_grid.max()]
    
    # Plot heatmap
    im = ax.imshow(r2_surface, extent=extent, origin='lower', aspect='auto', 
                   cmap='viridis', alpha=0.8)
    
    # R^2 contours
    contour_levels = [0.90, 0.95, 0.98]
    cs = ax.contour(a_grid, b_grid, r2_surface, levels=contour_levels, 
                    colors='white', linewidths=1, alpha=0.7)
    ax.clabel(cs, inline=True, fontsize=9, fmt=r'$%.2f$')
    
    # Analytical optimum
    ax.plot(a_opt, b_opt, 'rx', markersize=8, markeredgewidth=1.5, 
            label=f'Optimum: a={a_opt:.4f}, b={b_opt:.4f}, $R^2$={r2_opt:.6f}')
    
    
        # Add pivot rule bounds if available
    if pivot_bounds is not None:
        a_min_bound, a_max_bound = pivot_bounds[0], pivot_bounds[1]
        
        # Vertical lines for slope bounds (a) 
        ax.axvline(x=a_min_bound, color='blue', linestyle='--', linewidth=1.5, 
                   alpha=0.7)#
        ax.axvline(x=a_max_bound, color='blue', linestyle='--', linewidth=1.5, 
                   alpha=0.7, label='Pivot Bounds')  
        
        # Horizontal lines for intercept bounds (b)
        if len(pivot_bounds) >= 4:
            b_min_bound, b_max_bound = pivot_bounds[2], pivot_bounds[3]
            ax.axhline(y=b_min_bound, color='blue', linestyle='--', linewidth=1.5, 
                       alpha=0.7)
            ax.axhline(y=b_max_bound, color='blue', linestyle='--', linewidth=1.5, 
                       alpha=0.7)
            
    # Labels and title
    ax.set_xlabel('Slope parameter (a)', fontsize=12)
    ax.set_ylabel('Intercept parameter (b)', fontsize=12)
    ax.set_title(f'{data_label} - Parameter Space Heatmap', 
                 fontsize=14)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('$R^2$ Value', rotation=90, labelpad=15)
    
    # Legend
    ax.legend(loc='lower right', fontsize=8)
    
    # Grid
    ax.grid(True, alpha=0.3)
    
    # Save figure
    output_path = Path(output_dir) / f'{data_label}_parameter_space_heatmap.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved heatmap: {output_path}")
    plt.close(fig)
    
    # Flatten grids and save to .dat CSV file
    a_flat = a_grid.flatten()
    b_flat = b_grid.flatten()
    r2_flat = r2_surface.flatten()
    
    data_dict = {'a': a_flat, 'b': b_flat, 'R2': r2_flat}
    
    # execution directory
    csv_path = Path(output_dir) / f'{data_label}_parameter_space.dat'
    
    df = pd.DataFrame({'a': a_flat, 'b': b_flat, 'R2': r2_flat})
    df.to_csv(csv_path, index=False, float_format='%.6f')
    
    print(f"Saved CSV data: {csv_path}")



# Execute Pipelines
if __name__ == "__main__":
    # Get script execution directory
    script_dir = Path(__file__).parent
    
    # Experimental Data (WLS with 1/sigma^2 weights)
    print("EXPERIMENTAL DATA ANALYSIS")
    
    exp_weights = 1.0 / (EXP_SIGMA ** 2)
    
    # Run WLS analytical solver
    a_exp_opt, b_exp_opt, r2_exp_opt = analytical_wls(X_DATA, EXP_Y_DATA, exp_weights)
    print(f"Analytical Optimum: a={a_exp_opt:.6f}, b={b_exp_opt:.6f}")
    print(f"Maximum $R^2$: {r2_exp_opt:.6f}")
    
    # Generate dynamic grid for experimental data
    A_GRID_EXP, B_GRID_EXP = generate_parameter_space(a_exp_opt, b_exp_opt)
    
    # Run R^2 grid calculator
    r2_exp_surface = compute_r2_surface(X_DATA, EXP_Y_DATA, exp_weights, A_GRID_EXP, B_GRID_EXP)
    
    # Extract pivot bounds
    pivot_bounds_exp = pivot_rule_bounds(X_DATA, EXP_Y_DATA, EXP_SIGMA, exp_weights)
    print(f"Pivot Rule Bounds:")
    print(f"a_min={pivot_bounds_exp[0]:.6f}, b_min={pivot_bounds_exp[2]:.6f}")
    print(f"a_max={pivot_bounds_exp[1]:.6f}, b_max={pivot_bounds_exp[3]:.6f}")
    
    # Send to plotter
    analyze_and_plot(A_GRID_EXP, B_GRID_EXP, r2_exp_surface, a_exp_opt, b_exp_opt, 
                     r2_exp_opt, pivot_bounds=pivot_bounds_exp,  
                     data_label="Experimental", output_dir=script_dir)
    
    # CFD Data (OLS with weights = 1.0)
    print("CFD DATA ANALYSIS")
    
    cfd_weights = np.ones_like(CFD_Y_DATA)
    cfd_sigma = np.zeros_like(CFD_Y_DATA)
    
    # Run analytical solver (OLS equivalent with uniform weights)
    a_cfd_opt, b_cfd_opt, r2_cfd_opt = analytical_wls(X_DATA, CFD_Y_DATA, cfd_weights)
    print(f"Analytical Optimum: a={a_cfd_opt:.6f}, b={b_cfd_opt:.6f}")
    print(f"Maximum $R^2$: {r2_cfd_opt:.6f}")
    
    # Generate dynamic grid for CFD data
    A_GRID_CFD, B_GRID_CFD = generate_parameter_space(a_cfd_opt, b_cfd_opt)
    
    # Run R^2 grid calculator
    r2_cfd_surface = compute_r2_surface(X_DATA, CFD_Y_DATA, cfd_weights, A_GRID_CFD, B_GRID_CFD)
    
    # Send to plotter (skipping pivot bounds for CFD)
    analyze_and_plot(A_GRID_CFD, B_GRID_CFD, r2_cfd_surface, a_cfd_opt, b_cfd_opt, 
                     r2_cfd_opt, pivot_bounds=None, data_label="CFD", output_dir=script_dir)
    
    print("ANALYSIS COMPLETE")