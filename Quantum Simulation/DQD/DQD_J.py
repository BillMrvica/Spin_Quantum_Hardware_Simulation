import numpy as np
import matplotlib.pyplot as plt

def calculate_sw_terms(epsilon, U, t, dEz):
    """
    Calculates the alpha and beta terms from the Schrieffer-Wolff transformation.

    Args:
        epsilon (float or np.ndarray): Detuning energy in µeV.
        U (float): Charging energy in µeV.
        t (float): Interdot tunnel coupling in µeV.
        dEz (float or np.ndarray): Difference in Zeeman splitting in µeV.

    Returns:
        tuple: A tuple containing (alpha_dEz, alpha_minus_dEz, beta).
    """
    # Add a very small number to denominators to prevent division by zero
    safe_epsilon = 1e-9
    
    # Calculate the energy shift terms (alpha)
    alpha_dEz = t**2 / (U - epsilon - dEz/2 + safe_epsilon) + \
                t**2 / (U + epsilon - dEz/2 + safe_epsilon)
                
    alpha_minus_dEz = t**2 / (U - epsilon + dEz/2 + safe_epsilon) + \
                      t**2 / (U + epsilon + dEz/2 + safe_epsilon)
                      
    # Calculate the effective exchange coupling (beta)
    beta = 0.5 * (alpha_dEz + alpha_minus_dEz)
    
    return alpha_dEz, alpha_minus_dEz, beta

# --- Simulation Parameters ---
# Use the same base values as the previous examples
U = 1000.0  # µeV (Charging energy)
t = 100.0   # µeV (Tunnel coupling)

# --- Create Figure with Two Subplots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
plt.style.use('seaborn-v0_8-whitegrid')

# --- Plot 1: Varying Epsilon ---

# Fixed dEz for the first plot
dEz_fixed = 14.0  # µeV

# Epsilon range from 0 to 400 µeV
epsilon_values = np.linspace(0, 400, 500)

# Calculate terms for the epsilon range
alpha_dEz_vs_eps, alpha_minus_dEz_vs_eps, beta_vs_eps = \
    calculate_sw_terms(epsilon_values, U, t, dEz_fixed)

# Plot the results on the first subplot
ax1.plot(epsilon_values, alpha_dEz_vs_eps, lw=2.5, label=r'$\alpha(\Delta E_z)$')
ax1.plot(epsilon_values, alpha_minus_dEz_vs_eps, lw=2.5, label=r'$\alpha(-\Delta E_z)$')
ax1.plot(epsilon_values, beta_vs_eps, lw=3, linestyle='--', color='k', label=r'$\beta$ (Effective Exchange J)')

ax1.set_title(f'Parameters vs. Detuning (ε)\n(Fixed $\Delta E_z$ = {dEz_fixed} µeV)', fontsize=15)
ax1.set_xlabel('Detuning ε (µeV)', fontsize=12)
ax1.set_ylabel('Energy (µeV)', fontsize=12)
ax1.legend(fontsize=11)
ax1.grid(True)

# --- Plot 2: Varying dEz ---

# Fixed epsilon for the second plot
epsilon_fixed = 0.0  # µeV (symmetric point)

# dEz range from 10^-4 to 10^2 µeV on a logarithmic scale
dEz_values = np.logspace(-4, 2, 500)

# Calculate terms for the dEz range
alpha_dEz_vs_dEz, alpha_minus_dEz_vs_dEz, beta_vs_dEz = \
    calculate_sw_terms(epsilon_fixed, U, t, dEz_values)

# Plot the results on the second subplot
ax2.plot(dEz_values, alpha_dEz_vs_dEz, lw=2.5, label=r'$\alpha(\Delta E_z)$')
ax2.plot(dEz_values, alpha_minus_dEz_vs_dEz, lw=2.5, label=r'$\alpha(-\Delta E_z)$')
ax2.plot(dEz_values, beta_vs_dEz, lw=3, linestyle='--', color='k', label=r'$\beta$ (Effective Exchange J)')

ax2.set_xscale('log')
ax2.set_title(f'Parameters vs. Zeeman Difference ($\Delta E_z$)\n(Fixed ε = {epsilon_fixed} µeV)', fontsize=15)
ax2.set_xlabel('Zeeman Difference $\Delta E_z$ (µeV)', fontsize=12)
ax2.set_ylabel('Energy (µeV)', fontsize=12)
ax2.legend(fontsize=11)
ax2.grid(True, which="both", ls="-")

# --- Final Display ---
plt.tight_layout(pad=3.0)
plt.show()