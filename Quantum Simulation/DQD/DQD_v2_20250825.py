import numpy as np
import matplotlib.pyplot as plt
from qutip import Qobj

def H_sw_transformation(epsilon, U, t, dEz, Ez): 
    """
    Creates an effective 4x4 Hamiltonian based on a Schrieffer-Wolff transformation.
    The basis is {|↓↓⟩, |↓↑⟩, |↑↓⟩, |↑↑⟩}.

    Args:
        epsilon (float): Detuning energy in µeV.
        U (float): Charging energy in µeV.
        t (float): Interdot tunnel coupling in µeV.
        dEz (float): Difference in Zeeman splitting in µeV.
        Ez (float): Average Zeeman energy in µeV.

    Returns:
        Qobj: The effective Hamiltonian matrix as a QuTiP object.
    """
    # Add a very small number to denominators to prevent division by zero
    safe_epsilon = 1e-9

    alpha_dEz = t**2/(U - epsilon - dEz/2 + safe_epsilon) + t**2/(U + epsilon - dEz/2 + safe_epsilon)
    alpha_minus_dEz = t**2/(U - epsilon + dEz/2 + safe_epsilon) + t**2/(U + epsilon + dEz/2 + safe_epsilon)
    beta = 0.5*(alpha_dEz + alpha_minus_dEz)

    H_matrix = np.array([
        #    |↓↓⟩            |↓↑⟩                 |↑↓⟩                 |↑↑⟩      
        [    -Ez,              0,                     0,                    0      ], # |↓↓⟩
        [      0, -dEz/2 - alpha_minus_dEz,          beta,                 0      ], # |↓↑⟩
        [      0,           beta,           dEz/2 - alpha_dEz,          0      ], # |↑↓⟩
        [      0,              0,                     0,                   Ez     ]  # |↑↑⟩
    ], dtype=np.complex128)
    
    return Qobj(H_matrix)

# --- Simulation Parameters ---
# Using a larger `t` makes the effects of the SW transformation more visible
U = 1000.0   # µeV (Charging energy)
t = 100.0    # µeV (Tunnel coupling, increased for visibility)
Ez = 100.0   # µeV (Average Zeeman energy)

# --- Set Fixed and Sweeping Parameters ---
# Fix epsilon at 0 as requested
epsilon_fixed = 0.0

# Define the range for dEz to be swept, on a logarithmic scale
dEz_values = np.logspace(-4, 2, 500)

# --- State Tracking Initialization ---
# Get the initial eigenvalues and eigenvectors at the start of the dEz range
H_initial = H_sw_transformation(epsilon_fixed, U, t, dEz_values[0], Ez)
initial_evals, initial_evecs = H_initial.eigenstates()

# Prepare a list of lists to store the energy of each tracked state
tracked_energies = [[energy] for energy in initial_evals]
previous_evecs = initial_evecs

# --- Main State Tracking Loop ---
# Iterate through the rest of the dEz values
for dEz_val in dEz_values[1:]:
    H = H_sw_transformation(epsilon_fixed, U, t, dEz_val, Ez)
    evals, evecs = H.eigenstates()
    
    # Create an overlap matrix between previous and current eigenvectors
    overlap_matrix = np.abs([v_old.overlap(v_new) for v_old in previous_evecs for v_new in evecs]).reshape(len(evals), len(evals))
    
    # Find the permutation of new states that best matches the old states
    permutation = np.argmax(overlap_matrix, axis=1)
    
    # Reorder the results according to the tracked states
    reordered_evals = evals[permutation]
    reordered_evecs = [evecs[i] for i in permutation]
    
    # Append the new energies to the correct tracking list
    for i in range(len(reordered_evals)):
        tracked_energies[i].append(reordered_evals[i])
        
    previous_evecs = reordered_evecs

# --- Plotting the Final Result ---
plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(10, 7))

# Define labels based on the character of the states. 
# The order corresponds to the initial energy levels from lowest to highest.
labels = [
    "|↓↓⟩",
    "|↓↑⟩",
    "|↑↓⟩",
    "|↑↑⟩"
]

for i, energies in enumerate(tracked_energies):
    plt.plot(dEz_values, energies, lw=2.5, label=labels[i])

# --- Formatting the plot ---
plt.xscale('log')
plt.xlabel('Zeeman Difference $\Delta E_z$ (µeV)', fontsize=14)
plt.ylabel('Energy (µeV)', fontsize=14)
plt.title(f'Energy Levels vs. Zeeman Difference at ε = {epsilon_fixed} µeV', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12, title="Eigenstate")
plt.grid(True, which="both", ls="-")
plt.show()