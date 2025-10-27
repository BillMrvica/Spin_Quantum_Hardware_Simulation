import numpy as np
import matplotlib.pyplot as plt
from qutip import Qobj, basis

def create_hamiltonian(epsilon, U, t, dEz, Ez):
    """
    Creates the Hamiltonian from Eq. (2) of the paper.
    The basis is {|↑↑⟩, |↑↓⟩, |↓↑⟩, |↓↓⟩, S(0,2), S(2,0)}.

    Args:
        epsilon (float): Detuning energy in µeV.
        U (float): Charging energy in µeV.
        t (float): Interdot tunnel coupling in µeV.
        dEz (float): Difference in Zeeman splitting in µeV.
        Ez (float): Average Zeeman energy in µeV.

    Returns:
        Qobj: The Hamiltonian matrix as a QuTiP object.
    """
    # Hamiltonian matrix directly from Eq. (2) in the paper
    # Note: The basis in the paper is {T+, |↑↓>, |↓↑>, T-, S(0,2), S(2,0)}
    # where T+ = |↑↑> and T- = |↓↓>.
    H_matrix = np.array([
        #    |↓↓⟩    |↓↑⟩      |↑↓⟩      |↑↑⟩      S(0,2)       S(2,0)
        [    -Ez,      0,        0,        0,          0,           0         ], # |↑↑⟩
        [      0, -dEz/2,        0,        0,          t,           t         ], # |↑↓⟩
        [      0,      0,    dEz/2,        0,         -t,          -t         ], # |↓↑⟩
        [      0,      0,        0,       Ez,          0,           0         ], # |↓↓⟩
        [      0,      t,       -t,        0,  U - epsilon,         0         ], # S(0,2)
        [      0,      t,       -t,        0,          0,     U + epsilon     ]  # S(2,0)
    ], dtype=np.complex128)
    
    return Qobj(H_matrix)

# Parameters exactly as specified in the caption of Fig. 1(a)
U = 1000.0   # µeV
t = 10.0     # µeV
dEz = 14.0   # µeV
Ez = 100.0   # µeV

# Range of detuning values to plot
epsilon_values = np.linspace(-1500, 1500, 400)

# --- State Tracking Initialization ---
# Get the initial eigenvalues and eigenvectors at the start of the range
H_initial = create_hamiltonian(epsilon_values[0], U, t, dEz, Ez)
initial_evals, initial_evecs = H_initial.eigenstates()

# Prepare a list of lists to store the energy of each tracked state
tracked_energies = [[energy] for energy in initial_evals]
previous_evecs = initial_evecs

# --- Main State Tracking Loop ---
# Iterate through the rest of the detuning values
for eps in epsilon_values[1:]:
    H = create_hamiltonian(eps, U, t, dEz, Ez)
    evals, evecs = H.eigenstates()
    
    # Create an overlap matrix between previous and current eigenvectors
    # overlap_matrix[i, j] is the overlap of the i-th old state with the j-th new state
    overlap_matrix = np.abs([v_old.overlap(v_new) for v_old in previous_evecs for v_new in evecs]).reshape(len(evals), len(evals))
    
    # Find the permutation of new states that best matches the old states
    # This assumes no two states swap their order in a single small step
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
plt.figure(figsize=(12, 8))

# Define labels based on the character of the states at epsilon = 0
# This order corresponds to the initial energy levels from lowest to highest
labels = [
    "|↓↑>",
    "|↓↓>",
    "|↑↓>",
    "S(2,0)",
    "|↑↑>",
    "S(0,2)"
]

for i, energies in enumerate(tracked_energies):
    plt.plot(epsilon_values, energies, lw=2.5, label=labels[i])

# --- Formatting the plot to match Figure 1(a) ---
plt.xlabel('Detuning ε (µeV)', fontsize=14)
plt.ylabel('Energy (µeV)', fontsize=14)
plt.title('Recreation of Energy Level Diagram via State Tracking', fontsize=16)
plt.ylim(-200, 200)
plt.xlim(-1500, 1500)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12, title="State Character (at ε ≈ 0)")
plt.show()