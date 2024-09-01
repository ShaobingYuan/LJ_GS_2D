# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 19:33:29 2024

@author: shaob
"""

import numpy as np
from numba import njit
import matplotlib.pyplot as plt

@njit
def generate_initial_configuration(N):
    # Generate N cities randomly in a 1x1 square box
    initial_coordinates = np.random.rand(N, 2)  # Each city is represented by its (x, y) coordinates
    return initial_coordinates

@njit
def LJ_potential(sigma, coordinates, i, j):
    distance = np.sqrt(np.sum((coordinates[i] - coordinates[j]) ** 2))  # Calculate Euclidean distance
    
    # Find the shortest distance under periodic boundary condition
    periodic_displacement = np.array([[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [1, -1], [-1, 1], [-1, -1]])
    for k in range(8):
        new_distance = np.sqrt(np.sum((coordinates[i] + periodic_displacement[k] - coordinates[j]) ** 2))
        if new_distance < distance:
            distance = new_distance
    
    r = distance / sigma  # Calculate the renormalized distance
    V_LJ = 4 * ( r ** (- 12) - r ** (- 6))  # Calculate the interaction energy
    return V_LJ

@njit
def trial_move(coordinates, step_length, i):
    trial_coordinates = coordinates.copy()
    displacement_x = step_length * np.random.rand()  # Displacement in the x direction
    displacement_y = step_length * np.random.rand()  # Displacement in the y direction
    trial_coordinates[i][0] += displacement_x
    trial_coordinates[i][0] %= 1
    trial_coordinates[i][1] += displacement_y
    trial_coordinates[i][1] %= 1
    return trial_coordinates

@njit
def energy_change(N, sigma, coordinates, trial_coordinates, i):
    # Calculate the energy change
    Delta_E = 0
    for j in range(i):
        Delta_E += LJ_potential(sigma, trial_coordinates, i, j)
        Delta_E -= LJ_potential(sigma, coordinates, i, j)
    for j in range(i + 1, N):
        Delta_E += LJ_potential(sigma, trial_coordinates, i, j)
        Delta_E -= LJ_potential(sigma, coordinates, i, j)
    
    return Delta_E

@njit
def accept_move(Delta_E, T):
    
    # If Delta_E is negative, always accept the move
    if Delta_E < 0:
        return True
    
    # If E is positive, accept the move with probability e^{-E/T}
    if np.exp(- Delta_E / T) > np.random.rand():  # Using NumPy's random number generator for speed and efficiency
        return True
    else:
        return False
    
@njit
def temperature(initial_temperature, cooling_rate, iteration):
    
    current_temperature = initial_temperature * ((1 - cooling_rate) ** iteration)
    
    return current_temperature

@njit
def find_GS(N, sigma, step_length, initial_temperature, cooling_rate):
    coordinates = generate_initial_configuration(N)  # Initial configuration
    moving = True  # Not yet converged
    iteration = 0  # Iteration number
    E = 0
    energy = []
    temp = []
    
    # Let it approach equilibrium
    T = initial_temperature
    for k in range(15*N):
        for i in range(N):
            trial_coordinates = trial_move(coordinates, sigma, i)
            Delta_E = energy_change(N, sigma, coordinates, trial_coordinates, i)
            if accept_move(Delta_E, T):
                coordinates = trial_coordinates.copy()
    for k in range(15*N):
        for i in range(N):
            trial_coordinates = trial_move(coordinates, 10 * step_length, i)
            Delta_E = energy_change(N, sigma, coordinates, trial_coordinates, i)
            if accept_move(Delta_E, T):
                coordinates = trial_coordinates.copy()
    for k in range(15*N):
        for i in range(N):
            trial_coordinates = trial_move(coordinates, 4 * step_length, i)
            Delta_E = energy_change(N, sigma, coordinates, trial_coordinates, i)
            if accept_move(Delta_E, T):
                coordinates = trial_coordinates.copy()

    # Cool it down
    while moving:
        moving = False
        T = temperature(initial_temperature, cooling_rate, iteration)
        for k in range(2):
            for i in range(N):
                trial_coordinates = trial_move(coordinates, step_length, i)
                Delta_E = energy_change(N, sigma, coordinates, trial_coordinates, i)
                if accept_move(Delta_E, T):
                    coordinates = trial_coordinates.copy()
                    E += Delta_E
                    moving = True
                energy.append(E)
                temp.append(T)
        iteration += 1
    
    return coordinates, energy, temp

def plot_GS(coordinates, N, sigma):
    # Create the plot
    plt.figure(figsize=(6, 6))  # Set the figure size to visualize the square region
    plt.axis("equal")  # Equal aspect ratio
    plt.xlim(-0.2, 1.2)  # Set the limits of the x-axis to match the square region
    plt.ylim(-0.2, 1.2)  # Set the limits of the y-axis to match the square region
    
    x, y = coordinates[:, 0], coordinates[:, 1]
    plt.scatter(x, y, zorder=2, s=10)
    
    periodic_displacement_1 = np.array([[-1, 0], [0, -1], [-1, 1], [-1, -1]])
    periodic_displacement_2 = np.array([[1, 0], [0, 1], [1, 1], [1, -1]])
    
    for dx, dy in periodic_displacement_1:
        plt.scatter(x + dx, y + dy, color='red', alpha=0.5, zorder=1, s=10)
    
    for dx, dy in periodic_displacement_2:
        plt.scatter(x + dx, y + dy, color='red', alpha=0.5, zorder=3, s=10)

    plt.title(f'Ground State Configuration with {N} particles, sigma = {sigma:.2g}')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.tick_params(direction='in')
    
    plt.axhline(0, color='grey', linestyle='--', lw=1)
    plt.axhline(1, color='grey', linestyle='--', lw=1)
    plt.axvline(0, color='grey', linestyle='--', lw=1)
    plt.axvline(1, color='grey', linestyle='--', lw=1)
    
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(f"GS_{N}_{sigma:.2g}.pdf", dpi=600)
    plt.show()
    
def plot_energy_change(energy, temp):
    fig, ax1 = plt.subplots()
    
    # Plot energy change
    ax1.plot(energy, 'b-', label='Relative Energy Change')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Relative Energy Change')
    ax1.tick_params(direction='in')

    
    # Plot temperature
    ax2 = ax1.twinx()
    ax2.plot(temp, 'r-', label='Temperature')
    ax2.set_ylabel('Temperature')
    ax2.tick_params(direction='in')
    
    fig.legend(bbox_to_anchor=(0.8,0.9))
    
    # Title
    plt.title('Relative Energy Change While Cooling down')
    
    # Show
    fig.tight_layout()
    plt.savefig(f"Energy_Change_{N}_{sigma:.2g}.pdf", dpi=600)
    plt.show()
    
def GS_average_energy(coordinates, N, sigma):
    E_0 = 0
    for i in range(N):
        for j in range(i + 1, N):
            E_0 += LJ_potential(sigma, coordinates, i, j)
    epsilon_0 = E_0 / N
    return epsilon_0
    
import time

N = 418
sigma = (2 ** (1 / 3)) * (3 ** (- 1 / 4)) * (N ** ( - 0.5))
step_length = 0.025 * sigma
initial_temperature = 0.1
cooling_rate = 0.05

time1 = time.time()
coordinates, energy, temp = find_GS(N, sigma, step_length, initial_temperature, cooling_rate)
time2 = time.time()

print("elapsed time = ", time2 - time1)
    
plot_GS(coordinates, N, sigma)
plot_energy_change(energy, temp)
print(GS_average_energy(coordinates, N, sigma))
