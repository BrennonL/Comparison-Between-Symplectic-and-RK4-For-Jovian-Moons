import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
from horizons_interface import get_reference_ephemerides
from numerical_methods import (
    RK4,
    galilean_derivative,
    stormer_verlet,
    states_time_hist,
)

SECONDS_PER_DAY = 86400
JPL_STEP_DAYS = 10
SIM_YEARS = 10
DAYS_PER_YEAR = 365.25
METHOD_SUBSTEPS_PER_JPL_STEP = 1000

MOON_NAMES = ("Io", "Europa", "Ganymede", "Callisto")

def format_system_state(label, system_state):
    lines = [f"{label}:"]
    for moon_name, state in zip(MOON_NAMES, system_state):
        lines.append(f"{moon_name}: {np.array(state)}")
    return "\n".join(lines)

def compute_position_error_histories(reference_ephemeris, method_ephemeris):
    error_histories = {moon_name: [] for moon_name in MOON_NAMES}

    for reference_state, method_state in zip(reference_ephemeris, method_ephemeris):
        for moon_name, reference_moon, method_moon in zip(MOON_NAMES, reference_state, method_state):
            reference_position = np.array(reference_moon[:3], dtype=float)
            method_position = np.array(method_moon[:3], dtype=float)
            error_histories[moon_name].append(np.linalg.norm(method_position - reference_position))

    return error_histories

def plot_position_errors(time_days, rk4_errors, symplectic_errors):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

    for ax, moon_name in zip(axes.flat, MOON_NAMES):
        ax.plot(time_days, rk4_errors[moon_name], label="RK4")
        ax.plot(time_days, symplectic_errors[moon_name], label="Stormer-Verlet")
        ax.set_title(moon_name)
        ax.set_xlabel("Time [days]")
        ax.set_ylabel("Position error [km]")
        ax.grid(True)
        ax.legend()

    fig.suptitle("Position Error Relative to JPL")
    plt.tight_layout()
    plt.savefig(f"image.png")
    # plt.show()

def main():
    start_time = perf_counter()
    ephemerides, time_vec, moon_mus, Jupter_mu = get_reference_ephemerides()
    system_states0 = [np.array(state, dtype=float) for state in ephemerides[0]]
    print(system_states0)

    jpl_step_seconds = JPL_STEP_DAYS * SECONDS_PER_DAY
    rk4_step_seconds = jpl_step_seconds // METHOD_SUBSTEPS_PER_JPL_STEP
    total_days = SIM_YEARS * DAYS_PER_YEAR
    total_seconds = int(total_days * SECONDS_PER_DAY)

    if jpl_step_seconds % rk4_step_seconds != 0:
        raise ValueError("RK4 step size must divide the JPL step size exactly.")

    if total_seconds % rk4_step_seconds != 0:
        raise ValueError("Total integration time must be an integer multiple of the RK4 step size.")

    sample_ratio = jpl_step_seconds // rk4_step_seconds
    rk4_ephemeris, rk4_time_vec = states_time_hist(
        system_states0,
        moon_mus,
        Jupter_mu,
        galilean_derivative,
        RK4,
        0,
        total_seconds,
        rk4_step_seconds,
    )
    symplectic_ephemeris, symplectic_time_vec = states_time_hist(
        system_states0,
        moon_mus,
        Jupter_mu,
        galilean_derivative,
        stormer_verlet,
        0,
        total_seconds,
        rk4_step_seconds,
    )

    aligned_final_days = int(total_days // JPL_STEP_DAYS) * JPL_STEP_DAYS
    aligned_final_seconds = aligned_final_days * SECONDS_PER_DAY
    aligned_final_step_index = aligned_final_seconds // rk4_step_seconds
    jpl_index = aligned_final_days // JPL_STEP_DAYS
    jpl_state_final = ephemerides[jpl_index]
    rk4_state_final = rk4_ephemeris[aligned_final_step_index]
    symplectic_state_final = symplectic_ephemeris[aligned_final_step_index]
    aligned_sample_count = min(len(ephemerides), len(rk4_ephemeris[::sample_ratio]), len(symplectic_ephemeris[::sample_ratio]))
    time_days = np.arange(aligned_sample_count) * JPL_STEP_DAYS
    jpl_sampled_ephemeris = ephemerides[:aligned_sample_count]
    rk4_sampled_ephemeris = rk4_ephemeris[::sample_ratio][:aligned_sample_count]
    symplectic_sampled_ephemeris = symplectic_ephemeris[::sample_ratio][:aligned_sample_count]
    rk4_errors = compute_position_error_histories(jpl_sampled_ephemeris, rk4_sampled_ephemeris)
    symplectic_errors = compute_position_error_histories(jpl_sampled_ephemeris, symplectic_sampled_ephemeris)

    elapsed_seconds = perf_counter() - start_time
    report = "\n\n".join([
        f"Requested comparison epoch: {SIM_YEARS} years ({total_days} days)\n"
        f"Displayed final-state epoch: {aligned_final_days} days\n"
        f"RK4 step size: {rk4_step_seconds / SECONDS_PER_DAY} days\n"
        f"Runtime before plot display: {elapsed_seconds:.2f} seconds",
        format_system_state("JPL final system state", jpl_state_final),
        format_system_state("RK4 final system state", rk4_state_final),
        format_system_state("Stormer-Verlet final system state", symplectic_state_final),
    ])
    print(report)

    plot_position_errors(time_days, rk4_errors, symplectic_errors)

if __name__ == "__main__":
    main()
