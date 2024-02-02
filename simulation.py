import csv
import numpy as np
import os

def compute_feedback_and_affine_terms(Ax, Az, B, Qx, Qz, R, µ, Σ, γ, λ):
    """
    Compute the feedback matrices Kx and Kz, and compute the affine term f using
    dynamics matrices Ax = A, Az = A + C, B, cost matrices Qx = Q, Qz = P + Q,
    R, noise statistics µ, Σ, γ, and risk-parameter λ.
    """
    T = len(Ax)

    Kx = np.zeros(T)
    Kz = np.zeros(T)
    Sx = np.zeros(T + 1)
    Sz = np.zeros(T + 1)
    f  = np.zeros(T)
    g  = np.zeros(T + 1)

    Sx[T] = Qx[T] + 4 * λ * Qx[T] * Σ * Qx[T]
    Sz[T] = Qz[T] + 4 * λ * Qx[T] * Σ * Qx[T]
    g[T]  = 4 * λ * Qx[T] * γ[T]
    for t in range(T - 1, -1, -1):
        Sx[t] = - Ax[t] * Sx[t + 1] * B[t] * (R[t] + B[t] * Sx[t + 1] * B[t]) ** -1 * B[t] * Sx[t + 1] * Ax[t] + Ax[t] * Sx[t + 1] * Ax[t] + Qx[T] + 4 * λ * Qx[T] * Σ * Qx[T]
        Sz[t] = - Az[t] * Sz[t + 1] * B[t] * (R[t] + B[t] * Sz[t + 1] * B[t]) ** -1 * B[t] * Sz[t + 1] * Az[t] + Az[t] * Sz[t + 1] * Az[t] + Qz[T] + 4 * λ * Qx[T] * Σ * Qx[T]

        Kx[t] = - (R[t] + B[t] * Sx[t + 1] * B[t]) ** -1 * B[t] * Sx[t + 1] *  A[t]
        Kz[t] = - (R[t] + B[t] * Sz[t + 1] * B[t]) ** -1 * B[t] * Sz[t + 1] * Az[t]

        g[t]  = (Az[t] + B[t] * Kz[t]) * (2 * Sz[t + 1] * µ + g[t + 1]) + 4 * λ * Qx[t] * γ[t]

        f[t]  = - (R[t] + B[t] * Sz[t + 1] * B[t]) ** -1 * B[t] * (Sz[t + 1] * µ + g[t + 1] / 2)

    return [Kx, Kz, f]

def compute_optimal_trajectory(x_0, w, A, B, C, Kx, Kz, f):
    """
    Compute the optimal trajectory given the optimal feedback matrices Kx and
    Kz, the affine feedback term f, the dynamics matrices A, B, and C, the
    initial states x_0, and the disturbance process w.
    """
    k = len(x_0)
    T = len(w)

    u  = np.zeros((T + 1, k))
    x  = np.zeros((T + 1, k))
    z  = np.zeros((T + 1))

    x[0] = x_0
    z[0] = np.mean(x[0])
    for t in range(T):
        for i in range(k):
            u[t, i]     = Kx[t] * x[t, i] + (Kz[t] - Kx[t]) * z[t] + f[t]
            x[t + 1, i] = A[t] * x[t, i] + B[t] * u[t, i] + C[t] * z[t] + w[t, i]
        z[t + 1] = np.mean(x[t + 1])

    return [u, x, z]

def write_data_to_csv(independent_var, independent_var_name, file_name, data):
    """
    Save the quantiles, mean, and standard deviation of the data to a csv file.
    """
    header_row = [independent_var_name] + [f'{x / 2:0.1f}'.rstrip('0').rstrip('.') for x in range(201)] + ['mean', 'stddev']
    table_rows = [independent_var] + [None] * 201 + [np.mean(data, axis=1), np.std(data, axis=1)]

    for i in range(201):
        table_rows[i + 1] = np.quantile(data, i / 200.0, axis=1)

    with open(file_name, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(header_row)
        writer.writerows(zip(*table_rows))

if __name__ == '__main__':
    ## Parameters.
    print('Setting parameters and generating data...')

    k = 250 # Number of subsystems.
    T = 50 # Number of time steps.
    ts = np.arange(T + 1)
    N = 10000 # Number of simulations.

    A = [1.1] * np.ones(T)
    B = [0.3] * np.ones(T)
    C = [0.2] * np.ones(T)
    P = [0.4] * np.ones(T + 1)
    Q = [0.8] * np.ones(T + 1)
    R = [1.2] * np.ones(T + 1)

    p = 0.25 # Bernoulli distribution parameter.
    shift = - p # Amount to shift Bernoulli distribution.
    scale = 10 # Amount to scale Bernoulli distribution after shift.
    µ = scale * (p + shift) # Mean of disturbance.
    Σ = p * (1 - p) * scale ** 2 # Covariance of disturbance.
    γ = p * (1 - p) * (1 - 2 * p) * Q * scale **3 # Weighted third moment of distribution.

    µ_x_0 = 10 # Mean of initial states.
    Σ_x_0 = 2 # Covariance of initial states.

    ## Generate randomness.
    np.random.seed(0) # Fixed seed for repeatability.

    # Generate a disturbance process for each simulation.
    ws = [None] * N
    for n in range(N):
        ws[n] = scale * (np.random.binomial(n=1, p=p, size=(T, k)) + shift)

    # Generate fixed initial states.
    x_0 = np.random.normal(μ_x_0, Σ_x_0, k)

    print('Done!')

    ## Run simulations.
    print('Running simulations...')

    # Simulate for select values of λ to find performance of controller vs. time.
    print('\tRunning simulations to generate statistics vs. t...')

    λ     = [0, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1]
    λ_num = len(λ)

    avg_state_energy   = np.zeros((T + 1, N))
    max_state_energy   = np.zeros((T + 1, N))
    avg_control_effort = np.zeros((T + 1, N))
    max_control_effort = np.zeros((T + 1, N))

    for i in range(λ_num):
        [Kx, Kz, f] = compute_feedback_and_affine_terms(A, A + C, B, Q, P + Q, R, µ, Σ, γ, λ[i])
        λ_str = str(λ[i]).replace('.', '_')

        for n in range(N):
            print(f'\t\tλ = {λ[i]:0.3f} ({i + 1} of {λ_num}): {100 * n / N:.2f}%', end='\r', flush=True)

            [u, x, z] = compute_optimal_trajectory(x_0, ws[n], A, B, C, Kx, Kz, f)
            state_energies  = x * Q[:, None] * x
            control_efforts = u * R[:, None] * u
            avg_state_energy[:, n]   = np.mean(state_energies,  axis=1)
            max_state_energy[:, n]   = np.max(state_energies,   axis=1)
            avg_control_effort[:, n] = np.mean(control_efforts, axis=1)
            max_control_effort[:, n] = np.max(control_efforts,  axis=1)

        os.makedirs(os.path.dirname(fr'Simulation Data/lambda_{λ_str}/'), exist_ok=True)
        write_data_to_csv(ts, 't', fr'Simulation Data/lambda_{λ_str}/time_varying_avg_state_energy.csv',   avg_state_energy)
        write_data_to_csv(ts, 't', fr'Simulation Data/lambda_{λ_str}/time_varying_max_state_energy.csv',   max_state_energy)
        write_data_to_csv(ts, 't', fr'Simulation Data/lambda_{λ_str}/time_varying_avg_control_effort.csv', avg_control_effort)
        write_data_to_csv(ts, 't', fr'Simulation Data/lambda_{λ_str}/time_varying_max_control_effort.csv', max_control_effort)

        print(f'\t\tλ = {λ[i]:0.3f} ({i + 1} of {λ_num}): 100.00%')

    print('\tDone!')

    # Simulate for range of λ to find performance of controller vs. λ.
    print('\tRunning simulations to generate time average of statistics vs. λ...')

    λ_min_pow = - 5
    λ_max_pow =   3
    λ_num     = (λ_max_pow - λ_min_pow + 1) * 20 # Gives 20 divisions of each power of 10 step.
    λ         = np.logspace(λ_min_pow, λ_max_pow, λ_num)

    avg_state_energy   = np.zeros((λ_num, N))
    max_state_energy   = np.zeros((λ_num, N))
    avg_control_effort = np.zeros((λ_num, N))
    max_control_effort = np.zeros((λ_num, N))

    for i in range(λ_num):
        [Kx, Kz, f] = compute_feedback_and_affine_terms(A, A + C, B, Q, P + Q, R, µ, Σ, γ, λ[i])

        for n in range(N):
            print(f'\t\tλ = {λ[i]:0.6f} ({i + 1} of {λ_num}): {100 * n / N:.2f}%', end='\r', flush=True)

            [u, x, z] = compute_optimal_trajectory(x_0, ws[n], A, B, C, Kx, Kz, f)
            state_energies  = (x * Q[:, None] * x)[1:]
            control_efforts = (u * R[:, None] * u)[:T]
            avg_state_energy[i, n]   = np.mean(np.mean(state_energies,  axis=1))
            max_state_energy[i, n]   = np.mean(np.max(state_energies,   axis=1))
            avg_control_effort[i, n] = np.mean(np.mean(control_efforts, axis=1))
            max_control_effort[i, n] = np.mean(np.max(control_efforts,  axis=1))

        print(f'\t\tλ = {λ[i]:0.6f} ({i + 1} of {λ_num}): 100.00%')

    write_data_to_csv(λ, 'lambda', 'Simulation Data/lambda_varying_avg_state_energy.csv',   avg_state_energy)
    write_data_to_csv(λ, 'lambda', 'Simulation Data/lambda_varying_max_state_energy.csv',   max_state_energy)
    write_data_to_csv(λ, 'lambda', 'Simulation Data/lambda_varying_avg_control_effort.csv', avg_control_effort)
    write_data_to_csv(λ, 'lambda', 'Simulation Data/lambda_varying_max_control_effort.csv', max_control_effort)

    print('\tDone!')

    print('Done!')