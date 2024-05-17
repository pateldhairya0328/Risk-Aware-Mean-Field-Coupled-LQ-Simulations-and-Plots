import csv
import numpy as np
import os
import sys
from tqdm import tqdm

def compute_feedback_and_affine_terms(T, Ax, Az, B, Qx, Qz, R, µ, Σ, γ, λ):
    """
    Compute the feedback matrices Kx and Kz, and compute the affine term f using
    dynamics matrices Ax = A, Az = A + C, B, cost matrices Qx = Q, Qz = P + Q,
    R, noise statistics µ, Σ, γ, and risk-parameter λ.
    """
    Kx = np.zeros(T)
    Kz = np.zeros(T)
    Sx = np.zeros(T + 1)
    Sz = np.zeros(T + 1)
    f  = np.zeros(T)
    g  = np.zeros(T + 1)

    Sx[T] = Qx + 4 * λ * Qx * Σ * Qx
    Sz[T] = Qz + 4 * λ * Qx * Σ * Qx
    g[T]  = 4 * λ * Qx * γ
    for t in range(T - 1, -1, -1):
        Sx[t] = - Ax * Sx[t + 1] * B * (R + B * Sx[t + 1] * B) ** -1 * B * Sx[t + 1] * Ax + Ax * Sx[t + 1] * Ax + Qx + 4 * λ * Qx * Σ * Qx
        Sz[t] = - Az * Sz[t + 1] * B * (R + B * Sz[t + 1] * B) ** -1 * B * Sz[t + 1] * Az + Az * Sz[t + 1] * Az + Qz + 4 * λ * Qx * Σ * Qx

        Kx[t] = - (R + B * Sx[t + 1] * B) ** -1 * B * Sx[t + 1] * Ax
        Kz[t] = - (R + B * Sz[t + 1] * B) ** -1 * B * Sz[t + 1] * Az

        g[t]  = (Az + B * Kz[t]) * (2 * Sz[t + 1] * µ + g[t + 1]) + 4 * λ * Qx * γ

        f[t]  = - (R + B * Sz[t + 1] * B) ** -1 * B * (Sz[t + 1] * µ + g[t + 1] / 2)

    return [Kx, Kz, f]

def compute_optimal_trajectory(k, T, x_0, w, A, B, C, Kx, Kz, f):
    """
    Compute the optimal trajectory given the optimal feedback matrices Kx and
    Kz, the affine feedback term f, the dynamics matrices A, B, and C, the
    initial states x_0, and the disturbance process w.
    """
    u  = np.zeros((T + 1, k))
    x  = np.zeros((T + 1, k))
    z  = np.zeros((T + 1))

    x[0] = x_0
    z[0] = np.mean(x[0])
    for t in range(T):
        for i in range(k):
            u[t, i]     = Kx[t] * x[t, i] + (Kz[t] - Kx[t]) * z[t] + f[t]
            x[t + 1, i] = A * x[t, i] + B * u[t, i] + C * z[t] + w[t, i]
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
    ## Set simulation parameters.
    print('Setting simulation parameters...', flush=True)

    k = 250 # Number of subsystems.
    T = 50 # Number of time steps.
    ts = np.arange(T + 1)
    N = 10000 # Number of simulations.

    A = 1.1
    B = 0.3
    C = 0.2
    P = 0.4
    Q = 0.8
    R = 1.2

    µ_x_0 = 10 # Mean of initial states.
    Σ_x_0 = 2 # Covariance of initial states.

    p = 0.25 # Bernoulli distribution parameter.
    shift = - p # Amount to shift Bernoulli distribution.
    scale = 10 # Amount to scale Bernoulli distribution after shift.
    µ = scale * (p + shift) # Mean of shifted and scaled Bernoulli distribution.
    Σ = p * (1 - p) * scale ** 2 # Covariance of shifted and scaled Bernoulli distribution.
    γ = p * (1 - p) * (1 - 2 * p) * Q * scale ** 3 # Weighted third moment of shifted and scaled Bernoulli distribution.

    # µ_normal = np.log(15 / 4) - np.log(3 - 2 * np.cbrt(3)) / 2 # Mean of normal distribution underlying log-normal distribution.
    # Σ_normal = np.log(-1 + np.cbrt(3) + 1 / np.cbrt(3)) # Variance of normal distribution underlying log-normal distribution.
    # shift = np.exp(µ_normal + Σ_normal / 2) # Amount to shift log-normal distribution.
    # µ = np.exp(µ_normal + Σ_normal / 2) - shift # Mean of shifted log-normal distribution.
    # Σ = (np.exp(Σ_normal) - 1) * np.exp(2 * µ_normal + Σ_normal) # Covariance of shifted log-normal distribution.
    # γ = (np.exp(Σ_normal) + 2) * np.square((np.exp(Σ_normal) - 1)) * np.exp(3 * µ_normal + 3 * Σ_normal / 2) * Q # Weighted third moment of shifted log-normal distribution.

    print('Done!', flush=True)

    ## Generate randomness.
    print('Generating simulation random data...', flush=True)

    np.random.seed(0) # Fixed seed for repeatability.

    # Generate fixed initial states.
    x_0 = np.random.normal(μ_x_0, Σ_x_0, k)

    # Generate a disturbance process for each simulation.
    vs = [None] * N
    ws = [None] * N
    for n in range(N):
        vs[n] = scale * (np.random.binomial(n=1, p=p, size=(T, k)) - p) # Generate Bernoulli noise.
        ws[n] = np.random.lognormal(mean=µ_normal, sigma=Σ_normal**0.5, size=(T, k)) - shift # Generate log-normal noise.

    print('Done!', flush=True)

    vs = np.array(vs).flatten()
    ws = np.array(ws).flatten()
    print(vs.mean(), vs.var(), np.power(vs - vs.mean(), 3).mean(), np.power(vs - vs.mean(), 4).mean())
    print(ws.mean(), ws.var(), np.power(ws - ws.mean(), 3).mean(), np.power(ws - ws.mean(), 4).mean())

    ## Run simulations.
    print('Running simulations...', flush=True)

    os.makedirs(fr'data', exist_ok=True)
    os.makedirs(fr'data/control_effort', exist_ok=True)
    os.makedirs(fr'data/control_effort/vs_lambda', exist_ok=True)
    os.makedirs(fr'data/control_effort/vs_time/avg', exist_ok=True)
    os.makedirs(fr'data/control_effort/vs_time/max', exist_ok=True)
    os.makedirs(fr'data/state_energy', exist_ok=True)
    os.makedirs(fr'data/state_energy/vs_lambda', exist_ok=True)
    os.makedirs(fr'data/state_energy/vs_time/avg', exist_ok=True)
    os.makedirs(fr'data/state_energy/vs_time/max', exist_ok=True)

    # Simulate for select values of λ to find performance of controller vs. time.
    λ     = [0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    λ_num = len(λ)

    avg_state_energy   = np.zeros((T + 1, N))
    max_state_energy   = np.zeros((T + 1, N))
    avg_control_effort = np.zeros((T + 1, N))
    max_control_effort = np.zeros((T + 1, N))

    for i in tqdm(range(λ_num), leave=False, desc='Generating statistics vs. t', disable=not sys.stdout.isatty()):
        [Kx, Kz, f] = compute_feedback_and_affine_terms(T, A, A + C, B, Q, P + Q, R, µ, Σ, γ, λ[i])

        for n in tqdm(range(N), leave=False, desc=f'λ = {λ[i]:0.3f}', disable=not sys.stdout.isatty()):
            [u, x, z] = compute_optimal_trajectory(k, T, x_0, ws[n], A, B, C, Kx, Kz, f)
            state_energies  = Q * x * x
            control_efforts = R * u * u
            avg_state_energy[:, n]   = np.mean(state_energies,  axis=1)
            max_state_energy[:, n]   = np.max(state_energies,   axis=1)
            avg_control_effort[:, n] = np.mean(control_efforts, axis=1)
            max_control_effort[:, n] = np.max(control_efforts,  axis=1)

        λ_str = str(λ[i]).replace('.', '_')
        write_data_to_csv(ts, 't', fr'data/state_energy/vs_time/avg/{λ_str}.csv',   avg_state_energy)
        write_data_to_csv(ts, 't', fr'data/state_energy/vs_time/max/{λ_str}.csv',   max_state_energy)
        write_data_to_csv(ts, 't', fr'data/control_effort/vs_time/avg/{λ_str}.csv', avg_control_effort)
        write_data_to_csv(ts, 't', fr'data/control_effort/vs_time/max/{λ_str}.csv', max_control_effort)

    λ_min_pow = - 5
    λ_max_pow =   3
    λ_num     = (λ_max_pow - λ_min_pow) * 20 + 1 # Gives 20 divisions of each power of 10 step.
    λ         = np.logspace(λ_min_pow, λ_max_pow, λ_num)

    avg_state_energy   = np.zeros((λ_num, N))
    max_state_energy   = np.zeros((λ_num, N))
    avg_control_effort = np.zeros((λ_num, N))
    max_control_effort = np.zeros((λ_num, N))

    for i in tqdm(range(λ_num), leave=False, desc='Generating time average of statistics vs. λ', disable=not sys.stdout.isatty()):
        [Kx, Kz, f] = compute_feedback_and_affine_terms(T, A, A + C, B, Q, P + Q, R, µ, Σ, γ, λ[i])

        for n in tqdm(range(N), leave=False, desc=f'λ = {λ[i]:0.3f}', disable=not sys.stdout.isatty()):
            [u, x, z] = compute_optimal_trajectory(k, T, x_0, ws[n], A, B, C, Kx, Kz, f)
            state_energies  = (Q * x * x)[1:]
            control_efforts = (R * u * u)[:T]
            avg_state_energy[i, n]   = np.mean(np.mean(state_energies,  axis=1))
            max_state_energy[i, n]   = np.mean(np.max(state_energies,   axis=1))
            avg_control_effort[i, n] = np.mean(np.mean(control_efforts, axis=1))
            max_control_effort[i, n] = np.mean(np.max(control_efforts,  axis=1))

    write_data_to_csv(λ, 'lambda', fr'data/state_energy/vs_lambda/avg.csv',   avg_state_energy)
    write_data_to_csv(λ, 'lambda', fr'data/state_energy/vs_lambda/max.csv',   max_state_energy)
    write_data_to_csv(λ, 'lambda', fr'data/control_effort/vs_lambda/avg.csv', avg_control_effort)
    write_data_to_csv(λ, 'lambda', fr'data/control_effort/vs_lambda/max.csv', max_control_effort)

    print('Done!', flush=True)