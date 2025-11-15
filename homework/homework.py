from scipy.integrate import quad  # Module d'intégration "quad"
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt


import numpy as np
import naca16_509_m06_clcd as naca
import stdatm as sa

def compute_induction_factors(
    v, cst_pitch, hub_radius, Rtot, n_points,
    beta_deg, w, omega, B, c, beta_pitch=None,
    a_init=None, A_init=None
):
    
    n = omega / (2 * np.pi)
    D = 2 * Rtot
    R = np.linspace(hub_radius, Rtot, n_points)
    a_list, A_list = [], []

    if a_init is None:
        a_init = np.full(n_points, 0.1)
    if A_init is None:
        A_init = np.full(n_points, 0.0)

    for j, r in enumerate(R):
        # --- beta local ---
        if cst_pitch:
            beta = np.radians(beta_deg)
        else:
            if beta_pitch is None:
                raise ValueError("beta_pitch must be provided when cst_pitch is False")
            dbeta = np.deg2rad(beta_pitch) - np.deg2rad(beta_deg)
            p_ref_0 = 2*np.pi*0.75*Rtot*np.tan(np.radians(beta_deg))
            beta = np.arctan(p_ref_0 / (2*np.pi*r)) + dbeta

        # conditions initiales (si NaN -> valeurs par défaut)
        a = a_init[j] if np.isfinite(a_init[j]) else 0.1
        A = A_init[j] if np.isfinite(A_init[j]) else 0.0

        max_iter = 200
        tol = 1e-3
        converged = False

        for _ in range(max_iter):
            Vax  = v * (1 + a)
            Vtan = omega * r * (1 - A)
            if Vtan <= 1e-8:
                a, A = np.nan, np.nan
                break

            phi = np.arctan2(Vax, Vtan)
            sigma = B * c / (2*np.pi*r)
            alpha = beta - phi
            Cl, Cd, _ = naca.naca16_509_m06(alpha)

            # polar simplifiée: CL = 2π α, CD = 0
            #Cl = 2*np.pi*alpha
            #Cd = 0.0

            Cn = Cl*np.cos(phi) - Cd*np.sin(phi)
            Ct = Cl*np.sin(phi) + Cd*np.cos(phi)

            den_a = 2 * (1 - np.cos(2*phi))
            den_A = 2 * np.sin(2*phi)

            den_a = np.sign(den_a) * max(abs(den_a), 1e-8)
            den_A = np.sign(den_A) * max(abs(den_A), 1e-8)

            a_new = (sigma * Cn * (1 + a)) / den_a
            A_new = (sigma * Ct * (1 - A)) / den_A

            a_next = (1 - w)*a + w*a_new
            A_next = (1 - w)*A + w*A_new

            if abs(a_next - a) < tol and abs(A_next - A) < tol:
                converged = True
                break

            a, A = a_next, A_next

        if not converged:
            a, A = 0,0
            


        a_list.append(a)
        A_list.append(A)

    return R, np.array(a_list), np.array(A_list)


def compute_forces(R, a_factors, A_factors, v, rho, omega):
    """
    Calcule en une seule fois :
      - la distribution radiale de poussée T_r(R)
      - la distribution radiale de couple Q_r(R)
      - la poussée totale T_total [N]
      - le couple total Q_total [N·m]
      - la puissance mécanique P_mech [W]
    à partir des facteurs d'induction a(r), A(r).
    """
    # Poussée élémentaire (même formule que dans Thrust)
    vals_T = 4 * np.pi * R * rho * v**2 * a_factors * (1 + a_factors)
    T_r = cumulative_trapezoid(vals_T, R, initial=0.0)
    T_total = np.trapz(vals_T, R)

    # Couple élémentaire (même formule que dans torque)
    vals_Q = 4 * np.pi * R**3 * rho * v * omega * A_factors * (1 + a_factors)
    Q_r = cumulative_trapezoid(vals_Q, R, initial=0.0)
    Q_total = np.trapz(vals_Q, R)

    # Puissance mécanique
    P_mech = Q_total * omega

    return T_r, Q_r, T_total, Q_total, P_mech


def compute_performance_curves(
    cst_pitch,
    rho,
    Rtot,
    hub_radius,
    n_radial,          # n_points pour le rayon
    beta_deg,
    w,
    omega,
    B,
    c,
    beta_pitch=None,
    J_min=0.0,
    J_max=5.0,
    n_J=50
):
    """
    Calcule en une seule fois, pour un range de J, les grandeurs suivantes :
        - T(J)      : poussée totale [N]
        - Q(J)      : couple total [N·m]
        - P_mech(J) : puissance mécanique [W]
        - K_T(J)    : coefficient de poussée
        - K_Q(J)    : coefficient de couple
        - K_P(J)    : coefficient de puissance
        - eta(J)    : rendement propulsif
    """

    n = omega / (2 * np.pi)
    D = 2 * Rtot

    # Discrétisation en J
    J_array = np.linspace(J_min, J_max, n_J)

    # Allocation des tableaux
    T_array   = np.full_like(J_array, np.nan, dtype=float)
    Q_array   = np.full_like(J_array, np.nan, dtype=float)
    P_array   = np.full_like(J_array, np.nan, dtype=float)
    KT_array  = np.full_like(J_array, np.nan, dtype=float)
    KQ_array  = np.full_like(J_array, np.nan, dtype=float)
    KP_array  = np.full_like(J_array, np.nan, dtype=float)
    eta_array = np.full_like(J_array, np.nan, dtype=float)

    eps_abs = 1e-4
    a_guess = None
    A_guess = None

    for i, J in enumerate(J_array):
        v = J * n * D

        R, a_factors, A_factors = compute_induction_factors(
            v, cst_pitch, hub_radius, Rtot, n_radial,
            beta_deg, w, omega, B, c,
            beta_pitch=beta_pitch,
            a_init=a_guess,
            A_init=A_guess
        )

        # Si tout est NaN, on ne réutilise pas ces valeurs pour la suite
        if (not np.isfinite(a_factors).any()) or (not np.isfinite(A_factors).any()):
            a_guess = None
            A_guess = None
            continue

        # Réutiliser la solution comme guess pour le prochain J
        a_guess = a_factors.copy()
        A_guess = A_factors.copy()

        # --- Forces + puissance ---
        _, _, T_total, Q_total, P_mech = compute_forces(
            R, a_factors, A_factors, v, rho, omega
        )

        T_array[i] = T_total
        Q_array[i] = Q_total
        P_array[i] = P_mech

        # --- Coefficients adimensionnels ---
        KT = T_total / (rho * n**2 * D**4)
        KQ = Q_total / (rho * n**2 * D**5)
        KP = 2 * np.pi * KQ

        KT_array[i] = KT
        KQ_array[i] = KQ
        KP_array[i] = KP

        # --- Rendement ---
        eps_rel = 1e-3 * max(abs(KP), 1e-12)
        eps = max(eps_abs, eps_rel)

        if (KP > eps) :
            eta_array[i] = J * KT / KP
        else:
            eta_array[i] = np.nan

    results = {
        "J"      : J_array,
        "T"      : T_array,
        "Q"      : Q_array,
        "P_mech" : P_array,
        "KT"     : KT_array,
        "KQ"     : KQ_array,
        "KP"     : KP_array,
        "eta"    : eta_array,
    }

    return results

import pandas as pd
import os

def save_results(results, filename="BEM_results.csv", folder="results"):
    """
    Sauvegarde les résultats BEM dans un fichier CSV.
    Si le dossier 'results/' n'existe pas, il est créé automatiquement.
    """
    # Créer le dossier si besoin
    os.makedirs(folder, exist_ok=True)
    filepath = os.path.join(folder, filename)

    # Transformer le dictionnaire en DataFrame
    df = pd.DataFrame({
        "J": results["J"],
        "KT": results["KT"],
        "KQ": results["KQ"],
        "KP": results["KP"],
        "eta": results["eta"],
        "T": results["T"],
        "Q": results["Q"],
        "P_mech": results["P_mech"]
    })

    # Sauvegarde
    df.to_csv(filepath, index=False)
    print(f"✅ Résultats sauvegardés dans : {filepath}")



    
