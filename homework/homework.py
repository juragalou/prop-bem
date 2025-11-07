from scipy.integrate import quad  # Module d'intégration "quad"
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt


import numpy as np
import naca16_509_m06_clcd as naca
import stdatm as sa


def compute_induction_factors(v, cst_pitch, hub_radius, Rtot, n_points, beta_deg, w, omega, B, c, beta_pitch=None):
    R = np.linspace(hub_radius, Rtot, n_points)
    a_list, A_list = [], []

    for r in R:
        # --- beta local ---
        if cst_pitch:
            beta = np.radians(beta_deg)
        else:
            if beta_pitch is None:
                raise ValueError("beta_pitch must be provided when cst_pitch is False")
            dbeta = np.deg2rad(beta_pitch) - np.deg2rad(beta_deg)
            p_ref_0 = 2*np.pi*0.75*Rtot*np.tan(np.radians(beta_deg))
            beta = np.arctan(p_ref_0 / (2*np.pi*r)) + dbeta

        a, A = 0.1, 0.0
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
            den_A = 2 * np.sin(2*phi)
            den_a = 2 * (1 - np.cos(2*phi))
            if abs(den_A) < 1e-8 or abs(den_a) < 1e-8:
                a, A = np.nan, np.nan
                break

            sigma = B * c / (2*np.pi*r)
            alpha = beta - phi
            Cl, Cd, _ = naca.naca16_509_m06(alpha)
            #if not np.isfinite(Cl) or not np.isfinite(Cd):
            #Cl = 2*np.pi* alpha
            #Cd = 0

            Cn = Cl*np.cos(phi) - Cd*np.sin(phi)
            Ct = Cl*np.sin(phi) + Cd*np.cos(phi)

            A_new = (sigma * Ct * (1 - A)) / den_A
            a_new = (sigma * Cn * (1 + a)) / den_a

            if not np.isfinite(a_new) or not np.isfinite(A_new):
                a, A = np.nan, np.nan
                break
            a_new = np.clip(a_new, -0.5, 1.0)
            A_new = np.clip(A_new, -0.5, 1.0)

            a_next = (1 - w)*a + w*a_new
            A_next = (1 - w)*A + w*A_new

            if (abs(a_next - a) < tol) and (abs(A_next - A) < tol):
                a, A = a_next, A_next
                converged = True
                break
            a, A = a_next, A_next

        if not converged and np.isfinite(a) and np.isfinite(A):
            a, A = np.nan, np.nan

        a_list.append(a)
        A_list.append(A)

    return R, np.array(a_list), np.array(A_list)



def Thrust(R, solutions_a, v, rho):
    """
    Calcule la poussée totale exercée par l’hélice ou la turbine.
    """
    #a_vec = np.asarray(solutions_a, float)
    #mask = np.isfinite(a_vec)
    #vals = 4*np.pi*R[mask]*rho*(v[mask]**2) * a_vec[mask] * (1.0 + a_vec[mask])

    
    vals =4 * np.pi * R * rho * v**2 * solutions_a * (1 + solutions_a)
    T_r = cumulative_trapezoid(vals, R, initial=0.0) 
    T_total = np.trapz(vals, R)
    return T_r,T_total

def torque(R, solutions_A, solutions_a,v, rho, omega):
    """
    Calcule le couple total exercé par l’hélice ou la turbine.
    """
    vals =4 * np.pi * R**3 * rho * v * omega * solutions_A * (1 + solutions_a)

    Q_r = cumulative_trapezoid(vals, R, initial=0.0) 
    Q_total = np.trapz(vals, R)
    return Q_r,Q_total

def mechanical_power(Q_total, omega):
    """
    Calcule la puissance mécanique fournie par l’hélice ou la turbine.
    """
    P_mech = Q_total * omega

    return P_mech


def K_t(cst_pitch, rho, Rtot, hub_radius, n_points, beta_deg, w, omega, B, c, beta_pitch=None):
    """
    Calcule le coefficient de poussée K_t.
    """
    n = omega / (2*np.pi)
    K_t_values = []
    for J in np.linspace(0,5, n_points):
        #v = J * (n * 2 * Rtot)
        v = J*n*2*Rtot


        R, a_factors, A_factors = compute_induction_factors(v, cst_pitch, hub_radius, Rtot, n_points, beta_deg, w, omega, B, c, beta_pitch=beta_pitch)
        
        _, T_total = Thrust(R, a_factors, v, rho)
        K_t = T_total / (rho * n**2 * (2 * Rtot)**4)
        K_t_values.append((K_t, J))
        
    return np.array(K_t_values)

def K_q(cst_pitch, rho, Rtot, hub_radius, n_points, beta_deg, w, omega, B, c, beta_pitch=None):
    """
    Calcule le coefficient de couple K_q.
    """
    n = omega / (2*np.pi)
    K_q_values = []
    for J in np.linspace(0,5, n_points):
        v = J*n*2*Rtot
        R, a_factors, A_factors = compute_induction_factors(v, cst_pitch, hub_radius, Rtot, n_points, beta_deg, w, omega, B, c, beta_pitch=beta_pitch)

        Q_r, Q_total = torque(R, A_factors, a_factors,v, rho, omega)
        K_q = Q_total / (rho * n**2 * (2 * Rtot)**5)
        K_q_values.append((K_q, J))
       
    return np.array(K_q_values)



def K_p(cst_pitch, rho, Rtot, hub_radius, n_points,
        beta_deg, w, omega, B, c, beta_pitch=None):
    n = omega / (2*np.pi)
    K_q_values = K_q(cst_pitch, rho, Rtot, hub_radius, n_points,
                     beta_deg, w, omega, B, c, beta_pitch=beta_pitch)
    J = K_q_values[:, 1]
    k_q = K_q_values[:, 0]

    # on garde seulement les points vraiment propulsifs
    valid = k_q > 0       # ou k_q >= 1e-5, à affiner
    k_q = k_q[valid]
    J   = J[valid]

    k_p = 2 * np.pi * np.abs(k_q)   # ICI pas de abs

    return np.column_stack([k_p, J])


def eta_curve(cst_pitch, rho, Rtot, hub_radius, n_points,
              beta_deg, w, omega, B, c, beta_pitch=None):
    n = omega / (2*np.pi)

    # --- Calcul des coefficients ---
    KT = K_t(cst_pitch, rho, Rtot, hub_radius, n_points,
             beta_deg, w, omega, B, c, beta_pitch=beta_pitch)
    KP = K_p(cst_pitch, rho, Rtot, hub_radius, n_points,
             beta_deg, w, omega, B, c, beta_pitch=beta_pitch)

    J_KT, kT = KT[:, 1], KT[:, 0]
    J_KP, kP = KP[:, 1], KP[:, 0]

    # --- Interpolation de kP sur la grille de J_KT ---
    kP_interp = np.interp(J_KT, J_KP, kP, left=np.nan, right=np.nan)

    # --- Calcul du rendement ---
    eta = np.full_like(kT, np.nan, dtype=float)
    eps_abs = 1e-4
    eps_rel = 1e-3 * np.nanmax(np.abs(kP_interp))
    eps = max(eps_abs, eps_rel)

    # 🔸 On garde seulement les points physiques
    valid = (kT > 0) & (kP_interp > eps)

    eta[valid] = J_KT[valid] * kT[valid] / kP_interp[valid]

    # --- Nettoyage et retour ---
    valid_idx = np.isfinite(eta)
    return np.column_stack([J_KT[valid_idx], eta[valid_idx]])


        