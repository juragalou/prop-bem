from scipy.integrate import quad  # Module d'intégration "quad"
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt


import numpy as np
import naca16_509_m06_clcd as naca
import stdatm as sa





def compute_induction_factors( u_0, cst_pitch,hub_radius, Rtot, n_points, beta_deg, w, omega, B, c, beta_pitch= None):

    """
    Calcule les facteurs d’induction axiaux (a) et tangentiels (A)
    le long du rayon d’une hélice ou turbine à partir du modèle BEM simple ok ok .
    """
    R = np.linspace(hub_radius, Rtot, n_points)


    solutions_a = []
    solutions_A = []

    a = 0.3
    A = 0.01

    for r in R:
        if cst_pitch:
            beta = np.radians(beta_deg)
        else:
            if beta_pitch is None:
                raise ValueError("beta_pitch must be provided when cst_pitch is False") 
            dbeta = np.deg2rad(beta_pitch) - np.deg2rad(beta_deg)
            p_ref_0 = 2 * np.pi * 0.75 * Rtot * np.tan( np.radians(beta_deg))
            beta = np.arctan( p_ref_0 / (2 * np.pi * r) ) + dbeta

        a_new, A_new = 10, 10



        while np.abs(a - a_new) > 1e-2 or np.abs(A - A_new) > 1e-2 :

            phi = np.arctan((u_0 *(1+a))/ ((1-A)*(omega * r)))
            sigma = B * c / (2 * np.pi * r)


            alpha = beta - phi 

            Cl, Cd, flag = naca.naca16_509_m06(alpha)
            #Cl = 2 * np.pi * alpha  # Lift coefficient (linear approximation)
            #Cd = 0


            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)


            a_new = (sigma * Cn * (1 + a)) / (2 * (1 - np.cos(2*phi)))
            A_new = (sigma * Ct * (1 - A)) / (2 * np.sin(2*phi))


            a = (1-w)*a + w*a_new
            A = (1-w)*A + w*A_new

        solutions_a.append(a)
        solutions_A.append(A)

    return R, np.array(solutions_a), np.array(solutions_A)

def Thrust(R, solutions_a, v, rho):
    """
    Calcule la poussée totale exercée par l’hélice ou la turbine.
    """
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


def K_t(n , Vmax, cst_pitch, rho, Rtot, hub_radius, n_points, beta_deg, w, omega, B, c, beta_pitch=None):
    """
    Calcule le coefficient de poussée K_t.
    """
    K_t_values = []
    for v in np.linspace(0, Vmax, n_points):
        J = v / (n * 2 * Rtot)
        R, a_factors, A_factors = compute_induction_factors(v, cst_pitch, hub_radius, Rtot, n_points, beta_deg, w, omega, B, c, beta_pitch=beta_pitch)
        
        T_r, T_total = Thrust(R, a_factors, v, rho)
        K_t = T_total / (rho * n**2 * (2 * Rtot)**4)
        K_t_values.append((K_t, J))
        
    return np.array(K_t_values)

def K_q(n , Vmax, cst_pitch, rho, Rtot, hub_radius, n_points, beta_deg, w, omega, B, c, beta_pitch=None):
    """
    Calcule le coefficient de couple K_q.
    """
    K_q_values = []
    for v in np.linspace(0, Vmax, n_points):
        J = v / (n * 2 * Rtot)
        R, a_factors, A_factors = compute_induction_factors(v, cst_pitch, hub_radius, Rtot, n_points, beta_deg, w, omega, B, c, beta_pitch=beta_pitch)

        Q_r, Q_total = torque(R, A_factors, a_factors,v, rho, omega)
        K_q = Q_total / (rho * n**2 * (2 * Rtot)**5)
        K_q_values.append((K_q, J))
       
    return np.array(K_q_values)



def K_p(n , Vmax, cst_pitch, rho, Rtot, hub_radius, n_points, beta_deg, w, omega, B, c, beta_pitch=None):
    """
    Calcule le coefficient de puissance K_p.
    """
    K_p_values = []
    for v in np.linspace(0, Vmax, n_points):
        J = v / (n * 2 * Rtot)
        R, a_factors, A_factors = compute_induction_factors(v, cst_pitch, hub_radius, Rtot, n_points, beta_deg, w, omega, B, c, beta_pitch=beta_pitch)
        Q_r, Q_total = torque(R, A_factors, a_factors, v, rho, omega)

        P_mech = mechanical_power(Q_total, omega)
        K_p = P_mech / (rho * n**3 * (2 * Rtot)**5)
        K_p_values.append((K_p, J))

    return np.array(K_p_values)

