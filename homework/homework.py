from scipy.integrate import quad  # Module d'intégration "quad"
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt


import numpy as np
import naca16_509_m06_clcd as naca

rho=1.225        
w=0.3     
Rtot=0.5        
u_0=10          
omega=2 * np.pi * 19
c=0.15          
B=2             
beta_deg=25     
n_points=100    

def compute_induction_factors( u_0):

    """
    Calcule les facteurs d’induction axiaux (a) et tangentiels (A)
    le long du rayon d’une hélice ou turbine à partir du modèle BEM simple.
    """
    beta = np.deg2rad(beta_deg)
    R = np.linspace(0.125, Rtot, n_points) 


    solutions_a = []
    solutions_A = []

    a = 0.3
    A = 0.01

    for r in R:

        a_new, A_new = 10, 10



        while np.abs(a - a_new) > 1e-2 or np.abs(A - A_new) > 1e-2 :

            phi = np.arctan((u_0 *(1+a))/ ((1-A)*(omega * r)))
            sigma = B * c / (2 * np.pi * r)


            alpha = beta - phi 

            #Cl, Cd, flag = naca.naca16_509_m06(alpha)
            Cl = 2 * np.pi * alpha  # Lift coefficient (linear approximation)
            Cd = 0


            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)


            a_new = (sigma * Cn * (1 + a)) / (2 * max(1 - np.cos(2*phi), 1e-8))
            A_new = (sigma * Ct * (1 - A)) / (2 * max(np.sin(2*phi), 1e-8))


            a = (1-w)*a + w*a_new
            A = (1-w)*A + w*A_new

        solutions_a.append(a)
        solutions_A.append(A)

    return R, np.array(solutions_a), np.array(solutions_A)

def Thrust(R, solutions_a):
    """
    Calcule la poussée totale exercée par l’hélice ou la turbine.
    """
    vals =4 * np.pi * R * rho * u_0**2 * solutions_a * (1 + solutions_a)
    T_r = cumulative_trapezoid(vals, R, initial=0.0) 
    T_total = np.trapz(vals, R)
    return T_r,T_total

def torque(R, solutions_A, solutions_a):
    """
    Calcule le couple total exercé par l’hélice ou la turbine.
    """
    vals =4 * np.pi * R**3 * rho * u_0 * omega * solutions_A * (1 + solutions_a)
    Q_r = cumulative_trapezoid(vals, R, initial=0.0) 
    Q_total = np.trapz(vals, R)
    return Q_r,Q_total

def K_t(n , Vmax):
    """
    Calcule le coefficient de poussée K_t.
    """
    K_t_values = []
    for v in np.linspace(0, Vmax, 100):
        J = v / (n * 2 * Rtot)
        R, a_factors, A_factors = compute_induction_factors( u_0=v)
        T_r, T_total = Thrust(R, a_factors)
        K_t = T_total / (rho * n**2 * (2 * Rtot)**4)
        K_t_values.append((K_t, J))
    return np.array(K_t_values)
# TODO: tester une nouvelle méthode pour le calcul de Kt






if __name__ == "__main__":
    R, a_factors, A_factors = compute_induction_factors(u_0=u_0)
    T_r, T_total = Thrust(R, a_factors)
    Q_r, Q_total = torque(R, A_factors,a_factors)
    K_t_values = K_t(n=omega/(2*np.pi), Vmax=30)
    print("Rayons (m):", R)
    print("Facteurs d’induction axiaux (a):", a_factors)
    print("Facteurs d’induction tangentiels (A):", A_factors)
    print ("Poussée totale (N):", T_r)
    print("Poussée totale (N):", T_total)
    print("Couple total (N·m):", Q_r)
    print("Couple total (N·m):", Q_total)
    print("Coefficient de poussée K_t en fonction de J:", K_t_values)

    J = K_t_values[:, 1]
    kT = K_t_values[:, 0]
    mask = ~np.isnan(kT)
    J = J[mask]
    kT = kT[mask]
    plt.figure(figsize=(7,4))
    plt.plot(J, kT, marker='o', color='tab:blue', label=r"$k_T=f(J)$")
    plt.title("Coefficient de poussée $k_T$ en fonction du rapport d’avance $J$")
    plt.xlabel("J = V / (nD)")
    plt.ylabel(r"$k_T$")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()





