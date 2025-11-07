import numpy as np
import matplotlib.pyplot as plt
import homework as hw
import stdatm as sa
from scipy.optimize import fsolve 

# === Choisir le cas à exécuter ===
EXERCICE = 1 # <-- mets 1 ou 2 ici

# === Paramètres communs ===
rho = 1.225        
w = 0.3     

#if EXERCICE == 1:

def exo1():
    print("=== Exercice 1 : prop simple ===")
    Rtot = 0.5        
    omega = 2 * np.pi * 20
    c = 0.15          
    B = 2             
    beta_deg = 25     
    n_points = 100
    hub_radius = 0.125 

    # --- Calcul des facteurs et des coefficients ---

    K_t_values = hw.K_t(
        cst_pitch=True, rho=rho,
        Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
        beta_deg=beta_deg, w=w, omega=omega, B=B, c=c
    )

    K_q_values = hw.K_q(
         cst_pitch=True, rho=rho,
        Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
        beta_deg=beta_deg, w=w, omega=omega, B=B, c=c
    )

    K_p_values = hw.K_p(
         cst_pitch=True, rho=rho,
        Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
        beta_deg=beta_deg, w=w, omega=omega, B=B, c=c
    )

    eta_values = hw.eta_curve(True, rho, Rtot, hub_radius, n_points, beta_deg, w, omega, B, c)


    # --- Extraction et nettoyage des données ---

    J_T = K_t_values[:, 1]
    kT = K_t_values[:, 0]
    mask_T = ~np.isnan(kT)
    J_T = J_T[mask_T]
    kT = kT[mask_T]

    J_Q = K_q_values[:, 1]
    kQ = K_q_values[:, 0]
    mask_Q = ~np.isnan(kQ)
    J_Q = J_Q[mask_Q]
    kQ = kQ[mask_Q]

    J_P = K_p_values[:, 1]
    kP = K_p_values[:, 0]
    mask_P = ~np.isnan(kP)
    J_P = J_P[mask_P]
    kP = kP[mask_P]

    J_eta = eta_values[:,0]
    eta = eta_values[:,1]
    


    



    

    # --- Tracé des courbes ---
    plt.figure(figsize=(8,5))
    #plt.plot(J_T, kT, marker='o', markersize=3, label=r"$k_T(J)$")
    #plt.plot(J_Q, kQ, marker='s', markersize=3, label=r"$k_Q(J)$")
    #plt.plot(J_P,kP, marker='^', markersize=3, label=r"$k_P(J)$")
    plt.plot(J_eta, eta, marker='x', markersize=3, label=r"$\eta(J)$")

    # --- Mise en forme du graphe ---
    plt.title("Exercice 1 — Coefficients $k_T(J)$ et $k_Q(J)$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel("Coefficient")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()


#if EXERCICE == 2:
def exo2():

    # === Paramètres communs ===
    print("=== Exercice 2 : Hamilton-Standard ===")
    Rtot = 3.4 / 2
    hub_radius = 0.45 / 2
    B = 4
    c = 0.25
    beta_deg = 15
    w = 0.3
    rho = sa.stdatm(1000)[2]  # densité à 1000 m
    n_points = 50
    rpm_engine = 3000
    rpm_prop = 0.477 * rpm_engine
    n = rpm_prop / 60.0
    omega = 2*np.pi*n

    # --- Liste des beta_pitch à comparer ---
    beta_pitch_list = [10,20,30,40,50,60]
    """
    # --- Création du graphe ---
    plt.figure(figsize=(8,5))



    for beta_pitch in beta_pitch_list:
        
        # Calcul de kT(J)
        K_t_values = hw.K_t(
            cst_pitch=False, rho=rho,
            Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
            beta_deg=beta_deg, w=w, omega=omega, B=B, c=c,
            beta_pitch=beta_pitch
        )

        # Nettoyage des NaN
        J = K_t_values[:, 1]
        kT = K_t_values[:, 0]
        mask = np.isfinite(J) & np.isfinite(kT)
        J, kT = J[mask], kT[mask]

        # --- Tracé ---
        plt.plot(J, kT, marker='o', markersize=3, label=fr"$\beta_{{pitch}}={beta_pitch}^\circ$")

    # --- Mise en forme du graphe ---
    plt.title("Hamilton-Standard — Coefficient $k_T(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$k_T$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(8,5))

    for beta_pitch in beta_pitch_list:

      
        # Calcul de kT(J)
        K_p_values = hw.K_p(
            cst_pitch=False, rho=rho,
            Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
            beta_deg=beta_deg, w=w, omega=omega, B=B, c=c,
            beta_pitch=beta_pitch
        )

        # Nettoyage des NaN
        J = K_p_values[:, 1]
        kP = K_p_values[:, 0]
        mask = np.isfinite(J) & np.isfinite(kP) 
        J, kP = J[mask], kP[mask]

        # --- Tracé ---
        plt.plot(J, kP, marker='o', markersize=3, label=fr"$\beta_{{pitch}}={beta_pitch}^\circ$")

    # --- Mise en forme du graphe ---
    plt.title("Hamilton-Standard — Coefficient $k_P(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$k_P$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()
    
    
    plt.figure(figsize=(8,5))

    for beta_pitch in beta_pitch_list:


      
        # Calcul de kT(J)
        K_q_values = hw.K_q(
            cst_pitch=False, rho=rho,
            Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
            beta_deg=beta_deg, w=w, omega=omega, B=B, c=c,
            beta_pitch=beta_pitch
        )

        # Nettoyage des NaN
        J = K_q_values[:, 1]
        kQ = K_q_values[:, 0]
        mask = np.isfinite(J) & np.isfinite(kQ)
        J, kQ = J[mask], kQ[mask]


        # --- Tracé ---
        plt.plot(J, kQ, marker='o', markersize=3, label=fr"$\beta_{{pitch}}={beta_pitch}^\circ$")

    # --- Mise en forme du graphe ---
    plt.title("Hamilton-Standard — Coefficient $k_Q(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$k_Q$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()
    """
    # --- Tracé du rendement eta(J) pour chaque beta_pitch ---

# --- Tracé du rendement eta(J) pour chaque beta_pitch ---
    plt.figure(figsize=(8,5))
    for beta_pitch in beta_pitch_list:
        # Calcul de kT et kP
        eta_values = hw.eta_curve(cst_pitch=False, rho=rho,
            Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
            beta_deg=beta_deg, w=w, omega=omega, B=B, c=c,
            beta_pitch=beta_pitch
        )

        # Extraction des valeurs
        J = eta_values[:, 0]
        eta = eta_values[:, 1]



        plt.plot(J, eta, marker='x', markersize=3, label=fr"$\eta,\ \beta_{{pitch}}={beta_pitch}^\circ$")

    plt.title("Hamilton-Standard — Rendement $\\eta(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$\eta$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()


    


#if EXERCICE == 3:
def exo3(z ,P_engine):


    Rtot = 3.4 / 2
    hub_radius = 0.45 / 2
    B = 4
    c = 0.25
    beta_deg = 15
    w = 0.3
    omega = 2 * np.pi * 3000/60 
    n_points = 100
    rho = sa.stdatm(z)[2] 
    M = 8430 * 0.45359237    #[Kg]
    A_wing = 21.83 #[m^2]
    C_D0 = 0.0163 
    g = 9.81 #[m/s^2]
    e = 0.8
    b = 11.28 #[m]
    theta = np.deg2rad(5)
    AR = b**2/A_wing
    K = 1/(np.pi * AR * e)


    def power(u_0, beta_pitch):
        R, a_factors, A_factors = hw.compute_induction_factors(u_0=u_0, cst_pitch=False, hub_radius=hub_radius, Rtot=Rtot,n_points=n_points, beta_deg=beta_deg, w=w, omega=omega, B=B, c=c, beta_pitch= beta_pitch)
        D = 0.5 * rho  * u_0 **2 * A_wing * ( C_D0 + K* ((M*g*np.cos(theta))/(0.5*rho * u_0**2*A_wing))**2)
        _,T = hw.Thrust(R, a_factors, u_0, rho)
        _,Q = hw.torque(R, A_factors, a_factors,u_0,rho, omega )
        P_mech = hw.mechanical_power(Q, omega)

        return T, D, P_mech


    def equation(x):
        u_0, beta_pitch = x
        T, D, P_mech = power(u_0, beta_pitch)
        equ1 = T - (D + M * g * np.sin(theta))
        equ2 = P_engine - P_mech

        return[equ1,equ2]

    sol, infodict, ier, msg= fsolve(equation ,x0 = [100.0, 20.0], full_output = True)
    if ier != 1:
        print("FSOLVE n'a pas convergé:", msg)
    u0_sol, beta_sol = sol

    # sanity check
    T, D, P_mech = power(u0_sol, beta_sol)
    print(f"u0 = {u0_sol:.2f} m/s | beta = {beta_sol:.2f}°")
    print(f"T = {T:.0f} N | D + Mg sinθ = {(D + M*g*np.sin(theta)):.0f} N")
    print(f"P_mech = {P_mech/1e3:.1f} kW | P_engine = {P_engine/1e3:.1f} kW")

    return u0_sol , beta_sol 

print(exo2())