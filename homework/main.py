import numpy as np
import matplotlib.pyplot as plt
import homework as hw
import stdatm as sa

# === Choisir le cas à exécuter ===
EXERCICE = 2 # <-- mets 1 ou 2 ici

# === Paramètres communs ===
rho = 1.225        
w = 0.3     

if EXERCICE == 1:
    print("=== Exercice 1 : prop simple ===")
    Rtot = 0.5        
    u_0 = 10          
    omega = 2 * np.pi * 20
    c = 0.15          
    B = 2             
    beta_deg = 25     
    n_points = 100    
    hub_radius = 0.125 

    R, a_factors, A_factors = hw.compute_induction_factors(u_0=u_0, cst_pitch=True, hub_radius=hub_radius, Rtot=Rtot, n_points=n_points, beta_deg=beta_deg, w=w, omega=omega, B=B, c=c, beta_pitch=0)
    T_r, T_total = hw.Thrust(R, a_factors, u_0, rho)
    Q_r, Q_total = hw.torque(R, A_factors,a_factors,u_0, rho, omega)
    K_t_values = hw.K_t(n=omega/(2*np.pi), Vmax=30, cst_pitch=True, rho=rho, Rtot=Rtot, hub_radius=hub_radius, n_points=n_points, beta_deg=beta_deg, w=w, omega=omega, B=B, c=c)
    K_q_values = hw.K_q(n=omega/(2*np.pi), Vmax=30, cst_pitch=True, rho=rho, Rtot=Rtot, hub_radius=hub_radius, n_points=n_points, beta_deg=beta_deg, w=w, omega=omega, B=B, c=c)

if EXERCICE == 2:

    # === Paramètres communs ===
    print("=== Exercice 2 : Hamilton-Standard ===")
    Rtot = 3.4 / 2
    hub_radius = 0.45 / 2
    B = 4
    c = 0.25
    beta_deg = 15
    w = 0.3
    u_0 = 10
    omega = 2 * np.pi * 20
    rho = sa.stdatm(1000)[2]  # densité à 1000 m
    n_points = 100

    # --- Liste des beta_pitch à comparer ---
    beta_pitch_list = [15,20,25,30,35,40]

    # --- Création du graphe ---
    plt.figure(figsize=(8,5))

    for beta_pitch in beta_pitch_list:
        # Calcul pour chaque beta_pitch
        R, a_factors, A_factors = hw.compute_induction_factors(
            u_0=u_0, cst_pitch=False,
            hub_radius=hub_radius, Rtot=Rtot, n_points=n_points,
            beta_deg=beta_deg, w=w, omega=omega, B=B, c=c,
            beta_pitch=beta_pitch
        )

      
        # Calcul de kT(J)
        K_t_values = hw.K_t(
            n=omega/(2*np.pi), Vmax=163, cst_pitch=False, rho=rho,
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

if EXERCICE == 3:
    print("=== Exercice 3 ")

    
