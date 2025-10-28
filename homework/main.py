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

    # --- Calcul des facteurs et des coefficients ---
    R, a_factors, A_factors = hw.compute_induction_factors(
        u_0=u_0, cst_pitch=True, hub_radius=hub_radius, Rtot=Rtot,
        n_points=n_points, beta_deg=beta_deg, w=w, omega=omega, B=B, c=c, beta_pitch=0
    )

    T_r, T_total = hw.Thrust(R, a_factors, u_0, rho)
    Q_r, Q_total = hw.torque(R, A_factors, a_factors, u_0, rho, omega)

    K_t_values = hw.K_t(
        n=omega/(2*np.pi), Vmax=30, cst_pitch=True, rho=rho,
        Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
        beta_deg=beta_deg, w=w, omega=omega, B=B, c=c
    )

    K_q_values = hw.K_q(
        n=omega/(2*np.pi), Vmax=30, cst_pitch=True, rho=rho,
        Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
        beta_deg=beta_deg, w=w, omega=omega, B=B, c=c
    )

    K_p_values = hw.K_p(
        n=omega/(2*np.pi), Vmax=30, cst_pitch=True, rho=rho,
        Rtot=Rtot, hub_radius=hub_radius, n_points=n_points,
        beta_deg=beta_deg, w=w, omega=omega, B=B, c=c
    )


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

    #print(kP)
    eta = J_P * kT / kP

    neg_indices = np.where(eta < 0)[0]  # indices où kP < 0
    if len(neg_indices) > 0:
        idx = neg_indices[0]           # premier indice négatif
        eta = eta[:idx]
        J_P = J_P[:idx]


    

    # --- Tracé des courbes ---
    plt.figure(figsize=(8,5))
    #plt.plot(J_T, kT, marker='o', markersize=3, label=r"$k_T(J)$")
    #plt.plot(J_Q, kQ, marker='s', markersize=3, label=r"$k_Q(J)$")
    #plt.plot(J_P,kP, marker='^', markersize=3, label=r"$k_P(J)$")
    #plt.plot(J_P, eta, marker='x', markersize=3, label=r"$\eta(J)$")

    # --- Mise en forme du graphe ---
    plt.title("Exercice 1 — Coefficients $k_T(J)$ et $k_Q(J)$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel("Coefficient")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()


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


