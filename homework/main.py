import numpy as np
import matplotlib.pyplot as plt
import homework as hw
import stdatm as sa
import os
import pandas as pd

rho = 1.225
w   = 0.3

def exo1():
    print("=== Exercice 1 : prop simple ===")

    # --- Paramètres géométriques et opérationnels ---
    Rtot       = 0.5        # [m]
    hub_radius = 0.125      # [m]
    B          = 2
    c          = 0.15       # [m]
    beta_deg   = 25         # [deg]
    n_radial   = 100        # points BEM
    omega      = 2 * np.pi * 20   # [rad/s]

    # --- Calcul des courbes de performance avec la grosse fonction ---
    res = hw.compute_performance_curves(
        cst_pitch=True,
        rho=rho,
        Rtot=Rtot,
        hub_radius=hub_radius,
        n_radial=n_radial,
        beta_deg=beta_deg,
        w=w,
        omega=omega,
        B=B,
        c=c,
        beta_pitch=None,   # pas constant, donc None
        J_min=0.0,
        J_max=1.5,
        n_J=100
    )

    J   = res["J"]
    kT  = res["KT"]
    kQ  = res["KQ"]
    kP  = res["KP"]
    eta = res["eta"]

    # --- Masque pour enlever les NaN éventuels ---
    mask = np.isfinite(J) & np.isfinite(kT) & np.isfinite(kQ) & np.isfinite(kP)

    J_plot   = J[mask]
    kT_plot  = kT[mask]
    kQ_plot  = kQ[mask]
    kP_plot  = kP[mask]
    eta_plot = eta[mask]   # il peut rester des NaN si KP≤0, on peut remasquer si besoin

    # Si tu veux filtrer aussi eta :
    mask_eta = np.isfinite(J_plot) & np.isfinite(eta_plot)
    J_eta    = J_plot[mask_eta]
    eta_eta  = eta_plot[mask_eta]

    # --- Tracé des coefficients KT, KQ, KP ---
    plt.figure(figsize=(8,5))
    plt.plot(J_plot, kT_plot, marker='o', markersize=3, label=r"$k_T(J)$")
    plt.plot(J_plot, kQ_plot, marker='s', markersize=3, label=r"$k_Q(J)$")
    plt.plot(J_plot, kP_plot, marker='^', markersize=3, label=r"$k_P(J)$")

    plt.title("Exercice 1 — Coefficients $k_T(J)$, $k_Q(J)$ et $k_P(J)$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel("Coefficient")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # --- Tracé séparé du rendement η(J) ---
    plt.figure(figsize=(8,5))
    plt.plot(J_eta, eta_eta, marker='x', markersize=3, label=r"$\eta(J)$")

    plt.title("Exercice 1 — Rendement $\\eta(J)$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$\eta$")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()

def exo2():
    # === Paramètres communs ===
    print("=== Exercice 2 : Hamilton-Standard ===")
    Rtot       = 3.4 / 2
    hub_radius = 0.45 / 2
    B          = 4
    c          = 0.25
    beta_deg   = 15
    w          = 0.3
    rho        = sa.stdatm(1000)[2]  # densité à 1000 m
    n_points   = 50                  # points radiaux pour le BEM
    rpm_engine = 3000
    rpm_prop   = 0.477 * rpm_engine
    n          = rpm_prop / 60.0
    omega      = 2 * np.pi * n

    res60 = hw.compute_performance_curves(
    cst_pitch=False,
    rho=rho,
    Rtot=Rtot,
    hub_radius=hub_radius,
    n_radial=n_points,
    beta_deg=beta_deg,
    w=w,
    omega=omega,
    B=B,
    c=c,
    beta_pitch=60,
    J_min=0.1,
    J_max=5.0,
    n_J=50
)

 


    # --- Liste des beta_pitch à comparer ---
    beta_pitch_list = [10, 20, 30, 40, 50, 60]

    # --- Plage de J pour les courbes de perf ---
    J_min, J_max, n_J = 0.1, 5.0, 60

    # --- On calcule UNE FOIS toutes les courbes par beta_pitch ---
    curves_by_beta = {}
    for beta_pitch in beta_pitch_list:
        res = hw.compute_performance_curves(
            cst_pitch=False,
            rho=rho,
            Rtot=Rtot,
            hub_radius=hub_radius,
            n_radial=n_points,
            beta_deg=beta_deg,
            w=w,
            omega=omega,
            B=B,
            c=c,
            beta_pitch=beta_pitch,
            J_min=J_min,
            J_max=J_max,
            n_J=n_J
        )
        curves_by_beta[beta_pitch] = res

        hw.save_results(res, filename=f"beta{beta_pitch}.csv")


    # ===========================
    # 1) Graphe k_T(J)
    # ===========================
    plt.figure(figsize=(8, 5))
    for beta_pitch in beta_pitch_list:
        res = curves_by_beta[beta_pitch]
  

        J  = res["J"]
        kT = res["KT"]
        eta = res["eta"]

        # On garde uniquement les points valides et où η >= -0.11
        mask = np.isfinite(J) & np.isfinite(kT) & np.isfinite(eta) & (eta >= -0.11)
        J_plot  = J[mask]
        kT_plot = kT[mask]




        plt.plot(J_plot, kT_plot, marker='o', markersize=3,
                 label=fr"$\beta_{{pitch}}={beta_pitch}^\circ$")

    plt.title("Hamilton-Standard — Coefficient $k_T(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$k_T$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()

    # ===========================
    # 2) Graphe k_P(J)
    # ===========================
    plt.figure(figsize=(8, 5))
    for beta_pitch in beta_pitch_list:

        res = curves_by_beta[beta_pitch]
        J  = res["J"]
        kP = res["KP"]
        eta = res["eta"]

        mask = np.isfinite(J) & np.isfinite(kP) & np.isfinite(eta) & (eta >= -0.11)
        J_plot  = J[mask]
        kP_plot = kP[mask]

        plt.plot(J_plot, kP_plot, marker='o', markersize=3,
                label=fr"$\beta_{{pitch}}={beta_pitch}^\circ$")


    plt.title("Hamilton-Standard — Coefficient $k_P(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$k_P$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()

    # ===========================
    # 3) Graphe k_Q(J)
    # ===========================
    plt.figure(figsize=(8, 5))
    for beta_pitch in beta_pitch_list:
        res = curves_by_beta[beta_pitch]
        J  = res["J"]
        kQ = res["KQ"]
        eta = res["eta"]

        mask = np.isfinite(J) & np.isfinite(kQ) & np.isfinite(eta) & (eta >= -0.11)
        J_plot  = J[mask]
        kQ_plot = kQ[mask]

        plt.plot(J_plot, kQ_plot, marker='o', markersize=3,
                label=fr"$\beta_{{pitch}}={beta_pitch}^\circ$")


    plt.title("Hamilton-Standard — Coefficient $k_Q(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$k_Q$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()

    # ===========================
    # 4) Graphe η(J)
    # ===========================
    plt.figure(figsize=(8, 5))
    for beta_pitch in beta_pitch_list:
        res = curves_by_beta[beta_pitch]
        J   = res["J"]
        eta = res["eta"]

        mask = np.isfinite(J) & np.isfinite(eta) & (eta >= -0.11)
        J_plot   = J[mask]
        eta_plot = eta[mask]

        plt.plot(J_plot, eta_plot, marker='x', markersize=3,
                 label=fr"$\eta,\ \beta_{{pitch}}={beta_pitch}^\circ$")

    plt.title("Hamilton-Standard — Rendement $\\eta(J)$ pour différents $\\beta_{pitch}$")
    plt.xlabel(r"$J = V / (nD)$")
    plt.ylabel(r"$\eta$")
    plt.grid(True, alpha=0.4)
    plt.legend(title=r"$\beta_{pitch}$ (°)")
    plt.tight_layout()
    plt.show()





    


    


#if EXERCICE == 3:
def exo3(Z, P_engine_bhp):
    """
    Params:
      Z : altitude en ft
      P_engine : puissance moteur en bhp

    Trouve (u0, beta_pitch) tels que :
      1) T(u0,beta) = D(u0) + M g sin(theta)
      2) P_mech(u0,beta) = P_engine
    en s'appuyant sur homework.compute_induction_factors et homework.compute_forces.
    P_engine attendu en [W].
    """
    import numpy as np
    from scipy.optimize import fsolve
    import homework as hw
    import stdatm as sa

    # --- Paramètres prop / avion ---
    Rtot       = 3.4 / 2          # [m];s
    hub_radius = 0.45 / 2         # [m]
    gear_ratio = 0.477
    B          = 4
    c          = 0.25             # [m]
    beta_deg   = 15               # [deg] (pas de référence, pour le profil)
    w          = 0.3              # relaxation itérative BEM
    omega      = 2 * np.pi * (3000 * gear_ratio / 60)  # [rad/s]
    n_points   = 100
    z = Z * 0.3048  # altitude en m
    P_engine_W = P_engine_bhp * 745.7

    rho     = sa.stdatm(z)[2]
    M       = 8430 * 0.45359237   # [kg]
    A_wing  = 21.83               # [m^2]
    C_D0    = 0.0163
    g       = 9.81                # [m/s^2]
    e       = 0.8
    b       = 11.28               # [m]
    theta   = 0 #np.deg2rad(5)
    AR      = b**2 / A_wing
    K       = 1.0 / (np.pi * AR * e)

    def power(u_0, beta_pitch_deg):
        # Facteurs d'induction (attention: l'argument s'appelle v dans homework.py)
        R, a_factors, A_factors = hw.compute_induction_factors(
            v=u_0,
            cst_pitch=False,
            hub_radius=hub_radius,
            Rtot=Rtot,
            n_points=n_points,
            beta_deg=beta_deg,
            w=w,
            omega=omega,
            B=B,
            c=c,
            beta_pitch=beta_pitch_deg
        )

        # Poussée / Couple / Puissance à partir des facteurs
        # compute_forces retourne : T_r, Q_r, T_total, Q_total, P_mech
        _, _, T_total, Q_total, P_mech = hw.compute_forces(
            R, a_factors, A_factors, u_0, rho, omega
        )

        # Traînée avion (parasit + induite)
        D = 0.5 * rho * u_0**2 * A_wing * (
            C_D0 + K * ((M * g * np.cos(theta)) / (0.5 * rho * u_0**2 * A_wing))**2
        )

        return T_total, D, P_mech

    def equations(x):
        u_0, beta_pitch = x  # u_0 [m/s], beta_pitch [deg]
        T, D, P_mech = power(u_0, beta_pitch)
        eq1 = T - (D + M * g * np.sin(theta))
        eq2 = P_engine_W - P_mech
        return [eq1, eq2]

    # Point de départ raisonnable
    x0 = [200.0, 40.0]  # u0 ~ 100 m/s, beta_pitch ~ 20 deg

    sol, info, ier, msg = fsolve(equations, x0=x0, full_output=True)
    if ier != 1:
        print("⚠️ FSOLVE n'a pas convergé :", msg)

    u0_sol, beta_sol = sol

    # Sanity check / log
    T, D, P_mech = power(u0_sol, beta_sol)
    print(f"u0 = {u0_sol:.2f} m/s | beta = {beta_sol:.2f}°")
    print(f"T = {T:.0f} N | D + Mg sinθ = {(D + M*g*np.sin(theta)):.0f} N")
    print(f"P_mech = {P_mech/1e3:.1f} kW | P_engine = {P_engine_W/1e3:.1f} kW")

    return float(u0_sol), float(beta_sol)





def beta_pitch_optimal(
    Z, M,
    rpm_engine=5000,
    gear_ratio=0.477,
    Rtot=3.4/2,                  # m (diamètre 3.4 m → rayon 1.7 m)
    results_dir="results",       # CSV attendus: results/beta{10|20|...}.csv
    beta_pitch_list=(10,20,30,40,50,60),
    require_eta_positive=True,   # ignorer les points non propulsifs (eta<0)
    enforce_power=True           # respecter la contrainte P_mech <= P_engine
):
    """
    Choisit le beta_pitch qui maximise eta(J) au point de vol (z,M),
    en lisant les fichiers déjà calculés dans `results/`.

    Retourne un dict avec le pas choisi et les valeurs interpolées à J.
    """
    # --- 1) Vitesse de vol et J ---

    z = Z * 0.3048
    P, T_atm, rho = sa.stdatm(z)[:3]
    a = (rho * 287.05 * T_atm) ** 0.5              # vitesse du son [m/s]
    V = M * a                                      # vitesse avion [m/s]
    n = (gear_ratio * rpm_engine) / 60.0           # tr/s de l'hélice
    D = 2.0 * Rtot
    if n*D <= 1e-12:
        raise ValueError("n*D ≈ 0 : vérifie gear_ratio, rpm_engine et Rtot.")
    J = V / (n * D)



    best = {
        "beta_pitch_deg": None,
        "eta": -np.inf,
        "T": np.nan,
        "Q": np.nan,
        "P_mech": np.nan,
        "J": float(J),
        "V": float(V),
        "rho": float(rho),
        "n": float(n),
        "omega": float(2*np.pi*n),
    }

    # --- 2) Boucle sur les fichiers results/betaXX.csv ---
    for beta in beta_pitch_list:
        path = os.path.join(results_dir, f"beta{beta}.csv")
        if not os.path.exists(path):
            continue

        df = pd.read_csv(path)

        # Colonnes attendues : J, eta, (optionnellement T, Q, P_mech)
        if "J" not in df or "eta" not in df:
            continue

        J_vec   = df["J"].to_numpy(float)
        eta_vec = df["eta"].to_numpy(float)

        # filtre physique
        mask = np.isfinite(J_vec) & np.isfinite(eta_vec)
        if require_eta_positive:
            mask &= (eta_vec >= 0.0)
        if not np.any(mask):
            continue

        J_ok   = J_vec[mask]
        eta_ok = eta_vec[mask]

        # on ne prend que si J est couvert par la plage du CSV
        if not (J_ok.min() <= J <= J_ok.max()):
            continue

        # --- 3) Interpolation à J ---
        eta_J = float(np.interp(J, J_ok, eta_ok))

        # T, Q, P_mech si dispo
        T_J = Q_J = P_J = np.nan
        if "T" in df:
            T_J = float(np.interp(J, J_ok, df["T"].to_numpy(float)[mask]))
        if "Q" in df:
            Q_J = float(np.interp(J, J_ok, df["Q"].to_numpy(float)[mask]))
        if "P_mech" in df:
            P_J = float(np.interp(J, J_ok, df["P_mech"].to_numpy(float)[mask]))

        

    

        # --- 4) Maximisation de eta ---
        if np.isfinite(eta_J) and (eta_J > best["eta"]):
            best.update({
                "beta_pitch_deg": int(beta),
                "eta": eta_J,
                "T": T_J,
                "Q": Q_J,
                "P_mech": P_J,
            })

    if best["beta_pitch_deg"] is None:
        raise RuntimeError(
            f"Aucun fichier dans '{results_dir}' ne couvre J={J:.2f} "
            f"(ou tous dépassent P_engine si enforce_power=True)."
        )
    print("z =", Z)
    

    return best

# print(exo2())

# print(beta_pitch_optimal(
#     5000,  0.5 ,
#     rpm_engine=3000,
#     gear_ratio=0.477,
#     Rtot=3.4/2,                  # m (diamètre 3.4 m → rayon 1.7 m)
#     results_dir="results",       # CSV attendus: results/beta{10|20|...}.csv
#     beta_pitch_list=(10,20,30,40,50,60),
#     require_eta_positive=True,   # ignorer les points non propulsifs (eta<0)
#     enforce_power=True           # respecter la contrainte P_mech <= P_engine
# ))

print(exo3(5000, 1450))