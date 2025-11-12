import numpy as np
import matplotlib.pyplot as plt
import homework as hw
import stdatm as sa

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



def choose_pitch_and_performance(rpm_engine=3000, Mach=0.5, alt_ft=20000.0,
                                 gear_ratio=0.477, beta_pitch_list=(10,20,30,40,50,60),
                                 Rtot=3.4/2, hub_radius=0.45/2,
                                 B=4, c=0.25, beta_deg=15, w=0.3):
    # --- Atmosphère & vitesse ---
    z = alt_ft * 0.3048
    T, P, rho, a = sa.stdatm(z)   # si ta stdatm ne renvoie pas 'a', calcule-le: a = np.sqrt(1.4*287.05*T)
    V = Mach * a

    # --- Régime hélice ---
    n = gear_ratio * rpm_engine / 60.0
    omega = 2*np.pi*n
    D = 2*Rtot
    J_flight = V / (n*D)

    # --- Choix du pas: maximise eta(J_flight) ---
    best_eta = -np.inf
    best_beta = None
    for beta_pitch in beta_pitch_list:
        curve = hw.eta_curve(False, rho, Rtot, hub_radius, 80, beta_deg, w, omega, B, c, beta_pitch=beta_pitch)
        if curve.size == 0:
            continue
        Jc, etac = curve[:,0], curve[:,1]
        # on ne garde que la zone autour de J_flight
        if (J_flight < Jc.min()) or (J_flight > Jc.max()):
            continue
        eta_at_J = np.interp(J_flight, Jc, etac)
        if np.isfinite(eta_at_J) and eta_at_J > best_eta:
            best_eta = eta_at_J
            best_beta = beta_pitch

    if best_beta is None:
        raise RuntimeError(f"J_flight={J_flight:.2f} est hors des courbes; augmente la plage J ou la liste des pas.")

    # --- Perf à ce réglage (T, Q, P) ---
    R, a_f, A_f, phi, lam1, lam2 = hw.compute_induction_factors(
        v=V, cst_pitch=False, hub_radius=hub_radius, Rtot=Rtot, n_points=100,
        beta_deg=beta_deg, w=w, omega=omega, B=B, c=c, beta_pitch=best_beta
    )
    _, _, _, _, T_total, Q_total = hw.thrust_torque_lock(R, a_f, A_f, phi, lam1, lam2, v=V, rho=rho, B=B, c=c)
    P_mech = Q_total * omega
    eta_check = (T_total * V) / P_mech if P_mech > 1e-9 else np.nan

    return {
        "beta_pitch_deg": best_beta,
        "eta_at_J": best_eta,
        "J_flight": J_flight,
        "V": V, "rho": rho, "n": n, "omega": omega,
        "T": T_total, "Q": Q_total, "P_mech": P_mech, "eta_check": eta_check
    }


print(choose_pitch_and_performance(rpm_engine=3000, Mach=0.5, alt_ft=20000.0,
                                 gear_ratio=0.477, beta_pitch_list=(10,20,30,40,50,60),
                                 Rtot=3.4/2, hub_radius=0.45/2,
                                 B=4, c=0.25, beta_deg=15, w=0.3))