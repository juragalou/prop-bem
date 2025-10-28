import numpy as np
import matplotlib.pyplot as plt
import homework as hw

"données globales"


rho=1.225        
w=0.3     
Rtot=0.5        
u_0=10          
omega=2 * np.pi * 20
c=0.15          
B=2             
beta_deg=25     
n_points=100    


R, a_factors, A_factors = hw.compute_induction_factors(u_0=u_0)
T_r, T_total = hw.Thrust(R, a_factors, u_0)
Q_r, Q_total = hw.torque(R, A_factors,a_factors,u_0)
K_t_values = hw.K_t(n=omega/(2*np.pi), Vmax=30)
K_q_values = hw.K_q(n=omega/(2*np.pi), Vmax=30)
print("Rayons (m):", R)
print("Facteurs d’induction axiaux (a):", a_factors)
print("Facteurs d’induction tangentiels (A):", A_factors)
print ("Poussée  (N):", T_r)
print("Poussée totale (N):", T_total)
print("Couple  (N·m):", Q_r)
print("Couple total (N·m):", Q_total)
print("Coefficient de poussée K_t en fonction de J:", K_t_values)
print("Coefficient de couple K_q en fonction de J:", K_q_values)


   # --- Extraction et nettoyage des données kT ---
J_T = K_t_values[:, 1]
kT = K_t_values[:, 0]
mask_T = ~np.isnan(kT)
J_T = J_T[mask_T]
kT = kT[mask_T]

# --- Extraction et nettoyage des données kQ ---
J_Q = K_q_values[:, 1]
kQ = K_q_values[:, 0]
mask_Q = ~np.isnan(kQ)
J_Q = J_Q[mask_Q]
kQ = kQ[mask_Q]

# --- Tracé sur une seule figure ---
plt.figure(figsize=(7,4))

plt.plot(J_T, kT, marker='o', color='tab:blue', label=r"$k_T(J)$")
plt.plot(J_Q, kQ, marker='s', color='tab:orange', label=r"$k_Q(J)$")

plt.title("Coefficients $k_T$ et $k_Q$ en fonction du rapport d’avance $J$")
plt.xlabel(r"$J = V / (nD)$")
plt.ylabel("Coefficient")
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()