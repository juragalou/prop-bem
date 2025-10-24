from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

p = Path(__file__).parent / "Verification.txt"
data = np.genfromtxt(p, names=True, dtype=float, encoding="utf-8")  # whitespace-separated, handles "NaN"

J = data["J"]
kT = data["kT"]
mask = ~np.isnan(J) & ~np.isnan(kT)

plt.figure()
plt.plot(J[mask], kT[mask], "-o", markersize=3)
plt.xlabel("J")
plt.ylabel("kT")
plt.title("kT en fonction de J")
plt.grid(True)
plt.tight_layout()
plt.show()