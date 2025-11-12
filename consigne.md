## 🧩 **Section 5.1 – High speed performance tests**

### 1️⃣ Objectif

Tu dois utiliser **ton modèle BEM** pour **reproduire les mesures en vol** données dans le tableau  
(RPM = 3000, diverses altitudes, puissances moteur, vitesses mesurées, etc.).

---

### 2️⃣ Ce que tu dois faire

#### (a) Calculer la **vitesse de vol théorique** et le **pas optimal** (βₚᵢₜ꜀ʰ)

Pour chaque ligne du tableau :

- convertir la **puissance moteur** \( P_{\text{engine}} \) (bhp → W : multiplier par 745.7)  
- convertir l’**altitude** (ft → m)  
- puis, avec ton modèle (ta fonction `beta_pitch_optimal` par ex.) :
  - calculer la **vitesse de vol** \( V \) correspondante (c’est le \( J = \frac{V}{nD} \) qui équilibre poussée et traînée)
  - déterminer le **βₚᵢₜ꜀ʰ optimal** (celui qui maximise η tout en respectant la puissance moteur)

---

#### (b) Comparer tes résultats à la **vitesse réelle mesurée**

Tu feras un petit tableau comparatif :

| Altitude [ft] | P_engine [bhp] | V_mesurée [mph] | V_modèle [mph] | Différence [%] | β_pitch_optimal [°] |
|---------------|----------------|------------------|----------------|----------------|----------------------|
| ... | ... | ... | ... | ... | ... |

💬 Discute les écarts observés : résistance parasite, incertitudes sur la densité, simplifications BEM, etc.

---

#### (c) Analyser les tendances

Dans ton rapport :

- Comment **β_pitch** évolue-t-il avec l’altitude ?  
  → Il augmente (air plus rare → besoin d’un pas plus grand).  
- Comment évolue le **rendement propulsif η** ?  
  → Il diminue légèrement à haute altitude.  
- L’**advance ratio (J)** et le **Mach d’extrémité de pale** (tip Mach) :  
  → À haute vitesse, les extrémités atteignent Mach ≈ 0.85 → pertes de compressibilité.  
- Parle de la **répartition de l’angle d’attaque le long de la pale** (tu peux l’illustrer avec ton code).

---

## 🧩 **Section 5.2 – Climb performance tests**

### 1️⃣ Objectif

Cette fois, ce n’est plus du vol horizontal mais du **vol en montée**.  
On te donne le **taux de montée vertical** (*rate of climb*) en ft/min.

---

### 2️⃣ Ce que tu dois faire

#### (a) Calculer pour chaque ligne :

- la **vitesse totale du vol** :  
  \[
  u_0 = \sqrt{u_{0,x}^2 + u_{0,z}^2}
  \]
- l’**angle de montée** :  
  \[
  \theta = \arcsin\left(\frac{u_{0,z}}{u_0}\right)
  \]
- puis refaire ton calcul BEM :
  - pour obtenir la **vitesse horizontale** \( u_{0,x} \)
  - le **βₚᵢₜ꜀ʰ optimal**
  - la **poussée**, **couple** et **puissance mécanique**

---

#### (b) Calculer le **temps et la masse de carburant** consommée

Pour aller de 0 ft → 40000 ft :

- Prends les vitesses et puissances par segments (chaque ligne = un segment d’altitude)
- Calcule le temps de montée pour chaque intervalle
- Multiplie par le débit massique de carburant :  
  \[
  \dot{m}_f = \frac{P_{\text{engine}}}{\eta_{\text{prop}} \times LHV}
  \]
- Fais la somme pour obtenir la masse totale consommée.

---

#### (c) Discuter :

- Comment la vitesse et le rendement varient en montée ?  
- Pourquoi la poussée excédentaire diminue avec l’altitude ?  
- L’effet du régime de compresseur (“low blower” vs “high blower”) ?

---

## 🧩 **Section 5.3 – Take-off**

### 1️⃣ Contexte

Au décollage :

- \( \text{RPM} = 3000 \)
- \( P_{\text{engine}} = 1400\,\text{bhp} \)
- \( V_{\text{takeoff}} = 150\,\text{mph} = 67\,\text{m/s} \)
- \( \theta = 0° \)

---

### 2️⃣ Ce que tu dois calculer

#### (a) Coefficient de portance \( C_L \)

\[
C_L = \frac{2W}{\rho V^2 S}
\]

Compare-le à \( C_{L,\text{max}} \approx 1.5 \) (pour un avion de chasse à voilure laminaire).  
→ Cela te dira si la vitesse de décollage est réaliste.

---

#### (b) Poussée et βₚᵢₜ꜀ʰ à cette vitesse

Utilise ton modèle pour estimer :

- la poussée à 67 m/s  
- le pas optimal à cette vitesse  
- le rendement η  

---

#### (c) Accélération au décollage

\[
\frac{a}{g} = \frac{T - D}{W}
\]

→ Donne une estimation de l’accélération et du temps pour atteindre 150 mph.

---

## 💡 En résumé pratique

| Étape | Entrée | Sortie attendue | Outil |
|-------|---------|----------------|--------|
| **5.1** | (z, P_engine) | V_théorique, β_pitch_opt | `beta_pitch_optimal` |
| **5.2** | (z, P_engine, rate_of_climb) | θ, V, β_pitch_opt | Adaptation BEM + trig |
| **5.3** | (z=0, P_engine=1400 bhp, V=67 m/s) | T, β_pitch_opt, a/g | Ton modèle + formules simples |

---

