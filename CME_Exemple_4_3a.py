import numpy as np
import matplotlib.pyplot as plt
import control as ctrl
import sympy as sp

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14})

# -------------------------
# Dades exercici
# -------------------------
R = 60e-3           # m, radi
m = 500e-3          # kg, massa
J = 100e-3/2*R**2   # kg·m2, inèrcia
c = 10              # Ns/m, esmorteiment
k1 = 1              # N/m, rigidesa
k2 = 2              # N/m, rigidesa

x0 = 0.075          # N, pretensió corretja
M =  5e-3           # M, moment torsor

tau  = 30           # s, temps on apliquem moment
t_final = 50        # s, temps final simulació

# -------------------------
# Sistema equivalent
# -------------------------
m_eq = m + J / R**2
c_eq = c
k_eq = k1 + k2

x1_0 = 2*k1/(k1+k2)*x0
x2_0 = 2*k2/(k1+k2)*x0

# -------------------------
# Model
# -------------------------
model = ctrl.TransferFunction([1], [m_eq, c_eq, k_eq])

# -------------------------
# Temps
# -------------------------
t = np.linspace(0, t_final, 1000)

# -------------------------
# Entrades en temps (en Python és complicat fer-ho en espai 's')
# -------------------------
# Component per pretensió (constant des de t=0)
F0 = (-k1*x1_0 + k2*x2_0)
Q0 = F0 * np.ones_like(t)

# Component retardada (esglaó a t=30)
QM = (M/R) * (t >= tau)

# Entrada total
Q = Q0 + QM

# -------------------------
# Resposta
# -------------------------
t_out1, y1 = ctrl.forced_response(model, t, Q)

# Només component retardada
t_out2, y2 = ctrl.forced_response(model, t, QM)

# -------------------------
# Gràfica
# -------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# -------------------------
# Gràfica de la Q
# -------------------------
ax[0].plot(t_out1, Q, linewidth=1.5, label=r'Pretensió + $M(t)$')
ax[0].plot(t_out2, QM, linewidth=1.5, label=r'només $M(t)$')

ax[0].set_xlabel('t [s]')
ax[0].set_ylabel('Q [N]')
ax[0].legend()
ax[0].set_xlim([20, 50])
ax[0].set_ylim([-0.01, 0.3])
ax[0].grid(True)

# -------------------------
# Gràfica de la Posició
# -------------------------
ax[1].plot(t_out1, y1, linewidth=1.5, label='Pretensió inicial + moment retardat')
ax[1].plot(t_out2, y2, linewidth=1.5, label='Moment amb retard')

ax[1].set_xlabel('t [s]')
ax[1].set_ylabel('x [m]')
# ax[1].legend()
ax[1].set_xlim([20, 50])
ax[1].grid(True)

# -------------------------
# Layout i guardat
# -------------------------
fig.tight_layout()
fig.savefig('CME_Exemple_4_3_SOL0.pdf', bbox_inches='tight', transparent=True)

plt.show()

# -------------------------
# Teorema valor inicial i final
# -------------------------
s = sp.symbols('s')

F_0  = (-k1*x1_0 + k2*x2_0)/(s**2*m_eq + c_eq*s + k_eq)*(1/s)
F_tau = (sp.exp(-tau*s)*M/R)/(s**2*m_eq + c_eq*s + k_eq)*(1/s)

F = F_0 + F_tau

initial_value_tau = sp.limit(s*F_tau, s, sp.oo)
final_value_tau =   sp.limit(s*F_tau, s, 0)

print("Aplicant només moment M(t)")
print("-"*50)
print("Valor inicial sense CI:", float(initial_value_tau),"m")
print("Valor final sense CI  :", round(float(final_value_tau),4),"m \n")

initial_value = sp.limit(s*F, s, sp.oo)
final_value =   sp.limit(s*F, s, 0)

print("Aplicant condicions inicials i moment M(t)")
print("-"*50)
print("Valor inicial:", float(initial_value),"m")
print("Valor final:  ", round(float(final_value),4),"m \n")