import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

plt.rcParams.update({
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14})

# -------------------------
# DADES DEL PROBLEMA
# -------------------------

m2 = 100        # kg, massa bloc

md = 22.5       # kg, massa disc
R =  0.4        # m, radi del dissc
G  = 60e9       # N/m2, mòdul de rigidesa
Ip = 250e-9     # m4, inèrcia polar
L = 1.3         # m, longitud

k2 = 100e3      # N/m, rigidesa molla

ct = 60         # Nms/rad, esmorteïment torsional

M0 = 280        # Nm, moment aplicat al disc (aplicat i tret)

# %%
# -------------------------
# SISTEMA EQUIVALENT
# -------------------------
J = 0.5 * md * R**2 + m2 * R**2
kt1 = G * Ip / L
kt2 = k2 * R**2
kt = kt1 + kt2

# %%
# -------------------------
# CONDICIONS INICIALS
# -------------------------
theta0 = M0 / kt
theta_dot0 = 0

# %%

m = J
c = ct
k = kt

x0 = theta0
x_dot0 = theta_dot0 

# %% Mètode 1: analític
# -------------------------
# RESPOSTA TEMPORAL expressions analítiques
# -------------------------

# Paràmetres dinàmics
cc = 2 * np.sqrt(k * m)
zeta = c / cc
omega_n = np.sqrt(k/m)

print(f"zeta =  {zeta:.4f}")
print(f"c_c  =  {cc:.2f} Nms/rad \n")

t = np.linspace(0, 2.0, 2000)

if zeta < 1:
    print("Sistema subesmorteït")
    
    alpha = c / (2 * m)
    omega_d = np.sqrt(k/m - c**2 / (4 * m**2))
    
    # Constants de la solució
    C1 = x0
    C2 = (x_dot0 + alpha * x0) / omega_d
    
    
    q = np.exp(-alpha * t) * (C1 * np.cos(omega_d * t) +
                                 C2 * np.sin(omega_d * t))
    
    envolvent = np.sqrt(C1**2 + C2**2) * np.exp(-alpha * t)
    
elif zeta > 1:
    print("Sistema sobreesmorteït")
    
    term = np.sqrt(zeta**2 - 1)
    
    A = ( x_dot0 / omega_n) + x0 * ( zeta + term)
    B = (-x_dot0 / omega_n) + x0 * (-zeta + term)
    
    q = (np.exp(-zeta * omega_n * t) / (2 * term)) * (
        A * np.exp(omega_n * term * t) +
        B * np.exp(-omega_n * term * t))
    
else:
    print("Sistema críticament esmorteït")
    
    q = np.exp(-omega_n * t) * (x0 + (x_dot0 + omega_n * x0) * t)
    
# --- GRÀFICA ---
plt.figure()
plt.plot(t, q, label=r'$\theta(t)$')
if zeta < 1 : 
    plt.plot(t, envolvent, '--',  color = 'orange', linewidth = 1)
    plt.plot(t, -envolvent, '--', color = 'orange', linewidth = 1)
plt.xlabel(r"$t$ [s]")
plt.ylabel(r'$\theta(t)$ [rad]')
# plt.title("Resposta analítica θ(t)")
plt.grid()
plt.savefig('CME_Exemple_4_1_SOL2.pdf', bbox_inches='tight', transparent=True)
plt.show()

# %% Mètode 2: numèric
# --- FUNCIÓ DE TRANSFERÈNCIA Θ(s) ---
# Θ(s) = (m*s*x0 + c*x0) / (m*s^2 + c*s + k)

num = [m*x0, c*x0]
den = [m, c, k]

system = signal.TransferFunction(num, den)

# --- RESPOSTA TEMPORAL
t = np.linspace(0, 2.0, 2000)

t, q = signal.impulse(system, T=t)

# NOTA: impulse() aquí dona la resposta natural equivalent
# perquè ja hem incorporat les condicions inicials al numerador

# --- GRÀFICA ---
plt.figure()
plt.plot(t, q)
plt.xlabel(r"$t$ [s]")
plt.ylabel(r'$\theta(t)$ [rad]')
plt.title("Resposta numèrica θ(t)")
plt.grid()
plt.savefig('CME_Exemple_4_1_SOL2.pdf', bbox_inches='tight', transparent=True)
plt.show()