import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import StateSpace, lsim

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14})


# Paràmetres
m = 30e-6      # massa [kg]
c = 0.5e-3     # esmorteïment [Ns/m]
k = 1          # rigidesa [N/m]

x2I = 4e-6     # posició inicial massa m2
v2I = 4e-3     # velocitat inicial massa m2

t_final = 0.4  # temps simulacio [s]

# Matrius d'estat
A = np.array([[0,     0,     1,     0],
              [0,     0,     0,     1],
              [-2*k/m, k/m, -2*c/m, c/m],
              [k/m, -2*k/m, c/m, -2*c/m]])

B = np.array([[0],
              [0],
              [0],
              [1/m]])

C = np.array([[1, 0, 0, 0],
              [0, 1, 0, 0]])

D = np.array([[0],
              [0]])

# Sistema en espai d'estat
sys = StateSpace(A, B, C, D)

# Temps
t = np.arange(0, 0.4, 1e-4)

# Condicions inicials
x0 = np.array([0, x2I, 0, v2I])

# --- Resposta a condicions inicials ---
# Truquet: entrada zero
u0 = np.zeros_like(t)
T, Y_a, X_a = lsim(sys, U=u0, T=t, X0=x0)

# --- Força sinusoidal ---
f = 5e-6 * np.sin(210 * t)

# lsim espera entrada amb forma (N, inputs)
f = f.reshape(-1, 1)

T, Y_b, X_b = lsim(sys, U=f, T=t)

# --- Superposició ---
Z_1 = Y_a[:, 0] + Y_b[:, 0]
Z_2 = Y_a[:, 1] + Y_b[:, 1]

# --- Plots ---

# Força aplicada sinusoidal
plt.figure()
plt.plot(T, f, linewidth=2)
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'Força $f(t)$ [N]')
plt.title(r'Força $f(t)$')
plt.xlim([0,t_final])
plt.grid()

# Desplaçament del bloc x1
plt.figure()
plt.plot(T, Y_a[:, 0], '--', linewidth=2, label=r'$x_1$ inicial')
plt.plot(T, Y_b[:, 0], color = '#1f77b4', linewidth=2, label=r'$x_1$ força')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'posició $x_1(t)$ [m]')
plt.title('Resposta x_1(t)')
plt.legend()
plt.xlim([0,t_final])
plt.grid()

# Desplaçament del bloc x2
plt.figure()
plt.plot(T, Y_a[:, 1], '--', color = 'orange', linewidth=2, label=r'$x_2$ inicial')
plt.plot(T, Y_b[:, 1], color = 'orange', linewidth=2, label=r'$x_2$ força')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'posició $x_2(t)$ [m]')
plt.title('Resposta x_2(t)')
plt.legend()
plt.xlim([0,t_final])
plt.grid()

# Desplaçament total dels 2 blocs
plt.figure()
plt.plot(T, Z_1, linewidth=2, label=r'$x_1(t)$')
plt.plot(T, Z_2, color = 'orange', linewidth=2, label=r'$x_2(t)$')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'posicions $x_1(t)$ i $x_2(t)$ [m]')
plt.title('Resposta total')
plt.legend()
plt.xlim([0,t_final])
plt.grid()

plt.show()