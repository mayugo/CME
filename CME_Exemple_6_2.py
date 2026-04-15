import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import StateSpace, lsim

# %%
# Paràmetres

m = 13      # kg
c = 3       # Ns/m
k = 2       # N/m

x0 = [0.04, 0]   # [posició inicial, velocitat inicial]

# Temps
t = np.arange(0, 50, 0.01)

# %%
# Paràmetres resposta dinàmica

alpha = c / (2*m)
omg_d = np.sqrt(k/m - (c**2)/(4*m**2))

# Representació 1
C1 = x0[0]
C2 = (x0[1] + alpha*x0[0]) / omg_d

x1 = np.exp(-alpha*t) * (C1*np.cos(omg_d*t) + C2*np.sin(omg_d*t))

# Representació 2
X = np.sqrt(C1**2 + C2**2)
psi_d = np.arctan2(x0[0]*omg_d, (x0[1] + alpha*x0[0]))

Xt = X * np.exp(-alpha*t)
x2 = X * np.exp(-alpha*t) * np.sin(omg_d*t + psi_d)

# Gràfica analítica
plt.figure()
plt.plot(t, x1, label='x1')
plt.plot(t, Xt, '--', label='envolvent')
plt.xlabel('Temps [s]')
plt.ylabel('x(t) [m]')
plt.grid()

# %%
# Espai d’estat

A = np.array([[0, 1],
              [-k/m, -c/m]])
B = np.array([[0],
              [1/m]])
C = np.array([[1, 0]])
D = np.array([[0]])

system = StateSpace(A, B, C, D)

# Resposta amb condicions inicials
u = np.zeros_like(t)  # entrada nul·la (resposta lliure)
t_out, y_out, x_out = lsim(system, U=u, T=t, X0=x0)

# Gràfica espai d’estat
plt.figure()
plt.plot(t_out, y_out)
plt.xlabel('Temps [s]')
plt.ylabel('x(t) [m]')
plt.title('Resposta amb espai d’estat')
plt.grid()

plt.show()