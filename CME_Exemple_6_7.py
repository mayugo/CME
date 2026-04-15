import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# DADES
# ------------------------
omega1 = 66.01
omega2 = 107.78

phi1 = np.array([-0.247, 0.969])
phi2 = np.array([0.178, 0.984])

A = 1
t_f = 0.5
t = np.linspace(0, t_f, 1000)

# ------------------------
# MOVIMENTS
# ------------------------
x1 = A * phi1[0] * np.cos(omega1 * t)
theta1 = A * phi1[1] * np.cos(omega1 * t)

x2 = A * phi2[0] * np.cos(omega2 * t)
theta2 = A * phi2[1] * np.cos(omega2 * t)


# ------------------------
# Figures
# ------------------------

fig, axs = plt.subplots(1, 3, figsize=(12,5))

# --- (1) VECTORS PROPIS
axs[0].quiver(0, 0, phi1[0], phi1[1],
              angles='xy', scale_units='xy', scale=1,
              label='Mode 1')

axs[0].quiver(0, 0, phi2[0], phi2[1],
              angles='xy', scale_units='xy', scale=1,
              color='red', label='Mode 2')

axs[0].set_xlim(-1,1)
axs[0].set_ylim(-1,1)
axs[0].set_aspect('equal')
axs[0].axhline(0)
axs[0].axvline(0)
axs[0].grid()
axs[0].legend(loc='lower right')

# títol abaix
axs[0].text(0.5, -0.25, '(a) Vectors propis.',
            transform=axs[0].transAxes,
            ha='center')

# --- (2) Mode 1 temporal
axs[1].plot(t, x1, label=r'$x(t)$')
axs[1].plot(t, theta1, label=r'$\theta(t)$')
axs[1].grid()
axs[1].legend(loc='lower right')
axs[1].set_xlim([0,t_f])
axs[1].set_ylim([-1,1])
axs[1].set_xlabel('t [s]')

axs[1].text(0.5, -0.25, '(b) Resposa en mode 1.',
            transform=axs[1].transAxes,
            ha='center')

# --- (3) Mode 2 temporal
axs[2].plot(t, x2, label=r'$x(t)$')
axs[2].plot(t, theta2, label=r'$\theta(t)$')
axs[2].grid()
axs[2].legend(loc='lower right')
axs[2].set_xlim([0,t_f])
axs[2].set_ylim([-1,1])
axs[2].set_xlabel('t [s]')

axs[2].text(0.5, -0.25, '(c) Resposa en mode 2.',
            transform=axs[2].transAxes,
            ha='center')

plt.tight_layout()
plt.savefig('CME_Exemple_6_7_SOL.pdf', bbox_inches='tight', transparent=True)
plt.show()