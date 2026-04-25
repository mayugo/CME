import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14})

# ------------------------
# DADES
# ------------------------
omega1 = 2.33
omega2 = 7.88

phi1 = np.array([0.66, 0.7510])
phi2 = np.array([-0.3545, 0.935])

A = 1
t_f = 6.
t = np.linspace(0, t_f, 1000)

# ------------------------
# MOVIMENTS
# ------------------------
x11 = A * phi1[0] * np.cos(omega1 * t)
x21 = A * phi1[1] * np.cos(omega1 * t)

x12 = A * phi2[0] * np.cos(omega2 * t)
x22 = A * phi2[1] * np.cos(omega2 * t)


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
axs[1].plot(t, x11, label=r'$x_1(t)$')
axs[1].plot(t, x21, label=r'$x_2(t)$')
axs[1].axhline(0)
axs[1].grid()
axs[1].legend(loc='lower right')
axs[1].set_xlim([0,t_f])
axs[1].set_ylim([-1,1])
axs[1].set_xlabel('t [s]')

axs[1].text(0.5, -0.25, r"(b) Resposa en mode 1 a $\omega_1 =$" + str(omega1) + r" rad/s.",
            transform=axs[1].transAxes,
            ha='center')

# --- (3) Mode 2 temporal
axs[2].plot(t, x12, label=r'$x_1(t)$')
axs[2].plot(t, x22, label=r'$x_2(t)$')
axs[2].axhline(0)
axs[2].grid()
axs[2].legend(loc='lower right')
axs[2].set_xlim([0,t_f])
axs[2].set_ylim([-1,1])
axs[2].set_xlabel('t [s]')

axs[2].text(0.5, -0.25, r"(c) Resposa en mode 2 a $\omega_2 =$" + str(omega2) + r" rad/s.",
            transform=axs[2].transAxes,
            ha='center')

plt.tight_layout()
plt.savefig('CME_Exemple_6_6_SOL.pdf', bbox_inches='tight', transparent=True)
plt.show()