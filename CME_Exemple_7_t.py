import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14})

# --- DADES ---
m1 = 12
m2 = 10
mb = 16
Jb = 0.0533
d = 0.2
k1 = 3430
k2 = 2700

iopcio = 1

m1_ = [12, 13, 14, 15]  # kg
m1 = m1_[iopcio]

m2_ = [10, 9, 8, 11]  # kg
m2 = m2_[iopcio]

m3_ = [16, 17, 18, 19]  # kg
mb = m3_[iopcio]

AB_ = [0.2, 0.3, 0.4, 0.5]
AB = AB_[iopcio]


k1_ = [3430, 3530, 3630, 3730]  # N/m
k1 = k1_[iopcio]

k2_ = [2700, 2800, 2900, 3000]  # N/m
k2 = k2_[iopcio]

Jb = (1/12) * mb * (AB * AB)

c = 200  # N·s/m

F0_ = [35, 40, 45, 50]  # N
F0 = F0_[iopcio]

vtheta_ = [60, 45, 60, 45]  # degrees
vtheta = vtheta_[iopcio]

#%% ---  Determinació dels paràmetres del sistema equivalent ---

theta = 60 * np.pi / 180

J_eq =(m1*d**2*np.cos(theta)**2 +  
       m2*d**2*np.sin(theta)**2 + 
       Jb + mb*d**2/4)
c_eq = c*d**2*np.sin(theta)**2
k_eq = k1*d**2*np.cos(theta)**2 +  k2*d**2*np.sin(theta)**2

k_eq = 115.3
J_eq = 0.505333333


omega_n = np.sqrt(k_eq/J_eq)

omega_r = np.sqrt(k_eq/J_eq - 1/2*(c_eq/(J_eq))**2)
zeta = c_eq / (2 * np.sqrt(k_eq * J_eq))


t_final = 3

omega1 = omega_r
Q_eq = -3.5


print("\nParàmetres del sistema:")
print("-" * 50)
print(f"Inèrcia:        {J_eq:.3f} kg")
print(f"Rigidesa:       {k_eq:.1f} N/m")
print(f"Esmorteïment:   {c_eq:.2f} Ns/m")
print(f"zeta:           {zeta:.4f}")
print(f"omega_r:        {omega_r:.2f} rad/s")

print("-" * 50)

# =========================================================
# A) RESPOSTA FREQÜENCIAL ANALÍTICA
# =========================================================

omg = np.arange(0.01, 60, 0.01)

Omega = omg / omega_n

H_omg = 1 / np.sqrt((1 - Omega**2)**2 + (2*zeta*Omega)**2)
Theta_omg = np.arctan2(2*zeta*Omega, (1 - Omega**2))

fig1, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

axs[0].plot(omg, H_omg)
axs[0].set_ylabel(r"$H (j\omega)$ [-]")
axs[0].axvline(omega_n, color='tab:red', lw=3, ls='--', label=rf'\omega = {omega_n:.1f}$ rad/s')
axs[0].axvline(omega_r, color='orange', lw=3, ls='--', label=rf'\omega = {omega_r:.1f}$ rad/s')
axs[0].text(omega_n*1.05, min(H_omg), 
            r'$\Omega = 1$' + '\n' + rf'$\omega_n = {omega_n:.1f}\ \mathrm{{rad/s}}$',
            ha='left', va='bottom')
axs[0].text(omega_r/1.05, min(H_omg), 
            rf'$\Omega = {omega_r/omega_n:.2f}$' + '\n' +rf'$\;\omega_r = {omega_r:.1f}$ rad/s', 
            ha='right', va='bottom')

axs[0].text(omega_r/1.05, max(H_omg), rf' $\;H(\omega_r) = {max(H_omg):.2f}$', ha='right', va='bottom')
axs[0].set_xlim([0, 60])
axs[0].set_ylim([0, max(H_omg)*1.2])
axs[0].grid()

axs[1].plot(omg, Theta_omg * 180/np.pi) # Convertim a graus el desfasament
axs[1].axvline(omega_n, color='tab:red', lw=3, ls='--', label=rf'\omega = {omega_n:.1f}$ rad/s')
axs[1].axvline(omega_r, color='orange', lw=3, ls='--', label=rf'\omega = {omega_r:.1f}$ rad/s')
axs[1].set_xlabel(r"$\omega$ [rad/s]")
axs[1].set_ylabel(r"Desfase $\phi(j\omega)$ [$^\circ$]")
axs[1].set_xlim([0, 60])
axs[1].grid()

plt.tight_layout()
plt.savefig("CME_Exemple_7_t_SOL1.pdf", bbox_inches='tight'  )


#%% --- SISTEMA LTI ---
# G(s) = 1 / (J s^2 + c s + k)
sys = signal.lti([1.0], [J_eq, c_eq, k_eq])

# --- TEMPS ---
t = np.linspace(0, t_final, 3000)

# --- ENTRADES ---
u_sin = Q_eq * np.sin(omega1 * t)
u_step = Q_eq * np.ones_like(t)

# --- SIMULACIONS (lsim per tot) ---
_, x_sin, _  = signal.lsim(sys, U=u_sin, T=t)
_, x_step, _ = signal.lsim(sys, U=u_step, T=t)

# --- VALORS TEÒRICS ---
X  = Q_eq / np.sqrt((k_eq - J_eq*omega1**2)**2 + (c_eq*omega1)**2)
X0 = Q_eq / k_eq

print('X teòric: {:.4f} rad'.format(X))
print('X0 teòric: {:.4f} rad'.format(X0))


#%% --- PLOTS ---


fig1, ax1 = plt.subplots(figsize=(6, 5))

fact = 360/(2*np.pi)  # rad → graus

ax1.plot(t, x_sin*fact, label="resposta")

ax1.axhline(X*fact, color='orange', lw=1, ls='--',
            label=rf'$\Theta = {X*fact:.1f}$ graus')
ax1.axhline(-X*fact, color='orange', lw=1, ls='--')
ax1.axhline(0, color='k', lw=0.8)

ax1.set_xlabel(r"$t$ [s]")
ax1.set_ylabel(r"$\theta(t)$ [graus]")
ax1.set_xlim(0, t_final)

ax1.grid(alpha=0.3)
ax1.legend(loc='lower right')

ax1.set_title(
    rf"Resposta a $F(t)={F0}\sin({omega1:.1f} \cdot t)$ N",
    #y=-0.35,
    fontsize=16
)

plt.tight_layout()
plt.savefig('CME_Exemple_7_t_SOL2.pdf', bbox_inches='tight')



fig2, ax2 = plt.subplots(figsize=(6, 5))

ax2.plot(t, x_step*fact, label="resposta")

ax2.axhline(X0*fact, color='orange', lw=1, ls='--',
            label=rf'$\theta(t \rightarrow \infty) = {X0*fact:.1f}$ graus')
ax2.axhline(0, color='k', lw=0.8)

ax2.set_xlabel(r"$t$ [s]")
ax2.set_xlim(0, t_final)

ax2.grid(alpha=0.3)
ax2.legend(loc='right')

ax2.set_title(
    rf"Resposta a funció escaló $F(t)={F0}\,u(t)$ N",
    #y=-0.35,
    fontsize=16
)

plt.tight_layout()
plt.savefig('CME_Exemple_7_t_SOL3.pdf', bbox_inches='tight')

plt.show()