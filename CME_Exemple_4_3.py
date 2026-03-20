import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sympy as sp

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14})

# %%
# -------------------------
# Dades exercici
# -------------------------

R = 60e-3           # m, radi
m = 500e-3          # kg, massa
J = 100e-3/2*R**2   # kg·m2, inèrcia
c = 10              # Ns/m, esmorteiment
k1 = 1              # N/m, rigidesa
k2 = 2              # N/m, rigidesa

x0 = 0.05           # N, pretensió corretja
M = 10e-3           # M, moment torsor

t_jump  = 30        # s, temps on apliquem moment
t_final = 50        # s, temps final simulació

# %%
# -------------------------
# Model dinàmic
# -------------------------

def model(m,J,R,c,k1,k2,x0,M,t_jump,t_final):
    # -------------------------
    # Sistema equivalent
    # -------------------------
    m_eq = m + J/R**2
    c_eq = c
    k_eq = k1 + k2
    
    # Excitació
    # -------------------------
    x1_0 = 2*k1/(k1+k2)*x0
    x2_0 = 2*k2/(k1+k2)*x0

    Q_x0 = -k1*x1_0 + k2*x2_0   # força equivalent de les condicions inicials

    def torque(t):
        if t >= t_jump:
            return M/R
        else:
            return 0
    
    def Q(t):
        return torque(t) + Q_x0

    # Sistema d'EDO
    # -------------------------
    def system(t, y):
        x = y[0]
        v = y[1]
    
        dxdt = v
        dvdt = (Q(t) - c_eq*v - k_eq*x)/m_eq
    
        return [dxdt, dvdt]

    # Simulació
    # -------------------------
    t_span = (0,t_final)
    t_eval = np.linspace(0,t_final,4000)
            
    # Tram 1: abans del salt
    sol1 = solve_ivp(system, (0, t_jump), [0,0], t_eval=t_eval[t_eval <= t_jump])
    
    # condició inicial del tram 2 = final del tram 1
    y0_2 = sol1.y[:, -1]
    
    # Tram 2: després del salt
    sol2 = solve_ivp(system, (t_jump, t_final), y0_2, t_eval=t_eval[t_eval >= t_jump])
    
    # concatenar solució
    t_full = np.concatenate((sol1.t, sol2.t))
    x_full = np.concatenate((sol1.y[0], sol2.y[0]))
    v_full = np.concatenate((sol1.y[1], sol2.y[1]))
    
    Q_values = [Q(t) for t in t_full]
    
    return t_full,x_full,v_full, Q_values

# %%
# -------------------------
# Simulació i gràfic
# -------------------------

t_full, x_full,v_full, Q_values = model(m,J,R,c,k1,k2,x0,M,t_jump,t_final)
t_full2, x_full2,v_full2, Q_values = model(m*20,J*20,R,c,k1,k2,x0,M,t_jump,t_final)

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(8,8))

# Subplot 1: Q(t)
# -------------------------
axs[0].plot(t_full, Q_values, linestyle='--')
axs[0].set_ylabel(r'$Q(t) \, $[N]')
axs[0].grid(True)
axs[0].legend([r'$Q(t) = F_0 + M(t)/r \, $[N]'],loc='lower right')
axs[0].set_ylim(0, 0.3)

# Subplot 2: resposta x(t)
# -------------------------
axs[1].plot(t_full, x_full, label=r'Inèrcia baixa')
axs[1].plot(t_full2, x_full2, label=r'Inèrcia elevada')
axs[1].set_ylabel(r'$x\,$[m]')
axs[1].grid(True)
axs[1].legend()
axs[1].legend(loc='lower right')

# Subplot 3: resposta v(t)
# -------------------------
axs[2].plot(t_full, v_full, label=r'Inèrcia baixa')
axs[2].plot(t_full2, v_full2, label=r'Inèrcia elevada')
axs[2].set_xlabel(r'$t\,$[s]')
axs[2].set_ylabel(r'$\dot{x}\,$[m/s]')
axs[2].grid(True)

plt.xlim([20,t_final])
plt.tight_layout()
plt.show()

plt.tight_layout()
plt.savefig('CME_Exemple_4_3_SOL.pdf', bbox_inches='tight', transparent=True)
plt.show()
