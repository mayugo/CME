import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14})

# Paràmetres
M = 1.0
m = 0.5
L = 0.6

f = 1.0
g = 9.81

t_max = 16
t_for = 6

# Força aplicada al carro
def F(t):
    # return 1.0  # Newtons
    return f if t <= t_for else 0.0  # Newtons

# Equacions no lineals
def deriv(t, y):
    x, x_dot, theta, theta_dot = y
    
    # Matriu de masses
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    
    # Calcul acceleracions resolent el sistema
    # ddx*(M + m) + m*L*ddtheta*cos(theta) - m*L*theta_dot^2*sin(theta) = F
    # L*ddtheta + ddx*cos(theta) + g*sin(theta) = 0
    # resolvem com a sistema lineal per ddx i ddtheta
    
    A = np.array([[M + m, m*L*cos_theta],
                  [cos_theta, L]])
    b = np.array([F(t) + m*L*theta_dot**2*sin_theta,
                  -g*sin_theta])
    
    ddx, ddtheta = np.linalg.solve(A, b)
    
    return [x_dot, ddx, theta_dot, ddtheta]

# Condicions inicials
y0 = [0, 0, np.pi/4, 0]  # pendol inclinat 45 graus

# Temps
t_span = (0, t_max)
t_eval = np.linspace(0, t_max, 1000)

sol = solve_ivp(
    deriv, 
    t_span, 
    y0, 
    t_eval=t_eval, 
    method='RK45',   # mètode Runge-Kutta clàssic d'ordre 4-5
    rtol=1e-8, 
    atol=1e-10
)
# Càlcul energies
x = sol.y[0]
x_dot = sol.y[1]
theta = sol.y[2]
theta_dot = sol.y[3]

T = 0.5*M*x_dot**2 + 0.5*m*(x_dot**2 + (L*theta_dot)**2 + 2*x_dot*L*theta_dot*np.cos(theta))
V = m*g*L*(1 - np.cos(theta))
E = T + V

# Crear figure 2x2
fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# Primer subplot: x (posició carro)
axs[0, 0].plot(sol.t, sol.y[0], label='x (posició carro)', color='blue')
axs[0, 0].set_xlabel('')  # treu el label de l'eix X
# axs[0, 0].set_xlabel(r'$t$ [s]')
axs[0, 0].set_ylabel(r'$x(t)$ [m]')
# axs[0, 0].legend()
axs[0, 0].grid(True)
# axs[0, 0].set_title('x(t)')

# Segon subplot: θ (angle pendol)
axs[0, 1].plot(sol.t, sol.y[2], label='θ (angle pendol)', color='orange')
axs[0, 1].set_xlabel(r'$t$ [s]')
axs[0, 1].set_ylabel(r'$\theta(t)$ [rad]')
# axs[0, 1].legend()
axs[0, 1].grid(True)
# axs[0, 1].set_title('θ(t)')

# Tercer subplot: x_dot (velocitat carro)
axs[1, 0].plot(sol.t, sol.y[1], label='ẋ (velocitat carro)', color='blue')
axs[1, 0].set_xlabel('')  # treu el label de l'eix 
# axs[1, 0].set_xlabel(r'$t$ [s]')
axs[1, 0].set_ylabel(r'$\dot{x}(t)$ [m/s]')
# axs[1, 0].legend()
axs[1, 0].grid(True)
# axs[1, 0].set_title('ẋ(t)')


# Quart subplot: θ_dot (velocitat angular pendol)
axs[1, 1].plot(sol.t, sol.y[3], label='θ̇ (velocitat pendol)', color='orange')
axs[1, 1].set_xlabel(r'$t$ [s]')
axs[1, 1].set_ylabel(r'$\dot{\theta}(t)$ [rad/s]')
# axs[1, 1].legend()
axs[1, 1].grid(True)
# axs[1, 1].set_title('θ̇(t)')

for ax in axs.flat:   # axs.flat recorre tots els subplots del 2x2
    ax.axvline(x=t_for, color='gray', linestyle='--', linewidth=3)

plt.tight_layout()
plt.savefig('CME_Exemple_5_4_SOL1.pdf', bbox_inches='tight', transparent=True)
plt.show()

# Plot energies
plt.figure(figsize=(10,5))
plt.plot(sol.t, T, label=r'Energia cinètica $T$')
plt.plot(sol.t, V, label=r'Energia potencial $V$')
plt.plot(sol.t, E, label=r'Energia total $T+V$', linestyle='--', linewidth=3)
plt.axvline(x=t_for, color='gray', linestyle='--', linewidth=3)

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'energia [J]')
plt.legend()
plt.xlim([0, t_max])
plt.grid(True)
plt.savefig('CME_Exemple_5_4_SOL2.pdf', bbox_inches='tight', transparent=True)
plt.show()