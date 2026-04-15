import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------------------------
# 1️⃣ Paràmetres
# ---------------------------
m1, m2 = 1.0, 2.0
L1, L2 = 1.0, 2.0
g = 9.81

def M_ext(t):
    return 0.0

# ---------------------------
# 2️⃣ Equacions del doble pèndol
# ---------------------------
def double_pendulum_rhs(t, y):
    th1, th2, w1, w2 = y
    c = np.cos(th1 - th2)
    s = np.sin(th1 - th2)
    
    M11 = 1/3*m1*L1**2 + m2*L1**2
    M12 = 0.5*m2*L1*L2*c
    M21 = M12
    M22 = 1/3*m2*L2**2
    
    C1 = -0.5*m2*L1*L2*s*w2**2
    C2 = 0.5*m2*L1*L2*s*w1**2
    
    G1 = (0.5*m1 + m2)*g*L1*np.sin(th1)
    G2 = 0.5*m2*g*L2*np.sin(th2)
    
    M = np.array([[M11, M12],
                  [M21, M22]])
    rhs = np.array([M_ext(t) - C1 - G1,
                    -C2 - G2])
    
    acc = np.linalg.solve(M, rhs)
    return [w1, w2, acc[0], acc[1]]

# ---------------------------
# 3️⃣ Condicions inicials
# ---------------------------
y0 = [np.pi/6, np.pi/12, 0.0, 0.0]
t_span = (0, 10)
t_eval = np.linspace(0, 10, 500)

sol = solve_ivp(double_pendulum_rhs, t_span, y0, t_eval=t_eval, method='RK45')

# ---------------------------
# 4️⃣ Gràfic dels angles
# ---------------------------
plt.figure()
plt.plot(sol.t, sol.y[0], label=r'$\theta_1$')
plt.plot(sol.t, sol.y[1], label=r'$\theta_2$')
plt.xlabel('Temps [s]')
plt.ylabel('Angle [rad]')
plt.legend()
plt.title('Angles del doble pèndol de barres rígides')
plt.show()

# ---------------------------
# 4️⃣ Dibuix de posicions clau
# ---------------------------
num_posicions = 8  # quantes posicions volem veure
indices = np.linspace(0, len(t_eval)-1, num_posicions, dtype=int)

plt.figure(figsize=(8,6))
ax = plt.gca()
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-2.5, 0.5)
ax.set_aspect('equal')
ax.set_title("Posicions clau del doble pèndol de barres rígides")

colors = plt.cm.viridis(np.linspace(0,1,num_posicions))

for i, idx in enumerate(indices):
    th1 = sol.y[0, idx]
    th2 = sol.y[1, idx]
    
    x1 = L1*np.sin(th1)
    y1 = -L1*np.cos(th1)
    
    x2 = x1 + L2*np.sin(th2)
    y2 = y1 - L2*np.cos(th2)
    
    ax.plot([0, x1, x2], [0, y1, y2], 'o-', color=colors[i], lw=2, label=f't={t_eval[idx]:.1f}s')

ax.legend()
plt.show()