import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14
})

# %% 1. Definició de variables simbòliques

# Paràmetres del sistema
m1 = 1                 # massa m1 (kg)
m2 = 1                 # massa m2 (kg)
m3 = 1                 # massa m3 (kg)
m4 = 1                 # massa m4 (kg)

k1 = 1                 # rigidesa (N/m)

c_values = [0, 0.5, 5]                  # esmorteïment (Ns/m)
c_values = [0]                  # esmorteïment (Ns/m)

f0 = 10

t_final = 20 

# %% 2. Definició de variables simbòliques




s, t = sp.symbols('s t', real=True)

F_s = f0/s             # Funció esglaó en el domini de Laplace (F(t)=1)

t_eval = np.linspace(0, t_final, 1000)
results_x = []
results_v = []

# %% 2. Bucle de càlcul simbòlic i numèric
for c1 in c_values:
    
    m_eq = 8/3*m1+2*m2-3/2*m3+m4
    c_eq = 64/9*c1
    k_eq = 80/9*k1
    
    
    # Funció de transferència de la posició X(s)
    X_s    = F_s / (m_eq*s**2 + c_eq*s + k_eq)
    Xdot_s = F_s / (m_eq*s**2 + c_eq*s + k_eq)*s

    # Transformada inversa de Laplace per trobar x(t)
    x_t_expr = sp.inverse_laplace_transform(X_s, s, t)
    v_t_expr = sp.inverse_laplace_transform(Xdot_s, s, t)

    # Simplifiquem l'expressió per seguretat
    x_t_expr = sp.simplify(x_t_expr)
    v_t_expr = sp.simplify(v_t_expr)

    # Mostrem l'expressió per pantalla (equivalent a pretty en Matlab)
    print()
    print(f"--> Resultat x(t) per c = {c1} <--")
    print( "---------------------------------")
    sp.pprint(x_t_expr)
    print()

    # Convertim l'expressió simbòlica a una funció numèrica (lambdify)
    x_func = sp.lambdify(t, x_t_expr, modules=['numpy'])
    v_func = sp.lambdify(t, v_t_expr, modules=['numpy'])

    # Calculem valors numèrics
    x_vals = x_func(t_eval)
    v_vals = v_func(t_eval)
    results_x.append(x_vals)
    results_v.append(v_vals)

# %% 3. Representació respostes

plt.figure(figsize=(10, 4))

# Gràfica de la Posició
plt.subplot(1, 2, 1)
for i, c1 in enumerate(c_values):
    plt.plot(t_eval, results_x[i], label=fr'$c = {c1}\ \mathrm{{Ns/m}}$')
# plt.title(r'Posició $x(t)$')

plt.axhline(y=f0/k_eq, color='black', linestyle='--', label=r'$f_0/k_{eq}$')

plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$x_4\,[m]$')
plt.xlim([0,t_final])
plt.ylim([-f0/k_eq*0.05,f0/k_eq*2.8])

y0 = 0
y1 = f0 / k_eq
y2 = 2 * f0 / k_eq
plt.yticks(
    [y0, y1, y2],
    [0, r'$\frac{f_0}{k_{eq}}$', r'$2\frac{f_0}{k_{eq}}$'])

omega_n = np.sqrt(k_eq / m_eq)
T = 2*np.pi / omega_n
xticks = [T, 2*T, 3*T, 4*T]
lxticks = [ r'$T$', r'$2T$', r'$3T$', r'$4T$',]

plt.xticks(
    xticks, lxticks)

plt.grid(True)
plt.legend()

# Gràfica de la Velocitat
plt.subplot(1, 2, 2)
for i, c1 in enumerate(c_values):
    plt.plot(t_eval, results_v[i], label=fr'$c = {c1}\ \mathrm{{Ns/m}}$')
# plt.title(r'Força elàstica $F_r = k\,x(t)$')

plt.axhline(y=0, color='black', linestyle='--', label=r'$f_0/k_{eq}$')

plt.xlim([0,t_final])

y1 = -f0 / k_eq *omega_n
y0 = 0
y2 =  f0 / k_eq *omega_n
plt.yticks(
    [y1, y0, y2],
    [r'$-\frac{f_0}{\sqrt{k_{eq} m_{eq}}}$', 0,  r'$\frac{f_0}{\sqrt{k_{eq} m_{eq}}}$'])

plt.xticks(
    xticks, lxticks)

plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$\dot{x}_4\,[m/s]$', labelpad=-32)
plt.grid(True)

plt.tight_layout()
plt.savefig('CME_AC_4_9_SOL.pdf', bbox_inches='tight', transparent=True)
plt.show()

# plt.figure(figsize=(10, 5))

# # Gràfica de la Força elàstica
# plt.subplot(1, 3, 1)
# for i, c in enumerate(c_values):
#     plt.plot(t_eval, results_fr[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# # plt.title(r'Posició $x(t)$')
# plt.xlabel(r'$t\,[s]$')
# plt.ylabel(r'$F_r\,[N]$')
# plt.xlim([0,8])
# plt.ylim([-0.2,1.4])
# plt.grid(True)
# # plt.legend()

# # Gràfica de la Força esmorteïdor
# plt.subplot(1, 3, 2)
# for i, c in enumerate(c_values):
#     plt.plot(t_eval, results_fc[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# # plt.title(r'Força esmorteïdor $F_c = c\,\dot{x}(t)$')
# plt.xlabel(r'$t\,[s]$')
# plt.ylabel(r'$F_c\,[N]$')
# plt.xlim([0,8])
# plt.ylim([-0.2,1.4])
# plt.grid(True)
# # plt.legend()

# # Gràfica de la Força a la base
# plt.subplot(1, 3, 3)
# for i, c in enumerate(c_values):
#     plt.plot(t_eval, results_fb[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# # plt.title(r'Força a la base $F_r + F_c$')
# plt.xlabel(r'$t\,[s]$')
# plt.ylabel(r'$F_b\,[N]$')
# plt.xlim([0,8])
# plt.ylim([-0.2,1.4])
# plt.grid(True)
# # plt.legend()

# plt.tight_layout()
# plt.savefig('CME_AC_4_9_SOL.pdf', bbox_inches='tight', transparent=True)
# plt.show()
