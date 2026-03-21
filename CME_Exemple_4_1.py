import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14
})

# %% 1. Definició de variables simbòliques
s, t = sp.symbols('s t', real=True)

# Paràmetres del sistema
m = 1                 # massa (kg)
c_values = [2, 4, 5]  # esmorteeixements (Ns/m)
k = 4                 # rigidesa (N/m)
F_s = 1/s             # Funció esglaó en el domini de Laplace (F(t)=1)

t_eval = np.linspace(0, 8, 1000)
results_x = []
results_v = []
results_fr = []
results_fc = []
results_fb = []

# %% 2. Bucle de càlcul simbòlic i numèric
for c in c_values:
    # Funció de transferència de la posició X(s)
    X_s = F_s / (m*s**2 + c*s + k)
    Xdot_s = F_s / (m*s**2 + c*s + k) *s

    # Transformada inversa de Laplace per trobar x(t)
    x_t_expr = sp.inverse_laplace_transform(X_s, s, t)
    v_t_expr = sp.inverse_laplace_transform(Xdot_s, s, t)

    # Simplifiquem l'expressió per seguretat
    x_t_expr = sp.simplify(x_t_expr)
    v_t_expr = sp.simplify(v_t_expr)

    # Mostrem l'expressió per pantalla (equivalent a pretty en Matlab)
    print()
    print(f"--> Resultat x(t) per c = {c} <--")
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
    results_fr.append(x_vals*k)
    results_fc.append(v_vals*c)
    results_fb.append(x_vals*k + v_vals*c )

# %% 3. Representació respostes

plt.figure(figsize=(10, 5))

# Gràfica de la Posició
plt.subplot(1, 2, 1)
for i, c in enumerate(c_values):
    plt.plot(t_eval, results_x[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# plt.title(r'Posició $x(t)$')
plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$x\,[m]$')
plt.xlim([0,8])
plt.grid(True)
plt.legend()

# Gràfica de la Velocitat
plt.subplot(1, 2, 2)
for i, c in enumerate(c_values):
    plt.plot(t_eval, results_v[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# plt.title(r'Força elàstica $F_r = k\,x(t)$')
plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$\dot{x}\,[m/s]$')
plt.xlim([0,8])
plt.grid(True)
# plt.legend()

plt.tight_layout()
plt.savefig('CME_Exemple_4_1_SOL.pdf', bbox_inches='tight', transparent=True)
plt.show()

plt.figure(figsize=(10, 5))

# Gràfica de la Força elàstica
plt.subplot(1, 3, 1)
for i, c in enumerate(c_values):
    plt.plot(t_eval, results_fr[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# plt.title(r'Posició $x(t)$')
plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$F_r\,[N]$')
plt.xlim([0,8])
plt.ylim([-0.2,1.4])
plt.grid(True)
# plt.legend()

# Gràfica de la Força esmorteïdor
plt.subplot(1, 3, 2)
for i, c in enumerate(c_values):
    plt.plot(t_eval, results_fc[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# plt.title(r'Força esmorteïdor $F_c = c\,\dot{x}(t)$')
plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$F_c\,[N]$')
plt.xlim([0,8])
plt.ylim([-0.2,1.4])
plt.grid(True)
# plt.legend()

# Gràfica de la Força a la base
plt.subplot(1, 3, 3)
for i, c in enumerate(c_values):
    plt.plot(t_eval, results_fb[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
# plt.title(r'Força a la base $F_r + F_c$')
plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$F_b\,[N]$')
plt.xlim([0,8])
plt.ylim([-0.2,1.4])
plt.grid(True)
# plt.legend()

plt.tight_layout()
plt.savefig('CME_Exemple_4_1_SOL2.pdf', bbox_inches='tight', transparent=True)
plt.show()
