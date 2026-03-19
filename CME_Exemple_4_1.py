# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (manim)
#     language: python
#     name: manim
# ---

# %% colab={"base_uri": "https://localhost:8080/", "height": 976} executionInfo={"elapsed": 10147, "status": "ok", "timestamp": 1773084199968, "user": {"displayName": "Joan Andreu Mayugo Majo", "userId": "06103893660111504434"}, "user_tz": -60} id="2Aw5F2h_zrM6" outputId="10d14443-1281-4c0e-fcd4-c4ccfe419e39"
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

plt.rcParams.update({
    "text.usetex": True,      # Utilitza LaTeX per a tot el text
    "font.family": "serif",   # Fonts serif com LaTeX
    "font.size": 14
})

# %% colab={"base_uri": "https://localhost:8080/", "height": 976} executionInfo={"elapsed": 10147, "status": "ok", "timestamp": 1773084199968, "user": {"displayName": "Joan Andreu Mayugo Majo", "userId": "06103893660111504434"}, "user_tz": -60} id="2Aw5F2h_zrM6" outputId="10d14443-1281-4c0e-fcd4-c4ccfe419e39"
# 1. Definició de variables simbòliques
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

# %% colab={"base_uri": "https://localhost:8080/", "height": 976} executionInfo={"elapsed": 10147, "status": "ok", "timestamp": 1773084199968, "user": {"displayName": "Joan Andreu Mayugo Majo", "userId": "06103893660111504434"}, "user_tz": -60} id="2Aw5F2h_zrM6" outputId="10d14443-1281-4c0e-fcd4-c4ccfe419e39"
# 2. Bucle de càlcul simbòlic i numèric
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
    results_fr.append(x_vals * k)

# %% colab={"base_uri": "https://localhost:8080/", "height": 976} executionInfo={"elapsed": 10147, "status": "ok", "timestamp": 1773084199968, "user": {"displayName": "Joan Andreu Mayugo Majo", "userId": "06103893660111504434"}, "user_tz": -60} id="2Aw5F2h_zrM6" outputId="10d14443-1281-4c0e-fcd4-c4ccfe419e39"
plt.figure(figsize=(10, 5))

# Gràfica de la Posició
plt.subplot(1, 2, 1)
for i, c in enumerate(c_values):
    plt.plot(t_eval, results_x[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
plt.title(r'Posició $x(t)$')
plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$x\,[m]$')
plt.grid(True)
plt.legend()

# Gràfica de la Força elàstica
plt.subplot(1, 2, 2)
for i, c in enumerate(c_values):
    plt.plot(t_eval, results_v[i], label=fr'$c = {c}\ \mathrm{{Ns/m}}$')
plt.title(r'Velocitat $\dot{x}(t)$')
plt.xlabel(r'$t\,[s]$')
plt.ylabel(r'$\dot{x}\,[m/s]$')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig('grafica.pdf', bbox_inches='tight', transparent=True)
plt.show()
