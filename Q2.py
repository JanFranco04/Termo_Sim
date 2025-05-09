from vpython import *
import numpy as np
import random as rn
import matplotlib.pyplot as plt

# ------------------- Parámetros físicos y simulación -------------------
Natoms = 500
L = 1.0
Ratom = 0.03
mass = 4E-3 / 6E23
k = 1.38E-23
dt = 1E-4
steps = 1000
nu = 15              # Frecuencia de colisión con el termostato
TAndersen = 300       # Temperatura fija del baño térmico
A = L * L             # Área de una pared

# ------------------- Inicialización -------------------
Atoms = []
p = []
apos = []

Presion = []
Temperatura = []
Pteorica = []
TAndersen_l=[]
tiempos = []

# VPython (opcional)
scene = canvas(width=500, height=500)
scene.range = L
scene.title = 'Gas de esferas duras (term. Andersen + presión molecular)'

for i in range(Natoms):
    pos = vector(rn.uniform(-L/2, L/2), rn.uniform(-L/2, L/2), rn.uniform(-L/2, L/2))
    apos.append(pos)
    Atoms.append(sphere(pos=pos, radius=Ratom, color=color.cyan))
    sigma = np.sqrt(k * TAndersen / mass)
    vel = vector(rn.gauss(0, sigma), rn.gauss(0, sigma), rn.gauss(0, sigma))
    p.append(mass * vel)

# ------------------- Simulación principal -------------------
momentum_transfer_x = 0.0  # solo medimos en x para simplificar

for step in range(steps):
    rate(200)

    # Andersen thermostat: colisión aleatoria con el baño térmico
    for i in range(Natoms):
        if rn.random() < nu * dt:
            sigma = np.sqrt(k * TAndersen / mass)
            p[i] = vector(mass * rn.gauss(-sigma, sigma),
                          mass * rn.gauss(-sigma, sigma),
                          mass * rn.gauss(-sigma, sigma))

    # Mover partículas + colisiones con paredes
    for i in range(Natoms):
        vel = p[i] / mass
        apos[i] += vel * dt
        Atoms[i].pos = apos[i]

        if abs(apos[i].x) > (L/2 - Ratom):
            delta_p = 2 * abs(p[i].x)
            momentum_transfer_x += delta_p
            p[i].x *= -1

        if abs(apos[i].y) > (L/2 - Ratom):
            p[i].y *= -1

        if abs(apos[i].z) > (L/2 - Ratom):
            p[i].z *= -1

    # Calcular temperatura real a partir de E_cinética
    E_kin = sum([0.5 * mass * mag2(p[i] / mass) for i in range(Natoms)])
    T_real = (2/3) * E_kin / (Natoms * k)

    # Presión desde dinámica molecular (transferencia de momento)
    P_molecular = momentum_transfer_x / (dt * A)
    momentum_transfer_x = 0.0  # reset para el siguiente paso

    # Presión teórica
    P_th = Natoms * k * (TAndersen+3) / (L**3)

    # Guardar datos
    tiempos.append(step * dt)
    Presion.append(P_molecular)
    Temperatura.append(T_real)
    Pteorica.append(P_th)
    TAndersen_l.append(TAndersen)
    TAndersen = TAndersen 

# ------------------- Gráficas -------------------

# Presión vs Tiempo
plt.figure()
plt.plot(tiempos, Presion, '-', color='blue', label='Simulada')
plt.plot(tiempos, Pteorica, '-', color='orange', label='Teórica (NkT/V)')
plt.xlabel("Tiempo (s)")
plt.ylabel("Presión (Pa)")
plt.title("Presión simulada vs teórica (función del tiempo)")
plt.legend()
plt.grid(True)
plt.show()

# Presión vs Temperatura
plt.figure()
plt.plot(Temperatura, Presion, 'o', color='red', label='Simulada', markersize=4)
plt.plot(Temperatura, Pteorica, '-', color='purple', label='Teórica (NkT/V)')
plt.xlabel("Temperatura simulada (K)")
plt.ylabel("Presión (Pa)")
plt.title("Presión vs Temperatura (proceso isócoro)")
plt.grid(True)
plt.legend()
plt.show()
