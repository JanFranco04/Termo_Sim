from vpython import *
import numpy as np
import random as rn
import matplotlib.pyplot as plt

# Configuración
Natoms = 500
L = 1.0
Ratom = 0.01
mass = 4E-3 / 6E23
k = 1.38E-23
dt = 1E-4
steps = 700
nu = 50  # Frecuencia Andersen (aumenta para mejor acoplamiento)

# Inicialización
Atoms = []
p = []
apos = []
TAndersen = 300

Presion = []
Temperatura = []
TAndersen_l=[]
Pteorica = []

# VPython scene (opcional)
scene = canvas(width=500, height=500)
scene.range = L
scene.title = 'Gas con termostato de Andersen'

for i in range(Natoms):
    pos = vector(rn.uniform(-L/2,L/2), rn.uniform(-L/2,L/2), rn.uniform(-L/2,L/2))
    apos.append(pos)
    Atoms.append(sphere(pos=pos, radius=Ratom, color=color.cyan))
    sigma = np.sqrt(k * TAndersen / mass)
    vel = vector(rn.gauss(0, sigma), rn.gauss(0, sigma), rn.gauss(0, sigma))
    p.append(mass * vel)

for step in range(steps):
    rate(200)
    
    # Andersen Thermostat
    for i in range(Natoms):
        if rn.random() < nu * dt:
            sigma = np.sqrt(k * TAndersen / mass)
            p[i] = vector(mass * rn.gauss(0, sigma),
                          mass * rn.gauss(0, sigma),
                          mass * rn.gauss(0, sigma))
    
    # Mover partículas
    for i in range(Natoms):
        vel = p[i] / mass
        apos[i] += vel * dt
        Atoms[i].pos = apos[i]
        
    if abs(apos[i].x) > (L/2 - Ratom):
        p[i].x *= -1
    if abs(apos[i].y) > (L/2 - Ratom):
        p[i].y *= -1
    if abs(apos[i].z) > (L/2 - Ratom):
        p[i].z *= -1

    
    # Calcular presión desde el momento total transferido
    E_kin = 0
    for i in range(Natoms):
        v = p[i] / mass
        E_kin += 0.5 * mass * mag2(v)
    
    T_real = (2/3) * E_kin / (Natoms * k)
    P_real = Natoms * k * T_real / (L**3)
    P_th = Natoms * k * TAndersen / (L**3)
    
    Presion.append(P_real)
    Temperatura.append(T_real)
    Pteorica.append(P_th)

    TAndersen += 10  # cambiar temperatura del baño térmico
    TAndersen_l.append(TAndersen)

# Gráfico final
plt.figure()
plt.plot(TAndersen_l, Presion, 'o', color="red", label="Simulada", markersize=4)
plt.plot(TAndersen_l, Pteorica, '-', color="purple", label="P = NkT", markersize=5)
plt.xlabel("Temperatura simulada (K)")
plt.ylabel("Presión (Pa)")
plt.title("Presión vs Temperatura (isócoro)")
plt.grid(True)
plt.legend()
plt.show()
