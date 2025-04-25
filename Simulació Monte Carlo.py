import numpy as np
import matplotlib.pyplot as plt

N = 10**4  # nombre de partícules
d = 3      # dimensions
T = 300    # temperatura
M = 10**7  #nombre de passes
m = 1.67e-27      # massa d'una partícula (massa del protó)
k_B = 1.38e-23    # constant de Boltzmann

def gas_ideal(N, d, T, M):
    velocitats = np.random.normal(0, np.sqrt(k_B * T / m), (N, d))  # velocitats inicials

    def energia_particula(v):
        return 0.5 * m * np.sum(v**2) # energia cinètica d'una partícula

    energia_total = 0.5 * m * np.sum(velocitats**2) # energia cinètica total inicial
    energies = []
    for i in range(M):
        part_random = np.random.randint(0, N) # seleccionem una partícula aleatòria
        velocitat_vella = velocitats[part_random].copy()
        E_vella = energia_particula(velocitat_vella) # energia cinètica de la partícula seleccionada

        nova_velocitat = velocitat_vella + np.random.normal(0, np.sqrt(k_B * T / m), d) # nova velocitat aleatòria
        E_nova = energia_particula(nova_velocitat)

        delta_energia = E_nova - E_vella # variació d'energia deguda al canvi de velocitat

        if delta_energia < 0 or np.random.uniform(0, 1) < np.exp(-delta_energia / (k_B * T)):
            velocitats[part_random] = nova_velocitat #si la variació d'energia és negativa o es compleix la condició d'acceptació, actualitzem la velocitat
            energia_total += delta_energia #sumem la variació d'energia a l'energia total
    
        if i % 5000 == 0: # cada 5000 iteracions calculem l'energia total
            energies.append(energia_total)

   # Calculem la capacitat calorífica a cada pas
    energies = np.array(energies)
    E_mean = np.cumsum(energies) / np.arange(1, len(energies)+1)
    E2_mean = np.cumsum(energies**2) / np.arange(1, len(energies)+1)
    C_Vs = (E2_mean - E_mean**2) / (k_B * T**2)

    return C_Vs

C_Vs = gas_ideal(N, d, T, M)

# Valor teòric per comparar
C_V_teoric = (d / 2) * N * k_B

# Gràfic de la capacitat calorífica
plt.plot(C_Vs, label=r'$C_V$ estimada')
plt.axhline(C_V_teoric, color='r', linestyle='--', label=r'$C_V$ teòrica')
plt.xlabel('Iteració (cada 5000 passos)')
plt.ylabel(r'$C_V$ (J/K)')
plt.title('Estimació de la capacitat calorífica')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# Últim valor estimat de la capacitat calorífica
print(f"Capacitat calorífica teòrica: {C_V_teoric:.3e} J/K")
print(f"Capacitat calorífica simulada: {C_Vs[-1]:.3e} J/K")
print(f"Diferència relativa: {100 * abs(C_Vs[-1] - C_V_teoric) / C_V_teoric:.3f} %")