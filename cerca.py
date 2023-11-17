import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

N = 4  # Cambiar el nuemro de atomos, debe ser par porque se le agregan son anilos
V_1 = sp.symbols('V_1')
V_2 = sp.symbols('V_2')
g_1 = sp.symbols('g_1')

equations = []

# Inicializar las ecuaciones
for i in range(1, N + 1):
    if i == 1:
        equations.append(f'G{i}{N} = g_1 V_1 G{i+1}{N}')
    elif i == N:
        equations.append(f'G{i}{N} = g_1 + g_1 V_1 G{i-1}{N}')
    elif i == 2:
        equations.append(f'G{i}{N} = g_1 V_2 G{i+1}{N} + g_1 V_1 G{i-1}{N}')
    elif 2 <= i <= N-2:
        equations.append(f'G{i}{N} = g_1 V_2 G{i+1}{N} + g_1 V_2 G{i-1}{N}')
    else:
        equations.append(f'G{i}{N} = g_1 V_1 G{i+1}{N} + g_1 V_2 G{i-1}{N}')

# Imprimir la lista de ecuaciones
#print("Lista de ecuaciones:")
#for equation in equations:
#    print(equation)

# Parte nueva para crear la matriz
matrix = sp.zeros(N, N)

# Llenar la matriz con términos negativos fuera de la diagonal principal
for i in range(N):
    for j in range(N):
        if i == j:
            matrix[i, j] = 1
        elif i == j + 1:
            matrix[i, j] = -V_1 * g_1 if (j == 0 or j == N - 2) else -V_2 * g_1
        elif i + 1 == j:
            matrix[i, j] = -V_1 * g_1 if (i == 0 or i == N - 2) else -V_2 * g_1


# Parte nueva para calcular el determinante
determinante = matrix.det()

# Parte nueva para cambiar la primera columna
new_column = sp.Matrix([0]*(N-1) + [g_1])
matrix[:, 0] = new_column

# Calcular e imprimir el determinante de la matriz modificada
determinante_modificado = matrix.det()

# Parte nueva para cambiar la primera columna
new_first_column = sp.Matrix([1] + [-V_1 * g_1] + [0] * (N - 2))
matrix[:, 0] = new_first_column

# Parte nueva para cambiar la última columna
new_last_column = sp.Matrix([0] * (N - 1) + [g_1])
matrix[:, -1] = new_last_column

# Calcular e imprimir el determinante de la matriz con las columnas cambiadas
determinante_columnas_modificadas = matrix.det()

# Cálculo de G_{1N} y G_{NN}
G_1N = determinante_modificado / determinante
G_NN = determinante_columnas_modificadas / determinante

# Definir los parámetros
v = 1.0
gamma = 1.0
E_ini = -2.5
E_end = 2.5
N_pasos = 5000

E_step = (E_end - E_ini) / N_pasos
eta = E_step / 2 #eta = La mitad del paso

# Definir el rango de E
E_values = np.arange(E_ini, E_end, E_step)

# Inicializar la lista para almacenar los valores de tau(E)
tau_values = []

# Funciones para calcular V_1, V_2, g_1 y g_0 en función de E
def calculate_V1(g0, v):
    return (2 * g0**3 * v**3) / (g0 - g0**3 * v**2 + 2 * g0**2 * v)

def calculate_V2(g0, v):
    return (g0 * v - g0**3 * v**3) / (g0 - g0**3 * v**2 + 2 * g0**2 * v)

def calculate_g1(g0, v):
    return (g0 - g0**3 * v**2 + 2 * g0**2 * v) / (1 - 3 * g0**2 * v**2)

def calculate_g0(E, eta):
    return 1 / (E - 1j * eta)

# Calcular los valores de tau(E)
for i,E in enumerate(E_values):
    g0 = calculate_g0(E, eta)
    V1 = calculate_V1(g0, v)
    V2 = calculate_V2(g0, v)
    g1 = calculate_g1(g0, v)

    # Calcular los determinantes y sustituir los valores
    G_1N_value = (determinante_modificado / determinante).subs({V_1: V1, V_2: V2, g_1: g1})
    G_NN_value = (determinante_columnas_modificadas / determinante).subs({V_1: V1, V_2: V2, g_1: g1})

    tau_E = gamma**2*(abs(G_1N_value / ((1 + 1j * gamma * G_NN_value / 2)**2 + (G_1N_value)**2 * gamma**2 / 4))**2).as_real_imag()[0]

    tau_values.append( float( sp.simplify(tau_E) ) )
    
    
plt.plot(E_values, tau_values, label='$\\tau(E)$')
plt.xlabel('E')
plt.ylabel('$\\tau(E)$')
#plt.legend()
plt.grid(True)
plt.savefig('transmision.png')

with open('results.dat', 'w') as results:
  for energy, tau in zip(E_values,tau_values):
    results.write(f"{energy}\t{tau}\n")


















