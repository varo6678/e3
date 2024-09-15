# Imports necesarios
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import logging
import argparse

# Configuración del logger
def define_logger(logger_name='mna', logger_level='INFO'):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logger_level)
    ch = logging.StreamHandler()
    ch.setLevel(logger_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

informer = define_logger(logger_name='mna', logger_level='INFO')

def resolver_laplace_polar(Nr, Ntheta, R, T0, T1, tolerancia=1e-6, max_iter=10000, omega=1.0):
    # Crear la malla
    dr = R / (Nr - 1)
    dtheta = np.pi / (Ntheta - 1)
    r = np.linspace(0, R, Nr)
    theta = np.linspace(0, np.pi, Ntheta)
    R_grid, Theta_grid = np.meshgrid(r, theta, indexing='ij')

    # Inicializar la matriz de temperaturas
    u = np.zeros((Nr, Ntheta))

    # Condiciones de frontera
    u[-1, :] = T1  # Borde circular (r = R)
    u[:, 0] = T0   # Diámetro (θ = 0)
    u[:, -1] = T0  # Diámetro (θ = pi)

    # Iteraciones
    convergencia = False
    iter_count = 0

    while not convergencia and iter_count < max_iter:
        u_old = u.copy()
        iter_count += 1

        for i in range(1, Nr - 1):
            r_i = r[i]
            if r_i == 0:
                continue  # Evitar división por cero
            beta = (r_i * dtheta / dr) ** 2
            denom = 2 * (1 + beta)
            for j in range(1, Ntheta - 1):
                u_new = (1 / denom) * (u[i+1, j] + u[i-1, j] + beta * (u[i, j+1] + u[i, j-1]))
                u[i, j] = u[i, j] + omega * (u_new - u[i, j])

        # Manejar el centro (r = 0)
        u[0, :] = np.mean(u[1, :])  # Asumir simetría radial

        # Verificar convergencia
        max_diff = np.max(np.abs(u - u_old))
        if max_diff < tolerancia:
            convergencia = True
            informer.info(f'Convergencia alcanzada en {iter_count} iteraciones con diferencia máxima {max_diff:.2e}')

    if not convergencia:
        informer.warning(f'No se alcanzó la convergencia después de {max_iter} iteraciones')

    return r, theta, u

# Función para resolver la ecuación de Laplace en coordenadas cartesianas
def resolver_laplace_cartesiano(Nx, Ny, R, T0, T1, tolerancia=1e-6, max_iter=10000, omega=1.0):
    # Crear la malla
    dx = R / (Nx - 1)
    dy = R / (Ny - 1)
    x = np.linspace(0, R, Nx)
    y = np.linspace(0, R, Ny)
    X_grid, Y_grid = np.meshgrid(x, y, indexing='ij')

    # Inicializar la matriz de temperaturas
    u = np.zeros((Nx, Ny))

    # Aplicar condiciones de frontera
    # Borde circular (x^2 + y^2 = R^2)
    for i in range(Nx):
        for j in range(Ny):
            if x[i]**2 + y[j]**2 >= R**2:
                u[i, j] = T1

    # Diámetro (y = 0)
    u[:, 0] = T0

    # Iteraciones
    convergencia = False
    iter_count = 0

    while not convergencia and iter_count < max_iter:
        u_old = u.copy()
        iter_count += 1

        for i in range(1, Nx - 1):
            for j in range(1, Ny - 1):
                # Verificar si dentro.
                if x[i]**2 + y[j]**2 < R**2 and y[j] >= 0:
                    u_new = 0.25 * (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1])
                    u[i, j] = u[i, j] + omega * (u_new - u[i, j])

        # Verificar convergencia
        max_diff = np.max(np.abs(u - u_old))
        if max_diff < tolerancia:
            convergencia = True
            informer.info(f'Convergencia alcanzada en {iter_count} iteraciones con diferencia máxima {max_diff:.2e}')

    if not convergencia:
        informer.warning(f'No se alcanzó la convergencia después de {max_iter} iteraciones')

    return x, y, u

# Función para convertir la tabla a LaTeX
def convertir_tabla_a_latex(df: pd.DataFrame, ruta_salida: str):
    latex_code = df.to_latex(index=False)
    with open(ruta_salida, 'w') as f:
        f.write(latex_code)
    informer.info(f"Tabla en formato LaTeX guardada en {ruta_salida}")

if __name__ == '__main__':
    
    # Parámetros físicos y numéricos
    parser = argparse.ArgumentParser(description='Solución de la ecuación de Laplace en una lámina semicircular.')
    parser.add_argument('--verbosity', type=str, default='INFO', help='Nivel de verbosidad del logger.')
    argumentos_parseados = parser.parse_args()
    informer.setLevel(argumentos_parseados.verbosity)

    # Rutas de salida
    TEMATICA = 'ecuacion_laplace'
    FORMATO_GRAFICAS = '.png'
    OUTPUTS = 'OUTPUTS'
    RUTA_OUTPUTS = f"./{OUTPUTS}"
    RUTA_OUTPUTS_LATEX = f"./e3_latex/figuras"
    Path(RUTA_OUTPUTS).mkdir(parents=True, exist_ok=True)
    Path(RUTA_OUTPUTS_LATEX).mkdir(parents=True, exist_ok=True)

    # Parámetros del problema
    R = 1.0  # Radio de la lámina semicircular
    T0 = 0 # Temperatura en el diámetro
    T1 = 50   # Temperatura en el borde circular

    # Parámetros numéricos para coordenadas polares
    Nr = 50  # Aumentamos la resolución para mayor precisión
    Ntheta = 50
    tolerancia = 1e-6

    # Resolver en coordenadas polares
    r_polar, theta_polar, u_polar = resolver_laplace_polar(Nr, Ntheta, R, T0, T1, tolerancia=tolerancia)

    # Obtener las temperaturas en los puntos específicos (c)
    puntos_r = [0, R/4, R/2, 3*R/4, R]
    puntos_theta = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]

    datos_puntos = []

    for r_val in puntos_r:
        r_idx = np.argmin(np.abs(r_polar - r_val))
        for theta_val in puntos_theta:
            theta_idx = np.argmin(np.abs(theta_polar - theta_val))
            temp = u_polar[r_idx, theta_idx]
            datos_puntos.append([f"r={r_polar[r_idx]:.2f}, θ={theta_polar[theta_idx]:.2f}", temp])

    # Crear tabla de resultados
    tabla_polar = pd.DataFrame(datos_puntos, columns=['Punto (r, θ)', 'Temperatura (°C)'])

    # Guardar la tabla en CSV y LaTeX
    tabla_polar.to_csv(f"{RUTA_OUTPUTS}/tabla_polar.csv", index=False)
    convertir_tabla_a_latex(tabla_polar, f"{RUTA_OUTPUTS}/{TEMATICA}_tabla_polar.tex")


    # Visualización de la solución numérica en coordenadas polares
    plt.figure(figsize=(8, 6))
    R_grid, Theta_grid = np.meshgrid(r_polar, theta_polar, indexing='ij')
    X = R_grid * np.cos(Theta_grid)
    Y = R_grid * np.sin(Theta_grid)
    plt.contourf(X, Y, u_polar, levels=50, cmap='hot')
    plt.colorbar(label='Temperatura (°C)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Distribución de temperatura en coordenadas polares')
    plt.axis('equal')
    plt.savefig(f"{RUTA_OUTPUTS}/{TEMATICA}_polar{FORMATO_GRAFICAS}")
    plt.savefig(f"{RUTA_OUTPUTS_LATEX}/{TEMATICA}_polar{FORMATO_GRAFICAS}")

    # Parámetros numéricos para coordenadas cartesianas
    Nx = 50
    Ny = 50

    # Resolver en coordenadas cartesianas
    x_cart, y_cart, u_cart = resolver_laplace_cartesiano(Nx, Ny, R, T0, T1, tolerancia=tolerancia)

    # Obtener las temperaturas en los puntos específicos (d)
    puntos_x = [0, R/4, R/2, 3*R/4, R]
    puntos_y = [0, R/4, R/2, 3*R/4, R]

    datos_puntos_cart = []

    for x_val in puntos_x:
        x_idx = np.argmin(np.abs(x_cart - x_val))
        for y_val in puntos_y:
            y_idx = np.argmin(np.abs(y_cart - y_val))
            # Verificar si el punto está dentro del semicı́rculo
            if x_cart[x_idx]**2 + y_cart[y_idx]**2 <= R**2 and y_cart[y_idx] >= 0:
                temp = u_cart[x_idx, y_idx]
                datos_puntos_cart.append([f"x={x_cart[x_idx]:.2f}, y={y_cart[y_idx]:.2f}", temp])

    # Crear tabla de resultados
    tabla_cartesiano = pd.DataFrame(datos_puntos_cart, columns=['Punto (x, y)', 'Temperatura (°C)'])

    # Guardar la tabla en CSV y LaTeX
    tabla_cartesiano.to_csv(f"{RUTA_OUTPUTS}/tabla_cartesiano.csv", index=False)
    convertir_tabla_a_latex(tabla_cartesiano, f"{RUTA_OUTPUTS}/{TEMATICA}_tabla_cartesiano.tex")

    # Visualización de la solución numérica en coordenadas cartesianas
    plt.figure(figsize=(8, 6))
    X_grid, Y_grid = np.meshgrid(x_cart, y_cart, indexing='ij')
    plt.contourf(X_grid, Y_grid, u_cart, levels=50, cmap='hot')
    plt.colorbar(label='Temperatura (°C)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Distribución de temperatura en coordenadas cartesianas')
    plt.axis('equal')
    plt.savefig(f"{RUTA_OUTPUTS}/{TEMATICA}_cartesiano{FORMATO_GRAFICAS}")
    plt.savefig(f"{RUTA_OUTPUTS_LATEX}/{TEMATICA}_cartesiano{FORMATO_GRAFICAS}")

    informer.info("Cálculo completado y resultados guardados.")
