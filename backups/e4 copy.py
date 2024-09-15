# .

# Imports | external.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Imports | internal.
from e.src.lib.logger import (
    define_logger,
    dict_log_level,
    dict_level_log
)
from e.src.lib.placeholders import (
    ParametrosFisicos,
    ParametrosGeometricos,
    ParametrosComputacionales,
    NpuntosDireccion
)
from e.src.lib.constants import Rutas
from e.src.lib.parsers import define_parser

informer = define_logger(logger_name='mna', logger_level='INFO')

def solucion_analitica(x: float, t: float):
    """
    Calcula la solución analítica de la ecuación de onda usando el método de características.
    
    Args:
    - x: posición espacial.
    - t: tiempo.
    
    Returns:
    - u: valor de u(x, t) basado en la solución analítica.
    """
    # Características para el punto R (x - t, x + t).
    x_plus_t = x + t
    x_minus_t = x - t
    
    # Solución particular u_p(x) = (x^2)/2 - (x^4)/12.
    def u_p(x):
        return ((x**2) / 2 )- ((x**4) / 12)
    
    # Solución inicial en t=0, u(x, 0) = x(1 - x).
    def u_inicial(x):
        return x * (1 - x)
    
    # Método de D'Alembert.
    if 0 <= x_plus_t <= 1 and 0 <= x_minus_t <= 1:
        u_h = 0.5 * (u_inicial(x_plus_t) + u_inicial(x_minus_t))
    else:
        u_h = 0  # Fuera del dominio, se asume 0.
    
    # Solución total u(x, t).
    u = u_h + u_p(x)
    return u

def solve_wave_eq(Nx: int, Nt: int, L: float, T: float, cfl: float):
    """
    Resuelve la ecuación de onda hiperbólica con el esquema explícito en diferencias finitas.
    
    Args:
    - Nx: Número de puntos en la dirección espacial (x).
    - Nt: Número de puntos en la dirección temporal (t).
    - L: Longitud del dominio espacial.
    - T: Tiempo total a simular.
    - cfl: Número de Courant (CFL), define la relación entre dt y dx.
    
    Returns:
    - u: Matriz con las soluciones aproximadas.
    - x: Vector de posiciones espaciales.
    - t: Vector de tiempos.
    """
    # Discretización espacial y temporal
    dx = L / (Nx - 1)
    dt = cfl * dx  # Para mantener la estabilidad, dt <= dx/c
    x = np.linspace(0, L, Nx)
    t = np.linspace(0, T, Nt)
    
    # Coeficiente de estabilidad CFL
    r = (dt / dx)**2
    
    # Inicialización de la solución
    u = np.zeros((Nt, Nx))
    
    # Condiciones iniciales
    u[0, :] = x * (1 - x)  # u(x, 0) = x(1 - x)
    
    # Primera iteración: derivada temporal cero
    u[1, :] = u[0, :]  # u_t(x, 0) = 0 implica que u[1, :] = u[0, :]
    
    # Aplicar condiciones de frontera
    u[:, 0] = 0  # u(0, t) = 0
    u[:, -1] = 0  # u(1, t) = 0
    
    # Iteraciones en el tiempo
    for n in range(1, Nt-1):
        for i in range(1, Nx-1):
            u[n+1, i] = (2 * u[n, i] - u[n-1, i] + r * (u[n, i+1] - 2 * u[n, i] + u[n, i-1]) + dt**2 * (1 - x[i]**2))
    
    return u, x, t

def calcular_error(Nx, Nt, L, T, cfl, x_r, t_r, u_analitico):
    u_numerico, x, t = solve_wave_eq(Nx, Nt, L, T, cfl)
    dx = L / (Nx - 1)
    dt = cfl * dx
    i_R = int(x_r / dx)
    n_R = int(t_r / dt)
    u_numerico_R = u_numerico[n_R, i_R]
    return abs(u_analitico - u_numerico_R)


if __name__ == '__main__':
    
    parser = define_parser()
    argumentos_parseados = parser.parse_args()
    informer.setLevel(argumentos_parseados.verbosity)
    
    # Ruta de Outputs:
    OUTPUTS: str = 'OUTPUTS'
    RUTA_OUTPUTS: str = f"{Rutas.RUTA_PAQUETE}/../{OUTPUTS}"
    Path(RUTA_OUTPUTS).mkdir(parents=True, exist_ok=True)
    
    # Primero la solución analítica.
    # Punto R(0.3, 0.1)
    x_r = 0.3
    t_r = 0.1
    u_analitico = solucion_analitica(x_r, t_r)
    informer.info(f"Solución analítica en el punto R(0.3, 0.1): {u_analitico}")
    
    # Solución numérica.
    gparams = ParametrosGeometricos(T=0.1, L=1.0, Nx=10000, Nt=20000)
    cparams = ParametrosComputacionales(cfl=0.15)

    # Resolver la ecuación de onda
    u_numerico, x, t = solve_wave_eq(
        gparams.Nx, 
        gparams.Nt,
        gparams.L, 
        gparams.T, 
        cparams.cfl
    )
    
    dx = gparams.L / (gparams.Nx - 1)
    dt = cparams.cfl * dx
    
    i_R = int(0.3 / dx)
    n_R = int(0.1 / dt)
    
    # Solución numérica en el punto R
    u_numerico_R = u_numerico[n_R, i_R]
    print(f"Solución numérica en el punto R(0.3, 0.1): {u_numerico_R}")

    # Comparación.
    error = abs(u_analitico - u_numerico_R)
    print(f"Error absoluto en el punto R(0.3, 0.1): {error}")

    # # Valores de Nx y Nt para el análisis de convergencia
    # Nx_values = [50, 100, 200, 400, 1000]
    # errors = []

    # for Nx in Nx_values:
    #     Nt = 2 * Nx  # Relación fija entre Nt y Nx
    #     error = calcular_error(
    #         Nx, 
    #         Nt,
    #         gparams.L, 
    #         gparams.T, 
    #         cparams.cfl, 
    #         x_r, 
    #         t_r, 
    #         u_analitico
    #     )
    #     errors.append(error)

    # # Crear una tabla de resultados
    # df_errors = pd.DataFrame({
    #     'Nx': Nx_values,
    #     'Error absoluto': errors
    # })

    # # Mostrar la tabla
    # print(df_errors)

    # # Graficar el error vs Nx
    # plt.figure()
    # plt.loglog(Nx_values, errors, '-o')
    # plt.xlabel('Nx (Resolución espacial)')
    # plt.ylabel('Error absoluto')
    # plt.title('Análisis de convergencia')
    # plt.grid(True)
    # plt.show()


