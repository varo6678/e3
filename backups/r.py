import os
from pathlib import Path

from e.src.lib.placeholders import (
    ParametrosFisicos,
    NpuntosDireccion,
    ParametrosGeometricos,
    ParametrosComputacionales
)

from e import np, plt

from e.src.lib.constants import Rutas

from e.src.core._typing import (
    DictParametrosLike,
    List,
    Any,
    Callable
)

import logging
from logging import _nameToLevel

dict_log_level : DictParametrosLike = _nameToLevel

# Configurar logger
logger = logging.getLogger('mna')
logger.setLevel(dict_log_level['DEBUG'])

# Añadir un StreamHandler para mostrar los logs en la consola
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)

# Añadir un formato básico para los mensajes de log
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)

# Añadir el handler al logger
logger.addHandler(console_handler)

def DiferenciasFinitas2D(
    u: np.ndarray, 
    Nx: int, 
    Ny: int, 
    update_rule: Callable[[np.ndarray, int, int], float],  # Función para actualizar u[i,j]
    tol: float = 1e-6,
    max_iter: int = int(1e4)
    ) -> np.ndarray:
    
    for iteracion in range(max_iter):
        u_old = np.copy(u)
        
        # Iterar sobre los puntos internos.
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                # Llamada a la regla de actualización que depende de la ecuación.
                u[i, j] = update_rule(u_old, i, j)
        
        # Criterio de convergencia
        error = np.max(np.abs(u - u_old))
        if error < tol: 
            print(f"Convergencia alcanzada después de {iteracion} iteraciones.")
            return u, iteracion
        
    print(f"No se alcanzó la convergencia después de {max_iter} iteraciones.")
    return u, max_iter  


def gauss_seidel(u, Nx, Ny, tol=1e-6, max_iter=10000):
    
    for iteracion in range(max_iter):
        u_old = np.copy(u)
        
        # Iterar sobre los puntos internos.
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                # Esquema para la ecuacion de Laplace.
                u[i, j] = 0.25 * (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1])
        
        # Criterio de convergencia
        error = np.max(np.abs(u - u_old))
        if error < tol:
            logger.info(f"Convergencia alcanzada después de {iteracion} iteraciones.")
            return u, iteracion
        
    logger.info(f"No se alcanzó la convergencia después de {max_iter} iteraciones.")
    return u, max_iter

if __name__ == '__main__':
    
    # Temperaturas.
    ParametrosFisicos.T0 = 0
    ParametrosFisicos.T1 = 100

    Temperaturas = ParametrosFisicos.parametros_de_la_clase()

    # Npuntos en la direccion dada.
    NpuntosDireccion.Nx = 5
    NpuntosDireccion.Ny = 5

    Npuntos = NpuntosDireccion.parametros_de_la_clase()

    # Geometria del sistema rectangular | R = { (x, y) | 0 < x < 0.5, 0 < y < 0.5 }.
    ParametrosGeometricos.Lx = 0.5
    ParametrosGeometricos.Ly = 0.5

    Rparams = ParametrosGeometricos.parametros_de_la_clase()

    # Parametros computacionales.
    ParametrosComputacionales.max_iteraciones = int(1e4)
    ParametrosComputacionales.tolerancia = 1e-6

    CompParams = ParametrosComputacionales.parametros_de_la_clase()

    # Resumen de los parametros.

    lista_params_holders: List[List[Any]] = [
        list(Temperaturas.items()),
        list(Npuntos.items()),
        list(Rparams.items()),
        list(CompParams.items())
    ]

    logger.debug(f"Parametros:")
    for item in lista_params_holders:
        logger.debug(item)
        
    u = np.zeros((Npuntos['Nx'], Npuntos['Ny']))
    
    # Condiciones de contorno.

    # Lados a T0.
    u[:, 0] = Temperaturas['T0']
    u[0, :] = Temperaturas['T0']

    # Lados a T1.
    u[-1, :] = np.linspace(0, Temperaturas['T1'], num=Npuntos['Nx'])
    u[:, -1] = np.linspace(0, Temperaturas['T1'], num=Npuntos['Ny'])
    
    u = np.flipud(u)
    
    # Aplicar Gauss-Seidel
    u_final, iteraciones = gauss_seidel(u, Npuntos['Nx'], Npuntos['Ny'], CompParams['tolerancia'], CompParams['max_iteraciones'])

    logger.info(f"Solución después de {iteraciones} iteraciones:")
    logger.info(u_final)
    
    # Voy a crear una lista de resultados de los puntos interiores.
    logger.info(f"{'i':2} | {'wi':5}")
    logger.info(f"{'-'*20}")
    for i in range(1, Npuntos['Nx']-1):
        for j in range(1, Npuntos['Ny']-1):
            # if i != 0 and j != 0 and i != Npuntos['Nx']-1 and j != Npuntos['Ny']-1:
            punto_interior: float = u[i,j]
            logger.info(f"{i} | {u[i,j]:.2f}")