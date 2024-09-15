import numpy as np
import os
from pathlib import Path

from e.src.lib.placeholders import (
    ParametrosFisicos,
    NpuntosDireccion,
    ParametrosGeometricos,
    ParametrosComputacionales
)

from e.src.lib.constants import Rutas

from e.src.core._typing import (
    InputsLike,
    List,
    Any
)

import logging
from logging import _nameToLevel

dict_log_level : InputsLike = _nameToLevel

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

if __name__ == '__main__':
    
    # Crear parámetros físicos y otros.
    temperaturas = ParametrosFisicos(T0=0, T1=50)
    Npuntos = NpuntosDireccion(Nx=5, Ny=5)
    SCparams = ParametrosGeometricos(R=1.0, THETA=np.pi)
    CompParams = ParametrosComputacionales(max_iteraciones=int(1e4), tolerancia=1e-6)
    
    lista_params_holders: List[List[Any]] = [
        list(temperaturas.__dict__.items()),
        list(Npuntos.__dict__.items()), 
        list(SCparams.__dict__.items()), 
        list(CompParams.__dict__.items())
    ]
    
    logger.debug(f"Parametros:")
    for item in lista_params_holders:
        logger.debug(item)
        
    
    # Empiezo el problema.
    u = np.zeros(shape=(Npuntos.Nx, Npuntos.Ny))

    def respetar_condiciones_frontera(u: np.ndarray):
        
        # Dirichlet.
        u[:, 0] = temperaturas.T0   # En theta = 0
        u[:, -1] = temperaturas.T0  # En theta = pi
        u[-1, :] = temperaturas.T1  # En r = 1
        
        # Neumann.
        # u[0, :] = u[1, :]
        
        return u

    u = respetar_condiciones_frontera(u)
    
    
    def metodo_sor(
        u,
        Nx: int,
        Ny: int,
        Dx: float,
        Dy: float,
        omega: float,
        maximas_iteraciones: int,
        tolerancia: float,
    ) -> np.ndarray:
        
        # Iteraciones del método SOR
        for it in range(maximas_iteraciones):
            u_old = np.copy(u)
            
            # Actualizar solo los puntos interiores de la malla
            for i in range(1, Nx - 1):  # Excluir las condiciones en r = 0 y r = 1
                r_i = i * Dx
                for j in range(1, Ny - 1):  # Excluir las condiciones en theta = 0 y theta = pi
                    u_new = (1 / (2 * (1/Dx**2 + 1/(r_i**2 * Dy**2)))) * (
                        (u[i+1, j] + u[i-1, j]) / Dx**2 + 
                        (1 / (r_i**2)) * (u[i, j+1] + u[i, j-1]) / Dy**2
                    )
                    
                    # Sobrerrelajación
                    u[i, j] = (1 - omega) * u[i, j] + omega * u_new
            
            logger.info(f"\n")
            logger.info(f"Iteracion: {it}")
            logger.info(f"{u}")
            
            # Criterio de convergencia.
            diff = np.max(np.abs(u - u_old))
            if diff < tolerancia:
                print(f'Convergencia alcanzada en {it} iteraciones.')
                break
        else:
            print('No se alcanzó la convergencia.')

        return u

    solution = metodo_sor(
        u,
        Nx=Npuntos.Nx, 
        Ny=Npuntos.Ny,
        Dx=SCparams.R/Npuntos.Nx, 
        Dy=SCparams.THETA/Npuntos.Ny,
        omega=1.5,
        tolerancia=CompParams.tolerancia,
        maximas_iteraciones=CompParams.max_iteraciones
    )
    
