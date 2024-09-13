# /e/lib/clases.py

from ..core._abstractas import Parametros
from ..core._typing import (
    InputsLike
)

class ParametrosProblema(Parametros):
    
    """
    Ejemplo
    ---
    >>> Temperatura.T0 = 0
    >>> Temperatura.T1 = 50
    >>> NpuntosDireccion.Nx = 100
    >>> NpuntosDireccion.Ny = 100
    >>> SemiCirculoParametros.R = 1

    >>> inputs = {
    >>>     'T0' : Temperatura.T0,
    >>>     'T1' : Temperatura.T1,
    >>>     'Nx' : NpuntosDireccion.Nx,
    >>>     'Ny' : NpuntosDireccion.Ny,
    >>>     'R' : SemiCirculoParametros.R
    >>> }

    >>> params = ParametrosProblema(dict_parametros=inputs)
    >>> params.print_parametros
    """
    
    def __init__(self, dict_parametros: InputsLike) -> None:
        self.inputs = dict_parametros
        
    def __repr__(self) -> str:
        return f"ParametrosProblema({list(self.inputs.keys())})"
        
    @property
    def print_parametros(self):
        for parametro, valor in self.inputs.items():
            print(f"{parametro:3} | {valor:6}")   