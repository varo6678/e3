# /e/src/lib/placeholders.py

from ..core._typing import Optional
from ..core._abstractas import Parametros

__all__ = [
    'ParametrosFisicos',
    'NpuntosDireccion',
    'ParametrosGeometricos',
    'ParametrosComputacionales'
]

class ParametrosFisicos(Parametros):
    pass

class NpuntosDireccion(Parametros):
    pass

class ParametrosGeometricos(Parametros):
    pass

class ParametrosComputacionales(Parametros):
    pass