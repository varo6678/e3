# /e/src/lib/placeholders.py

from dataclasses import dataclass

@dataclass
class ParametrosFisicos:
    """Class for keeping track of an item in inventory."""
    T0: float
    T1: float

@dataclass
class NpuntosDireccion:
    Nx: int
    Ny: int

@dataclass
class ParametrosGeometricos:
    R: float
    
@dataclass
class ParametrosComputacionales:
    max_iteraciones: int
    tolerancia: float