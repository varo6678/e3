# /e/src/lib/constants.py

import os
import sys
from pathlib import Path
from dataclasses import dataclass

__all__ = [
    'Rutas'
]

@dataclass(frozen=True)
class Rutas:
    RUTA_PAQUETE: str = str(Path(__file__).resolve().parents[2])