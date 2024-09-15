# /e/src/lib/parsers.py

import argparse
from .logger import (
    dict_log_level,
    dict_level_log
)
from ..core._typing import Any

def define_parser(mensaje_descripcion: str = "Este script procesa datos para MNA.") -> Any:

    parser = argparse.ArgumentParser(description=mensaje_descripcion)

    parser.add_argument(
        "-vsy", "--verbosity",
        type=int,
        choices=[level for level in dict_log_level.values()],
        default='INFO',
        help=f"Nivel de verbosidad {list(dict_log_level.items())}"
    )

    parser.add_argument(
        "-sp", "--show_plots", 
        action="store_true",
        help="Muestra los plots del script."
    )

    parser.add_argument(
        "-pl", "--parallel", 
        action="store_true",
        help="Hace los calculos (los que procedan) en paralelo."
    )
    
    return parser