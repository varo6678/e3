# /e/src/lib/logger.py

import logging
from logging import _nameToLevel, _levelToName
from e.src.core._typing import (
    DictParametrosLike,
)

__all__ = [
    'dict_log_level',
    'dict_level_log',
    'define_logger'
]

dict_log_level: DictParametrosLike = _nameToLevel
dict_level_log: DictParametrosLike = _levelToName

def define_logger(logger_name: str, logger_level: str = 'DEBUG'):
    
    # Configurar logger.
    logger = logging.getLogger(logger_name)
    logger.setLevel(dict_log_level[logger_level])

    console_handler = logging.StreamHandler()

    # Añadir un formato básico para los mensajes de log.
    formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)

    # Añadir el handler al logger.
    logger.addHandler(console_handler)
    
    return logger