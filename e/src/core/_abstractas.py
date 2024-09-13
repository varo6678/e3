# /e/core/_abstractas.py

from abc import ABC, abstractmethod


class Parametros(ABC):
    
    def __init__(self) -> None:
        pass
    
    @property
    @abstractmethod
    def print_parametros(self):
        pass