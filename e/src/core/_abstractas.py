# /e/src/core/_abstractas.py

from abc import ABC, abstractmethod
from ._typing import (
    ListaStringsLike,
    DictParametrosLike
)


class Parametros(ABC):
    
    """
    # Explicacion
    Esta clase pretende facilitar el uso de guardado de parametros.
    
    ## Example
    >>> class SDEModelParameters(Parameters):
    >>>     mu = 2
    >>>     sigma = 1
    >>>     X0 = 1
    >>> 
    >>> params = SDEModelParameters()
    >>> print(params)  # Salida esperada: "Parameters: SDEModelParameters"
    >>> print(params.nombre) # Salida esperada: "SDEModelParameters"
    >>> print(params.parametros_de_la_clase())  # {'mu': 2, 'sigma': 1, 'X0': 1}
    """
    
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self) -> str:
        base_class_name = self.__class__.__bases__[0].__name__
        return f"{base_class_name}: {self.__class__.__name__}"

    @property
    def nombre(self) -> str:
        return self.__repr__().split(':')[1].strip()

    @classmethod
    def lista_de_funciones_prop_de_una_clase(cls) -> ListaStringsLike:
        return [p for p in dir(cls) if isinstance(getattr(cls, p), property)]

    @classmethod
    def parametros_de_la_clase(cls) -> DictParametrosLike:
        # Obtener las propiedades de la clase
        property_names = cls.lista_de_funciones_prop_de_una_clase()
        
        # Obtener todos los atributos que no sean métodos ni propiedades internas
        internal_variables_dict = {k: v for k, v in vars(cls).items() if not k.startswith("__")}
        
        # Excluir las propiedades y los atributos internos que comienzan con "_"
        store_keys = [] + property_names
        for key in internal_variables_dict.keys():
            if key.startswith('_'): 
                store_keys.append(key)
        for key in store_keys:
            if key in internal_variables_dict:
                del internal_variables_dict[key]
        
        return internal_variables_dict