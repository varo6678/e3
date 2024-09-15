# /e/src/core/_typing.py

from ... import np

from typing import (
    Dict, 
    Any,
    Union,
    Tuple,
    List,
    Optional,
    TypeVar,
    Generic,
    Callable
)

__all__ = [
    'Dict', 
    'Any',
    'Union',
    'Tuple',
    'List',
    'Optional',
    'TypeVar',
    'Generic',
    'Callable',
    'HiperparametrosLike',
    'ModuloLike',
    'CoordenadaLike',
    'CoordenadasLike',
    'DictOptionsLIke',
    'InputLike',
    'ResultadosLike',
    'ListaStringsLike',
    'DictParametrosLike',
    'ListaIntLike'
]

type ResultadosLike = Dict[str, Any]
type DictOptionsLIke = Dict[str, Any]
type DictParametrosLike = Dict[str, Any]
type HiperparametrosLike = Dict[str, int | float]
type ModuloLike = float
type CoordenadaLike = float | np.ndarray | int
type CoordenadasLike = Tuple[CoordenadaLike, ...]
type InputLike = Tuple[int | float]
type ListaStringsLike = List[str]
type ListaIntLike = List[int]
