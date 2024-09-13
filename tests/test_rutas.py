# .

import pytest

# Rutas de e.
from e.src.lib.constants import Rutas as RutasE

class TestRutas:
    
    @pytest.mark.parametrize("ruta, esperado", [
        (RutasE.RUTA_PAQUETE, "e")
    ])
    def test_de_rutas(self, ruta, esperado): 
        assert ruta.split("/")[-1] == esperado