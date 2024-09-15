# e

Este repositorio contiene scripts para ayudar a mejorar paulatinamente al paquete e.

## Autor

Álvaro Monforte Marín

## Recomendaciones

- Primero crear un entorno  
  `conda create -n env_e python=3.12`

- Despues activarlo  
  `conda activate env_e`

- Instalar dependencias  
  `pip install -r requirements.txt`

## Requisitos fuertes

- Python 3.12.4
- NumPy
- Matplotlib
- Pandas
- Jinja2
- pytest
- SciPy

## Requisitos débiles

- `requirements.txt`

## Estructura del Repositorio

- `backups`: Contiene copias y pruebas de los scripts.
- `bash_utils`: Contiene scripts de Bash para automatizar tareas.
- `OUTPUTS`: Contiene los resultados de los cálculos y gráficos generados por los scripts.
- `tests`: Contiene scripts de prueba para las funciones.
- `e` : contiene el engine basico del proyecto.
- `e*_latex` : contiene los scripts de latex para la generacion de los informes.
- `pyutils`: Contiene scripts de Python con funciones útiles.

## Cómo Ejecutar

1. Asegúrate de tener todas las bibliotecas requeridas instaladas.
2. p.e. `python e3_final -vsy 10`

## Qué Hace Cada Script

- e3_final.py: Resuelve en cierta medida el ejercicio 3.
- e4_final.py: Resuelve en cierta medida el ejercicio 4.
- e5_final.py: Resuelve en cierta medida el ejercicio 5.


## Outputs

Se encontraran en la carpeta autogenerada por el script si no existe: `OUTPUTS`