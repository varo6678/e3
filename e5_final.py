# Imports | external.
from pathlib import Path
from scipy.linalg import solve

# Imports | internal.
from e.src.lib.logger import (
    define_logger,
    dict_log_level,
    dict_level_log
)
from e.src.lib.placeholders import (
    ParametrosFisicos,
    ParametrosGeometricos,
    ParametrosComputacionales,
    NpuntosDireccion
)
from e.src.lib.constants import Rutas
from e.src.lib.parsers import define_parser
from e import plt, np, pd

informer = define_logger(logger_name='mna', logger_level='INFO')


# Definimos funciones para obtener las matrices de rigidez y el vector de carga para elementos lineales y cuadráticos
def stiffness_matrix_linear(n_elements):
    """Genera la matriz de rigidez para elementos finitos lineales con n elementos"""
    n_nodes = n_elements + 1
    K = np.zeros((n_nodes, n_nodes))
    h = 1.0 / n_elements
    
    # Ensamblar la matriz de rigidez
    for i in range(n_elements):
        K[i, i] += 1 / h
        K[i, i + 1] += -1 / h
        K[i + 1, i] += -1 / h
        K[i + 1, i + 1] += 1 / h
    
    return K

def load_vector_linear(n_elements):
    """Genera el vector de carga para elementos finitos lineales con n elementos"""
    n_nodes = n_elements + 1
    f = np.zeros(n_nodes)
    h = 1.0 / n_elements
    
    # Ensamblar el vector de carga
    for i in range(n_elements):
        x_i = i * h
        x_i1 = (i + 1) * h
        f[i] += h / 2 * (x_i + h / 2)
        f[i + 1] += h / 2 * (x_i1 + h / 2)
    
    return f

def apply_boundary_conditions(K, f):
    """Aplica las condiciones de contorno u(0)=0 y du/dx(1)=0"""
    K[0, :] = 0
    K[0, 0] = 1
    f[0] = 0
    
    # La condición du/dx(1) = 0 ya está implícita en el problema (sin flujo)
    return K, f

# Matriz de rigidez para funciones cuadráticas
def stiffness_matrix_quadratic(n_elements):
    """Genera la matriz de rigidez para elementos finitos cuadráticos con n elementos"""
    n_nodes = 2 * n_elements + 1  # 3 nodos por elemento menos 1
    K = np.zeros((n_nodes, n_nodes))
    h = 1.0 / n_elements

    # Ensamblar la matriz de rigidez
    for i in range(n_elements):
        # Indices de los nodos de cada elemento
        idx = [2 * i, 2 * i + 1, 2 * i + 2]
        Ke = np.array([
            [7, -8, 1],
            [-8, 16, -8],
            [1, -8, 7]
        ]) * (1 / (3 * h))

        # Agregar el Ke a la matriz global
        for a in range(3):
            for b in range(3):
                K[idx[a], idx[b]] += Ke[a, b]

    return K

def load_vector_quadratic(n_elements):
    """Genera el vector de carga para elementos finitos cuadráticos con n elementos"""
    n_nodes = 2 * n_elements + 1  # 3 nodos por elemento menos 1
    f = np.zeros(n_nodes)
    h = 1.0 / n_elements

    # Ensamblar el vector de carga
    for i in range(n_elements):
        # Indices de los nodos de cada elemento
        idx = [2 * i, 2 * i + 1, 2 * i + 2]
        fe = np.array([1, 4, 1]) * (h / 6)  # Vector de carga para el término fuente lineal en x

        # Agregar el fe al vector global
        for a in range(3):
            f[idx[a]] += fe[a] * (2 * i + a) * h / 2  # Peso con x medio en cada subintervalo

    return f

if __name__ == '__main__':
    
    # Ruta de Outputs:
    TEMATICA: str = 'galerkin'
    FORMATO_GRAFICAS: str = '.png'
    OUTPUTS: str = 'OUTPUTS'
    RUTA_OUTPUTS: str = f"{Rutas.RUTA_PAQUETE}/../{OUTPUTS}"
    RUTA_OUTPUTS_LATEX: str = f"{Rutas.RUTA_PAQUETE}/../e5_latex/figuras"
    Path(RUTA_OUTPUTS).mkdir(parents=True, exist_ok=True)

    # Resolver para el caso (a) - 4 elementos finitos, interpolación lineal
    n_elements = 4
    K = stiffness_matrix_linear(n_elements)
    f = load_vector_linear(n_elements)
    K_bc, f_bc = apply_boundary_conditions(K, f)

    # Resolución del sistema
    u = solve(K_bc, f_bc)

    # Graficamos la solución obtenida
    x = np.linspace(0, 1, n_elements + 1)

    # Guardamos la solución para la comparación posterior
    solutions = {"4_elements_linear": u}

    # Resolver para el caso (b) - 8 elementos finitos, interpolación lineal
    n_elements = 8
    K = stiffness_matrix_linear(n_elements)
    f = load_vector_linear(n_elements)
    K_bc, f_bc = apply_boundary_conditions(K, f)

    # Resolución del sistema
    u_8_elements = solve(K_bc, f_bc)

    # Graficamos la solución obtenida
    x_8_elements = np.linspace(0, 1, n_elements + 1)

    # Guardamos la solución para la comparación posterior
    solutions["8_elements_linear"] = u_8_elements

    # Resolver para el caso (c) - 2 elementos finitos, interpolación cuadrática
    n_elements = 2
    K_quad = stiffness_matrix_quadratic(n_elements)
    f_quad = load_vector_quadratic(n_elements)
    K_quad_bc, f_quad_bc = apply_boundary_conditions(K_quad, f_quad)

    # Resolución del sistema
    u_quad_2_elements = solve(K_quad_bc, f_quad_bc)

    # Graficamos la solución obtenida
    x_quad_2_elements = np.linspace(0, 1, 2 * n_elements + 1)

    # Guardamos la solución para la comparación posterior
    solutions["2_elements_quadratic"] = u_quad_2_elements

    # Resolver para el caso (d) - 4 elementos finitos, interpolación cuadrática
    n_elements = 4
    K_quad = stiffness_matrix_quadratic(n_elements)
    f_quad = load_vector_quadratic(n_elements)
    K_quad_bc, f_quad_bc = apply_boundary_conditions(K_quad, f_quad)

    # Resolución del sistema
    u_quad_4_elements = solve(K_quad_bc, f_quad_bc)

    # Graficamos la solución obtenida
    x_quad_4_elements = np.linspace(0, 1, 2 * n_elements + 1)

    # Guardamos la solución para la comparación posterior
    solutions["4_elements_quadratic"] = u_quad_4_elements

    # Definir la solución analítica de la ecuación diferencial
    def analytical_solution(x):
        """Solución analítica de la ecuación diferencial"""
        return -x**2 / 2 + 0.5 * (x + np.sin(x)) 

    # Definir la norma L2 para evaluar el error
    def l2_error(u_num, u_analytical, x):
        """Calcula el error en norma L2 entre la solución numérica y la analítica"""
        error = np.sqrt(np.sum((u_num - u_analytical)**2 * np.diff(x, append=x[-1])))
        return error

    # Crear el dominio fino para la solución analítica
    x_fine = np.linspace(0, 1, 1000)
    u_analytical = analytical_solution(x_fine)

    # Graficar las soluciones obtenidas junto con la solución analítica
    plt.figure(figsize=(10, 6))

    # Solución analítica
    plt.plot(x_fine, u_analytical, label="Solución Analítica", linestyle='--', color='black')

    # Solución numérica con 4 elementos lineales
    plt.plot(x, solutions["4_elements_linear"], label="4 elementos, lineal", marker='o')

    # Solución numérica con 8 elementos lineales
    plt.plot(x_8_elements, solutions["8_elements_linear"], label="8 elementos, lineal", marker='x')

    # Solución numérica con 2 elementos cuadráticos
    plt.plot(x_quad_2_elements, solutions["2_elements_quadratic"], label="2 elementos, cuadrático", marker='s')

    # Solución numérica con 4 elementos cuadráticos
    plt.plot(x_quad_4_elements, solutions["4_elements_quadratic"], label="4 elementos, cuadrático", marker='d')

    plt.title("Comparación de soluciones numéricas y analítica")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{RUTA_OUTPUTS}/{TEMATICA}_comparacion{FORMATO_GRAFICAS}")
    plt.savefig(f"{RUTA_OUTPUTS_LATEX}/{TEMATICA}_comparacion{FORMATO_GRAFICAS}")

    # Calcular errores en norma L2 para las diferentes soluciones numéricas
    u_analytical_4 = analytical_solution(x)
    u_analytical_8 = analytical_solution(x_8_elements)
    u_analytical_quad_2 = analytical_solution(x_quad_2_elements)
    u_analytical_quad_4 = analytical_solution(x_quad_4_elements)


    # Calcular errores
    errors = {
        "4 elementos lineal": l2_error(solutions["4_elements_linear"], u_analytical_4, x),
        "8 elementos lineal": l2_error(solutions["8_elements_linear"], u_analytical_8, x_8_elements),
        "2 elementos cuadrático": l2_error(solutions["2_elements_quadratic"], u_analytical_quad_2, x_quad_2_elements),
        "4 elementos cuadrático": l2_error(solutions["4_elements_quadratic"], u_analytical_quad_4, x_quad_4_elements),
    }

    # Mostrar los errores calculados
    error_df = pd.DataFrame(list(errors.items()), columns=["Método", "Error L2"])
    
    
    def convertir_tabla_a_latex(df: pd.DataFrame, ruta_salida: str):
        latex_code = df.to_latex(index=False)
        with open(ruta_salida, 'w') as f:
            f.write(latex_code)
        informer.info(f"Tabla en formato LaTeX guardada en {ruta_salida}")
        
    convertir_tabla_a_latex(error_df, f"{RUTA_OUTPUTS}/{TEMATICA}_tabla_errores.tex")

    # Graficar los errores para un análisis visual de la convergencia
    plt.figure(figsize=(8, 5))
    methods = list(errors.keys())
    error_values = list(errors.values())
    plt.bar(methods, error_values, color=['blue', 'orange', 'green', 'red'])
    plt.title("Errores L2 de los diferentes métodos")
    plt.ylabel("Error L2")
    plt.xticks(rotation=45, ha='right')
    plt.grid(True)
    plt.savefig(f"{RUTA_OUTPUTS}/{TEMATICA}_l2{FORMATO_GRAFICAS}")
    plt.savefig(f"{RUTA_OUTPUTS_LATEX}/{TEMATICA}_l2{FORMATO_GRAFICAS}")
    
    # Función para realizar el cálculo de la solución numérica y el error para diferentes números de elementos
    def compute_error_convergence(n_elements, interpolation_type="linear"):
        if interpolation_type == "linear":
            # Para elementos lineales
            K = stiffness_matrix_linear(n_elements)
            f = load_vector_linear(n_elements)
            K_bc, f_bc = apply_boundary_conditions(K, f)
            u_num = solve(K_bc, f_bc)
            x = np.linspace(0, 1, n_elements + 1)
        elif interpolation_type == "quadratic":
            # Para elementos cuadráticos
            K = stiffness_matrix_quadratic(n_elements)
            f = load_vector_quadratic(n_elements)
            K_bc, f_bc = apply_boundary_conditions(K, f)
            u_num = solve(K_bc, f_bc)
            x = np.linspace(0, 1, 2 * n_elements + 1)
        
        # Solución analítica en los puntos correspondientes
        u_analytical = analytical_solution(x)
        
        # Cálculo del error L2
        error = l2_error(u_num, u_analytical, x)
        
        return error

    # Número de elementos a probar para el estudio de convergencia
    element_counts = [2, 4, 8, 16, 32]
    errors_linear = []
    errors_quadratic = []

    # Calcular errores para interpolaciones lineales y cuadráticas
    for n_elements in element_counts:
        error_linear = compute_error_convergence(n_elements, interpolation_type="linear")
        error_quadratic = compute_error_convergence(n_elements, interpolation_type="quadratic")
        errors_linear.append(error_linear)
        errors_quadratic.append(error_quadratic)

    # Graficar la convergencia en un gráfico log-log
    plt.figure(figsize=(8, 6))
    plt.loglog(element_counts, errors_linear, label="Interpolación lineal", marker='o')
    plt.loglog(element_counts, errors_quadratic, label="Interpolación cuadrática", marker='s')
    plt.xlabel("Número de elementos")
    plt.ylabel("Error L2")
    plt.title("Estudio de convergencia")
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig(f"{RUTA_OUTPUTS}/{TEMATICA}_l2_convergencia{FORMATO_GRAFICAS}")
    plt.savefig(f"{RUTA_OUTPUTS_LATEX}/{TEMATICA}_l2_convergencia{FORMATO_GRAFICAS}")
