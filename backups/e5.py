import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x):
    B = 1 / np.cos(1)
    u_exact = B * np.sin(1 - x) - x
    return u_exact

def fem_galerkin_1d_linear(n_elements):
    n_nodes = n_elements + 1
    x = np.linspace(0, 1, n_nodes)
    
    # Inicializar matriz de rigidez y vector de cargas
    K = np.zeros((n_nodes, n_nodes))
    f = np.zeros(n_nodes)
    
    # Puntos y pesos de cuadratura de Gauss (2 puntos)
    xi_gauss = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    w_gauss = np.array([1, 1])
    
    for e in range(n_elements):
        nodes = [e, e + 1]
        x_e = x[nodes]
        h_e = x_e[1] - x_e[0]
        
        K_e = np.zeros((2, 2))
        f_e = np.zeros(2)
        
        for i in range(len(xi_gauss)):
            xi = xi_gauss[i]
            w = w_gauss[i]
            
            N = np.array([0.5 * (1 - xi), 0.5 * (1 + xi)])
            dN_dxi = np.array([-0.5, 0.5])
            
            dx_dxi = h_e / 2
            dN_dx = dN_dxi / dx_dxi
            
            x_gauss = N @ x_e
            
            K_e += np.outer(dN_dx, dN_dx) * dx_dxi * w  # Corregido
            f_e += -N * x_gauss * dx_dxi * w
        
        for i_local, i_global in enumerate(nodes):
            for j_local, j_global in enumerate(nodes):
                K[i_global, j_global] += K_e[i_local, j_local]
            f[i_global] += f_e[i_local]
    
    K = K[1:, 1:]
    f = f[1:]
    
    u_interior = np.linalg.solve(K, f)
    u = np.concatenate(([0], u_interior))
    
    return x, u

def fem_galerkin_1d_quadratic(n_elements):
    n_nodes = 2 * n_elements + 1
    x = np.linspace(0, 1, n_nodes)
    
    K = np.zeros((n_nodes, n_nodes))
    f = np.zeros(n_nodes)
    
    xi_gauss = np.array([-np.sqrt(3/5), 0.0, np.sqrt(3/5)])
    w_gauss = np.array([5/9, 8/9, 5/9])
    
    for e in range(n_elements):
        nodes = [2 * e, 2 * e + 1, 2 * e + 2]
        x_e = x[nodes]
        
        K_e = np.zeros((3, 3))
        f_e = np.zeros(3)
        
        for i in range(len(xi_gauss)):
            xi = xi_gauss[i]
            w = w_gauss[i]
            
            N = np.array([
                0.5 * xi * (xi - 1),
                (1 - xi**2),
                0.5 * xi * (xi + 1)
            ])
            
            dN_dxi = np.array([
                xi - 0.5,
                -2 * xi,
                xi + 0.5
            ])
            
            x_gauss = N @ x_e
            dx_dxi = dN_dxi @ x_e  # Corregido el cálculo del Jacobiano
            dN_dx = dN_dxi / dx_dxi
            
            K_e += np.outer(dN_dx, dN_dx) * dx_dxi * w
            f_e += -N * x_gauss * dx_dxi * w
        
        for i_local, i_global in enumerate(nodes):
            for j_local, j_global in enumerate(nodes):
                K[i_global, j_global] += K_e[i_local, j_local]
            f[i_global] += f_e[i_local]
    
    K = K[1:, 1:]
    f = f[1:]
    
    u_interior = np.linalg.solve(K, f)
    u = np.concatenate(([0], u_interior))
    
    return x, u


if __name__ == "__main__":
    # Casos a resolver
    cases = [
        {'n_elements': 4, 'order': 'linear', 'label': '4 elementos lineales'},
        {'n_elements': 8, 'order': 'linear', 'label': '8 elementos lineales'},
        {'n_elements': 2, 'order': 'quadratic', 'label': '2 elementos cuadráticos'},
        {'n_elements': 4, 'order': 'quadratic', 'label': '4 elementos cuadráticos'},
    ]
    
    x_plot = np.linspace(0, 1, 200)
    u_exact = exact_solution(x_plot)
    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, u_exact, 'k-', label='Solución exacta')
    
    errors = []
    elements = []
    
    for case in cases:
        if case['order'] == 'linear':
            x, u = fem_galerkin_1d_linear(case['n_elements'])
        elif case['order'] == 'quadratic':
            x, u = fem_galerkin_1d_quadratic(case['n_elements'])
        else:
            raise ValueError("Orden no soportado")
        
        # Calcular el error
        u_interp = np.interp(x_plot, x, u)
        error = np.abs(u_exact - u_interp)
        max_error = np.max(error)
        print(f"{case['label']}: error máximo = {max_error:.2e}")
        
        plt.plot(x, u, marker='o', label=case['label'])
        
        errors.append(max_error)
        elements.append(len(x))
    
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Solución de la Ecuación Diferencial por el Método de Galerkin')
    plt.legend()
    plt.grid(True)
    plt.savefig('galerkin_solution.png')
    plt.show()
    
    # Estudio de convergencia
    plt.figure()
    plt.loglog(elements, errors, 'o-')
    plt.xlabel('Número de nodos')
    plt.ylabel('Error máximo')
    plt.title('Estudio de convergencia')
    plt.grid(True, which="both", ls="--")
    plt.savefig('convergence_study.png')
    plt.show()
