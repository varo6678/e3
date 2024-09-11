import numpy as np
import matplotlib.pyplot as plt

def solve_laplace_SOR_semicircular_cartesian(N_x, N_y, R, T0, T1, tol=1e-6, max_iter=10000, omega=1.25):
    # Definir los pasos espaciales
    Delta_x = R / (N_x - 1)
    Delta_y = R / (N_y - 1)
    
    # Inicializar la malla con un gradiente lineal
    u = np.linspace(T0, T1, N_x).reshape(N_x, 1) * np.ones((1, N_y))
    
    # Generar las coordenadas x e y
    x = np.linspace(0, R, N_x)
    y = np.linspace(0, R, N_y)
    X, Y = np.meshgrid(x, y)
    
    # Crear una máscara para definir el dominio semicircular
    mask = (X**2 + Y**2 <= R**2)  # Puntos dentro del semicírculo

    # Condiciones de contorno
    u[:, 0] = T0  # En el diámetro (y = 0), T0
    u[-1, mask[-1, :]] = T1  # En el borde circular (x^2 + y^2 = R^2), T1
    
    # Iteraciones del método SOR
    for it in range(max_iter):
        u_old = np.copy(u)
        
        # Actualizar solo los puntos interiores de la malla dentro del semicírculo
        for i in range(1, N_x - 1):
            for j in range(1, N_y - 1):
                if mask[i, j]:  # Solo actualizar puntos dentro del semicírculo
                    u_new = (1 / (2 * (1 / Delta_x**2 + 1 / Delta_y**2))) * (
                        (u[i+1, j] + u[i-1, j]) / Delta_x**2 + 
                        (u[i, j+1] + u[i, j-1]) / Delta_y**2
                    )
                    
                    # Sobrerrelajación
                    u[i, j] = (1 - omega) * u[i, j] + omega * u_new
        
        # Criterio de convergencia
        diff = np.max(np.abs(u - u_old))
        if diff < tol:
            print(f'Convergencia alcanzada en {it} iteraciones.')
            break
    else:
        print('No se alcanzó la convergencia.')

    return u, mask, X, Y

# Parámetros del problema
N_x = 100  # Aumentar el número de puntos en la dirección x
N_y = 100  # Aumentar el número de puntos en la dirección y
R = 1.0    # Radio del semicírculo
T0 = 0.0   # Temperatura en el diámetro (y = 0)
T1 = 50.0  # Temperatura en el borde circular

# Resolver la ecuación de Laplace
solution, mask, X, Y = solve_laplace_SOR_semicircular_cartesian(N_x, N_y, R, T0, T1, omega=1.25, tol=1e-6)

# Aplicar la máscara a la solución (establecer NaN en puntos fuera del semicírculo)
solution_masked = np.where(mask, solution, np.nan)

# Gráfica de la solución en el semicírculo usando pcolormesh
plt.figure(figsize=(6, 5))
plt.pcolormesh(X, Y, solution_masked, cmap='hot', shading='auto')
plt.colorbar(label='Temperatura')
plt.title('Distribución de temperatura en un semicírculo en coordenadas cartesianas')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')  # Asegurar que los ejes tienen la misma escala
plt.show()
