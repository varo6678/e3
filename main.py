# main.py

# imports | externos.
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

def solve_laplace_SOR(N_r, N_theta, T0, T1, tol=1e-6, max_iter=10000, omega=1.5):
    # Definir los pasos espaciales
    Delta_r = 1.0 / (N_r - 1)
    Delta_theta = np.pi / (N_theta - 1)
    
    # Inicializar la malla
    u = np.zeros((N_r, N_theta))
    
    # Condiciones de frontera
    u[-1, :] = T1  # En r = 1 (Dirichlet)
    u[:, 0] = T0   # En theta = 0 (Dirichlet)
    u[:, -1] = T0  # En theta = pi (Dirichlet)

    # Corregir las esquinas: Prioridad a la condición de r = 1
    u[-1, 0] = T1   # Esquina en (r = 1, theta = 0)
    u[-1, -1] = T1  # Esquina en (r = 1, theta = pi)

    # Corregir las esquinas en r = 0 para theta = 0 y theta = pi
    u[0, 0] = T0  # Esquina en (r = 0, theta = 0)
    u[0, -1] = T0  # Esquina en (r = 0, theta = pi)
    
    # Iteraciones del método SOR
    for it in range(max_iter):
        u_old = np.copy(u)
        
        # Actualizar solo los puntos interiores de la malla
        for i in range(1, N_r - 1):  # Excluir las condiciones en r = 0 y r = 1
            r_i = i * Delta_r
            for j in range(1, N_theta - 1):  # Excluir las condiciones en theta = 0 y theta = pi
                u_new = (1 / (2 * (1/Delta_r**2 + 1/(r_i**2 * Delta_theta**2)))) * (
                    (u[i+1, j] + u[i-1, j]) / Delta_r**2 + 
                    (1 / (r_i**2)) * (u[i, j+1] + u[i, j-1]) / Delta_theta**2
                )
                
                # Sobrerrelajación
                u[i, j] = (1 - omega) * u[i, j] + omega * u_new
        
        # Condición de Neumann en r = 0 (simetría)
        u[0, :] = u[1, :]

        # Reaplicar las condiciones de Dirichlet en theta = 0 y theta = pi
        u[:, 0] = T0  # En theta = 0
        u[:, -1] = T0  # En theta = pi

        # Corregir nuevamente las esquinas explícitamente en cada iteración
        u[0, 0] = T0   # Esquina en (r = 0, theta = 0)
        u[0, -1] = T0  # Esquina en (r = 0, theta = pi)
        u[-1, 0] = T1  # Esquina en (r = 1, theta = 0)
        u[-1, -1] = T1 # Esquina en (r = 1, theta = pi)
        
        # Criterio de convergencia
        diff = np.max(np.abs(u - u_old))
        if diff < tol:
            print(f'Convergencia alcanzada en {it} iteraciones.')
            break
    else:
        print('No se alcanzó la convergencia.')

    return u


# Parámetros del problema
N_r = 50
N_theta = 50
T0 = 0.0   # Dirichlet en theta = 0 y theta = pi
T1 = 50.0  # Dirichlet en r = 1
tol = 1e-6 # Definir la tolerancia global p ara ser usada luego

# # Resolver la ecuación
# solution = solve_laplace_SOR(N_r, N_theta, T0, T1, tol=tol)
# print("Valores en r = 1:", solution[-1, :])
# print("Valores en rtheta = 0:", solution[:, 0])
# print("Valores en theta = pi", solution[:, -1])

# # Verificación de las condiciones de frontera excluyendo las esquinas
# assert np.allclose(solution[-1, 1:-1], T1, atol=tol), "La condición de Dirichlet en r = 1 no se cumple (excluyendo las esquinas)."  # Verificar fila final excluyendo las esquinas
# assert np.allclose(solution[1:-1, 0], T0, atol=tol), "La condición de Dirichlet en theta = 0 no se cumple (excluyendo las esquinas)."  # Verificar toda la columna en theta = 0, excluyendo esquinas
# assert np.allclose(solution[1:-1, -1], T0, atol=tol), "La condición de Dirichlet en theta = pi no se cumple (excluyendo las esquinas)."  # Verificar toda la columna en theta = pi, excluyendo esquinas

# # Verificación de las esquinas
# assert solution[0, 0] == T0, "La esquina (r = 0, theta = 0) no respeta la condición de Dirichlet."
# assert solution[0, -1] == T0, "La esquina (r = 0, theta = pi) no respeta la condición de Dirichlet."
# assert solution[-1, 0] == T1, "La esquina (r = 1, theta = 0) no respeta la condición de Dirichlet."
# assert solution[-1, -1] == T1, "La esquina (r = 1, theta = pi) no respeta la condición de Dirichlet."



#####################################################################################
# Mejorar el apartado de graficas personalizandolas.
#####################################################################################

# # Gráfica de la solución
# r = np.linspace(0, 1, N_r)
# theta = np.linspace(0, np.pi, N_theta)
# R, Theta = np.meshgrid(r, theta)

# # Convertir a coordenadas cartesianas para visualización
# X = R * np.cos(Theta)
# Y = R * np.sin(Theta)

# plt.figure(figsize=(6, 5))
# plt.contourf(X, Y, solution.T, 50, cmap='hot')
# plt.colorbar(label='Temperatura [°C]')
# plt.title('Distribución de temperatura en una lámina semicircular')
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.show()

#####################################################################################


# # Probar con diferentes valores de omega
# for omega in [1.0, 1.25, 1.5, 1.75, 1.9]:
#     print(f"\nResolviendo con omega = {omega}")
#     solution = solve_laplace_SOR(N_r, N_theta, T0, T1, omega=omega, tol=tol)

# # Resolver para una malla más densa
# solution_coarse = solve_laplace_SOR(N_r=50, N_theta=50, T0=0, T1=50, tol=tol)
# solution_fine = solve_laplace_SOR(N_r=100, N_theta=100, T0=0, T1=50, tol=tol)

# from scipy.interpolate import RegularGridInterpolator

# def check_residual(u, N_r, N_theta, Delta_r, Delta_theta):
#     residual = np.zeros_like(u)
#     for i in range(1, N_r-1):
#         r_i = i * Delta_r
#         for j in range(1, N_theta-1):
#             residual[i, j] = (
#                 (u[i+1, j] - 2*u[i, j] + u[i-1, j]) / Delta_r**2 +
#                 (1 / r_i) * (u[i+1, j] - u[i-1, j]) / (2 * Delta_r) +
#                 (1 / r_i**2) * (u[i, j+1] - 2*u[i, j] + u[i, j-1]) / Delta_theta**2
#             )
#     return residual

# # Comparar la solución interpolada con la solución fina en los puntos de la malla fina
# interp_coarse_solution = RegularGridInterpolator((np.linspace(0, 1, 50), np.linspace(0, np.pi, 50)), solution_coarse)

# # Generar los puntos en la malla fina
# r_fine = np.linspace(0, 1, 100)
# theta_fine = np.linspace(0, np.pi, 100)
# R_fine, Theta_fine = np.meshgrid(r_fine, theta_fine)
# points_fine = np.array([R_fine.flatten(), Theta_fine.flatten()]).T

# # Interpolar en los puntos de la malla fina
# fine_interpolated_solution = interp_coarse_solution(points_fine).reshape(100, 100)

# # Verificar que la interpolación de coarse a fine es consistente
# assert np.allclose(solution_fine, fine_interpolated_solution, atol=1e-2), "La solución no es consistente en los puntos de la malla fina."

# # Visualización
# plt.figure(figsize=(12, 5))

# plt.subplot(1, 3, 1)
# plt.contourf(R_fine, Theta_fine, solution_coarse, 50, cmap='hot')
# plt.title('Solución coarse (50x50)')

# plt.subplot(1, 3, 2)
# plt.contourf(R_fine, Theta_fine, fine_interpolated_solution, 50, cmap='hot')
# plt.title('Solución interpolada')

# plt.subplot(1, 3, 3)
# plt.contourf(R_fine, Theta_fine, solution_coarse - fine_interpolated_solution, 50, cmap='coolwarm')
# plt.title('Diferencia (coarse - interpolada)')
# plt.colorbar()

# plt.show()

# # Verificar que el residuo es pequeño
# residual = check_residual(solution_fine, N_r=100, N_theta=100, Delta_r=1.0/99, Delta_theta=np.pi/99)
# assert np.max(np.abs(residual)) < tol, f"El residuo es demasiado grande: {np.max(np.abs(residual))}"
