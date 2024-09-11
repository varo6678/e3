from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import numpy as np

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
tol = 1e-6 # Tolerancia

# Resolver para una malla gruesa y una más densa
solution_coarse = solve_laplace_SOR(N_r=50, N_theta=50, T0=T0, T1=T1, tol=tol)
solution_fine = solve_laplace_SOR(N_r=100, N_theta=100, T0=T0, T1=T1, tol=tol)

# Interpolación de coarse a fine
interp_coarse_to_fine = RegularGridInterpolator((np.linspace(0, 1, 50), np.linspace(0, np.pi, 50)), solution_coarse)
r_fine = np.linspace(0, 1, 100)
theta_fine = np.linspace(0, np.pi, 100)
R_fine, Theta_fine = np.meshgrid(r_fine, theta_fine)
points_fine = np.array([R_fine.flatten(), Theta_fine.flatten()]).T
fine_interpolated_from_coarse = interp_coarse_to_fine(points_fine).reshape(100, 100)

# Interpolación inversa de fine a coarse
interp_fine_to_coarse = RegularGridInterpolator((np.linspace(0, 1, 100), np.linspace(0, np.pi, 100)), solution_fine)
r_coarse = np.linspace(0, 1, 50)
theta_coarse = np.linspace(0, np.pi, 50)
R_coarse, Theta_coarse = np.meshgrid(r_coarse, theta_coarse)
points_coarse = np.array([R_coarse.flatten(), Theta_coarse.flatten()]).T
coarse_interpolated_from_fine = interp_fine_to_coarse(points_coarse).reshape(50, 50)

# Verificar las diferencias
diff_fine_to_coarse = np.abs(solution_coarse - coarse_interpolated_from_fine)
diff_coarse_to_fine = np.abs(solution_fine - fine_interpolated_from_coarse)

print(f"Diferencia máxima entre la solución coarse y la interpolación desde fine: {np.max(diff_fine_to_coarse)}")
print(f"Diferencia máxima entre la solución fine y la interpolación desde coarse: {np.max(diff_coarse_to_fine)}")


# Visualización de las diferencias excluyendo los bordes
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.contourf(R_coarse[2:-2, 2:-2], Theta_coarse[2:-2, 2:-2], np.abs(solution_coarse[2:-2, 2:-2] - coarse_interpolated_from_fine[2:-2, 2:-2]), 50, cmap='coolwarm')
plt.colorbar(label='Diferencia')
plt.title('Diferencia interna coarse vs interpolada desde fine')

plt.subplot(1, 2, 2)
plt.contourf(R_fine[2:-2, 2:-2], Theta_fine[2:-2, 2:-2], np.abs(solution_fine[2:-2, 2:-2] - fine_interpolated_from_coarse[2:-2, 2:-2]), 50, cmap='coolwarm')
plt.colorbar(label='Diferencia')
plt.title('Diferencia interna fine vs interpolada desde coarse')

plt.show()

# Comparar coarse con la interpolación desde fine
# Ajustar la tolerancia si las diferencias son mínimas
# Excluir las esquinas de las comparaciones
# Excluir esquinas y bordes de la comparación
# assert np.allclose(solution_coarse[2:-2, 2:-2], coarse_interpolated_from_fine[2:-2, 2:-2], atol=1e-2), "La interpolación desde fine no coincide con la malla coarse (excluyendo esquinas y bordes)."
# assert np.allclose(solution_fine[2:-2, 2:-2], fine_interpolated_from_coarse[2:-2, 2:-2], atol=1e-2), "La interpolación desde coarse no coincide con la malla fine (excluyendo esquinas y bordes)."

# assert np.allclose(solution_coarse[1:-1, 1:-1], coarse_interpolated_from_fine[1:-1, 1:-1], atol=1e-2), "La interpolación desde fine no coincide con la malla coarse (excluyendo esquinas)."
# assert np.allclose(solution_fine[1:-1, 1:-1], fine_interpolated_from_coarse[1:-1, 1:-1], atol=1e-2), "La interpolación desde coarse no coincide con la malla fine (excluyendo esquinas)."

# # Verificar las esquinas por separado
# assert solution_coarse[0, 0] == coarse_interpolated_from_fine[0, 0], "La esquina (0,0) no coincide"
# assert solution_coarse[0, -1] == coarse_interpolated_from_fine[0, -1], "La esquina (0,pi) no coincide"
# assert solution_coarse[-1, 0] == coarse_interpolated_from_fine[-1, 0], "La esquina (1,0) no coincide"
# assert solution_coarse[-1, -1] == coarse_interpolated_from_fine[-1, -1], "La esquina (1,pi) no coincide"

# # Verificar las zonas cercanas a los bordes con una tolerancia más amplia
# assert np.allclose(solution_coarse[1:-1, :], coarse_interpolated_from_fine[1:-1, :], atol=1e-1), "La interpolación desde fine no coincide en los bordes de la malla coarse."
# assert np.allclose(solution_fine[1:-1, :], fine_interpolated_from_coarse[1:-1, :], atol=1e-1), "La interpolación desde coarse no coincide en los bordes de la malla fine."

# Verificar las zonas internas
assert np.allclose(solution_coarse[2:-2, 2:-2], coarse_interpolated_from_fine[2:-2, 2:-2], atol=1e-2), "La interpolación desde fine no coincide con la malla coarse (excluyendo esquinas y bordes)."
assert np.allclose(solution_fine[2:-2, 2:-2], fine_interpolated_from_coarse[2:-2, 2:-2], atol=1e-2), "La interpolación desde coarse no coincide con la malla fine (excluyendo esquinas y bordes)."

# Verificar los bordes con mayor tolerancia
assert np.allclose(solution_coarse[1:-1, :], coarse_interpolated_from_fine[1:-1, :], atol=1e-1), "La interpolación desde fine no coincide en los bordes de la malla coarse."
assert np.allclose(solution_fine[1:-1, :], fine_interpolated_from_coarse[1:-1, :], atol=1e-1), "La interpolación desde coarse no coincide en los bordes de la malla fine."



# Verificar que el residuo es pequeño en la solución fina
def check_residual(u, N_r, N_theta, Delta_r, Delta_theta):
    residual = np.zeros_like(u)
    for i in range(1, N_r-1):
        r_i = i * Delta_r
        for j in range(1, N_theta-1):
            residual[i, j] = (
                (u[i+1, j] - 2*u[i, j] + u[i-1, j]) / Delta_r**2 +
                (1 / r_i) * (u[i+1, j] - u[i-1, j]) / (2 * Delta_r) +
                (1 / r_i**2) * (u[i, j+1] - 2*u[i, j] + u[i, j-1]) / Delta_theta**2
            )
    return residual

# Calcular el residuo para la malla fina
residual = check_residual(solution_fine, N_r=100, N_theta=100, Delta_r=1.0/99, Delta_theta=np.pi/99)
assert np.max(np.abs(residual)) < tol, f"El residuo es demasiado grande: {np.max(np.abs(residual))}"
