import os

def print_directory_tree(root_dir, indent=""):
    """
    Imprime la estructura del directorio en formato de árbol.
    
    :param root_dir: Directorio raíz desde donde empezar a imprimir.
    :param indent: Identación usada para el árbol.
    """
    # Listar los archivos y directorios en el directorio actual
    files_and_dirs = os.listdir(root_dir)
    
    # Iterar por todos los elementos listados
    for i, name in enumerate(files_and_dirs):
        path = os.path.join(root_dir, name)
        # Verificar si es el último elemento en el directorio actual
        is_last = (i == len(files_and_dirs) - 1)
        # Seleccionar los caracteres adecuados para la rama del árbol
        tree_char = "└── " if is_last else "├── "
        print(indent + tree_char + name)

        # Si es un directorio, hacer recursión
        if os.path.isdir(path):
            new_indent = indent + ("    " if is_last else "│   ")
            print_directory_tree(path, new_indent)

if __name__ == "__main__":
    # Ruta al paquete 'e'
    root_directory = "e"  # Cambia esto si el paquete está en otro lugar
    print(f"Estructura del paquete '{root_directory}':")
    print_directory_tree(root_directory)
