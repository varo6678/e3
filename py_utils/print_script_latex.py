import os

def escape_latex(text):
    return text.replace('_', r'\_')

def generate_latex_from_directory(root_dir, output_tex_file, language="python"):

    with open(output_tex_file, 'w') as tex_file:
        # Escribimos la cabecera del archivo LaTeX
        tex_file.write(r"\chapter{e}" + "\n")

        # Función para recorrer la estructura del paquete e
        def process_directory(directory, depth=0):
            for item in sorted(os.listdir(directory)):
                item_path = os.path.join(directory, item)

                # Si es un directorio, crear una nueva subsección
                if os.path.isdir(item_path):
                    tex_file.write(f"\n\\{'sub' * depth}section{{{escape_latex(item)}}}\n")
                    process_directory(item_path, depth + 1)
                
                # Si es un archivo .py, insertar su contenido en un bloque minted
                elif item.endswith('.py'):
                    tex_file.write(f"\n\\{'sub' * (depth + 1)}section{{{escape_latex(item)}}}\n")
                    tex_file.write(r"\begin{minted}[fontsize={\fontsize{5.5}{6.5}\selectfont}, breaklines]{python}" + "\n")
                    
                    # Leer y escribir el contenido del archivo de código
                    with open(item_path, 'r') as code_file:
                        tex_file.write(code_file.read())
                    
                    # Finalizamos el bloque minted con salto de línea
                    tex_file.write("\n" + r"\end{minted}" + "\n\n")
        
        # Procesamos el directorio raíz
        process_directory(root_dir)
    
    print(f"Archivo LaTeX generado: {output_tex_file}")


if __name__ == "__main__":
    # Ruta al paquete 'e'
    root_directory = "e"  # Cambia esto por la ruta a tu paquete

    # Nombre del archivo de salida .tex
    output_tex_file = "/home/ramanujan/e3/e4_latex/anexo/apendiceC.tex"

    # Generar el archivo LaTeX con el contenido del paquete e
    generate_latex_from_directory(root_directory, output_tex_file)
