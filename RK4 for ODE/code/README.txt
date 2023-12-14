
En la carpeta se presentan los archivos:

Código de fortran y ejecutable:
- main.f90
- main.exe

Archivos de entrada y salida:
- input_data.txt
- Results_data.txt

Código que crea gráficos (python) iterando el ejecutable main.exe de fortran
- plots.py

Se incluye este último dado que fue interesante ver como complementar estos dos lenguajes, 
fortran para los cálculos y python para las gráficas. Se encontró la forma de desde python modificar
el archivo de entrada txt y así correr el main.exe. De esta manera se puede correr para todos los betas
en un loop.