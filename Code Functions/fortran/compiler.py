import os

name = 'main'

command = f'gfortran {name}.f90 -o {name}'
os.system(command)
print('Compiled !!')