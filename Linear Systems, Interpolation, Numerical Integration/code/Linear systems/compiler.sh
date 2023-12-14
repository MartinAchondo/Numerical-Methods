#Shell Script for compile and run f90 code
name="p1"

s_command="gfortran "${name}".f90 -o "${name}" -L./lib/ -llapack -lblas"
$s_command

echo "Compiled !!"
echo "----------------"

s2_command="./"${name}".exe"
$s2_command
 