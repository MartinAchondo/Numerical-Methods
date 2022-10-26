#Shell Script for compile and run f90 code
name="main"

s_command="gfortran "${name}".f90 -o "${name}
$s_command

echo "Compiled !!"
echo "----------------"

./main.exe
 