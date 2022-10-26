
# Programa para graficar valores iterando el archivo exe creado por fortran

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200
#plt.rcParams["figure.figsize"] = (10,6)

import pandas as pd
import os
import time

tf = 13
alpha = 50
betas = [0.33, 0.5, 0.67, 0.91, 1, 1.33, 2]

file = 'Results_data.txt'
data = pd.DataFrame()

#Loop que corre el main.exe creado en fortran con los valores de beta
for beta in betas:
    with open('input_data.txt','w') as f:
        f.write(str(tf)+"\n")
        f.write(str(alpha)+"\n")
        f.write(str(beta))

    time.sleep(0.5)

    os.startfile('main.exe')
    time.sleep(0.5)

    datax = pd.read_csv("Results_Data.txt", header=None, delim_whitespace=True)

    data['x_beta'+str(beta)] = datax[2]
    data['y_beta'+str(beta)] = datax[1]

#Valores guardados en dataframe son graficados
for beta in betas:
    x = 'x_beta'+str(beta)
    y = 'y_beta'+str(beta)
    plt.plot(data[x],data[y], ls='-',linewidth=1.5, label=r'$\beta=$'+str(beta))

plt.grid()
plt.xlabel(r"$\tau'$")
plt.ylabel(r'$\mathcal{A}$')
plt.xlim([0,8])
#plt.text(6.05,310,r'$\alpha=$'+str(alpha))
plt.legend()
plt.show();