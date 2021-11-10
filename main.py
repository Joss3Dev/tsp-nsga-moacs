import os, random, math, statistics
from nsga import lectura_archivo, test_metricas, evaluar
from moacs import test_metricas_moacs, y_true_nsga_a_moacs, y_true_moacs_a_nsga,evaluar_moacs

#Primera instancia
file1 = os.path.join(os.getcwd(), "instancias", "tsp_KROAB100.TSP.TXT")
tupla = lectura_archivo(file1)
print("\nPrimera Instancia\n")
print('NSGA')
y_true, conjunto_frentes = test_metricas(tupla)
print('MOACS')
y_true, frentes_pareto, col = test_metricas_moacs(tupla, y_true_nsga_a_moacs(y_true))
print('Evaluacion')
print('NSGA')
evaluar(y_true_moacs_a_nsga(y_true), conjunto_frentes)
print('MOACS')
evaluar_moacs(y_true, frentes_pareto, col)

#Segunda instancia
file2 = os.path.join(os.getcwd(), "instancias", "tsp_kroac100.tsp.txt")
tupla = lectura_archivo(file2)
print("\nSecunda Instancia\n")
print('NSGA')
y_true, conjunto_frentes = test_metricas(tupla)
print('MOACS')
y_true, frentes_pareto, col = test_metricas_moacs(tupla, y_true_nsga_a_moacs(y_true))
print('Evaluacion')
print('NSGA')
evaluar(y_true_moacs_a_nsga(y_true), conjunto_frentes)
print('MOACS')
evaluar_moacs(y_true, frentes_pareto, col)
