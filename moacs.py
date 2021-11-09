import os
import pandas as pd
import numpy as np
from nsga import lectura_archivo


# Inicializar parametros
max_iter = 20
beta = 0.10         # Beta
m = 10              # Número de Hormigas
ciudades = []       # Las ciudades
q0 = 0.6            # Probabilidad de exploracion o explotacion
tau0 = 0         # Tau subcero
rho = 0.02

file1 = os.path.join(os.getcwd(), "instancias", "tsp_KROAB100.TSP.TXT")
file2 = os.path.join(os.getcwd(), "instancias", "tsp_kroac100.tsp.txt")
valores = lectura_archivo(file1)

costo1 = [[float(i) for i in j] for j in valores[2][:100]]                                 # Costo 1 entre las ciudades
costo2 = [[float(i) for i in j] for j in valores[2][100:]]                                 # Costo 2 entre las ciudades

n = len(costo1)
eta1 = [[0 for i in range(n)] for j in range(n)]    # Visibilidad de los arcos para el objetivo 1
eta2 = [[0 for i in range(n)] for j in range(n)]    # Visibilidad de los arcos para el objetivo 2
tau = [[0 for i in range(n)] for j in range(n)]     # Feromonas de los arcos

ciudades = [i for i in range(n)]
for ciudad in ciudades:
    for vecino in ciudades[ciudad+1:]:
        eta1[ciudad][vecino] = eta1[vecino][ciudad] = 1 / costo1[ciudad][vecino]
        eta2[ciudad][vecino] = eta2[vecino][ciudad] = 1 / costo2[ciudad][vecino]
        tau[ciudad][vecino] = 0


class Colonia:
    def __init__(self, tam=m):
        """Clase que representa la colonia de hormigas. tam es el numero de hormigas.
        Las hormigas son inicializadas aleatoriamente."""
        self.tam = tam
        self.miembros = [Hormiga(i+1) for i in range(tam)]
        self.pareto_set = self.miembros
        for miembro in self.miembros:
            miembro.construir_camino_sin_feromonas()
        self.actualizar_frente_pareto()
        self.actualizacion_global_feromonas()


    def encontrar_no_dominados(self, hormigas):
        """Encontrar las hormigas no dominadas en la iteracion actual."""
        pareto_set = []
        for hormiga1 in hormigas:
            hormiga1_dominada = False
            for hormiga2 in hormigas:
                if not hormiga1_dominada and hormiga2 != hormiga1:
                    if hormiga1.costos_camino[0] >= hormiga2.costos_camino[0] and hormiga1.costos_camino[1] >= hormiga2.costos_camino[1] and (hormiga1.costos_camino[0] > hormiga2.costos_camino[0] or hormiga1.costos_camino[1] > hormiga2.costos_camino[1]):
                        hormiga1_dominada = True
            if not hormiga1_dominada:
                pareto_set.append(hormiga1)
        if len(pareto_set) > 0:
            return pareto_set

    def actualizar_frente_pareto(self):
        no_dominados = self.encontrar_no_dominados(self.miembros)
        hormigas = set(no_dominados).union(set(self.pareto_set))
        self.pareto_set = copy.deepcopy(self.encontrar_no_dominados(list(hormigas)))

    def nuevos_caminos(self):
        """Envia a todas las hormigas de la colonia en busca de un nuevo camino."""
        for miembro in self.miembros:
            miembro.construir_nuevo_camino()

    def actualizacion_global_feromonas(self):
        """Actualizacion global de feromonas. Debe ser llamado despues de que se haya encontrado a la mejor hormiga."""


        """ AL TERMINAR LA ITERACION SE ACTUALIZA LOS TAU CON LA ECUACION GENERICA"""
        #  usar pareto set y encontrar el frente pareto
        #  if tau_0 < tau_0':
        #   reiniciar tau, osea la tabla de feromonas
        #  else:
        #   for elem in elementos en el conjunto pareto:
        #       i = elem.i
        #       j = elem.j
        #       tau(i, j) = (1 - ro)*tau(i, j) + ro / (evaluan_funcion_1(elemento_pareto) * evaluan_funcion_2(elemento_pareto))

        hormigas_optimas = len(self.pareto_set)
        promedio_costo_1 = sum([hormiga.costos_camino[0] for hormiga in self.pareto_set]) / hormigas_optimas
        promedio_costo_2 = sum([hormiga.costos_camino[1] for hormiga in self.pareto_set]) / hormigas_optimas
        tau0_prima = 1 / (promedio_costo_1 * promedio_costo_2)
        global tau0
        if tau0 < tau0_prima:
            global tau
            tau0 = tau0_prima
            tau = [[tau0_prima if i != j else 0 for i in range(n)] for j in range(n)]
        else:
            for hormiga in self.pareto_set:
                for idx in range(1, len(hormiga.camino)):
                    i = hormiga.camino[idx-1]
                    j = hormiga.camino[idx]
                    tau[i][j] = tau[j][i] = (1 - rho) * tau[i][j] + rho / (costo1[i][j] * costo2[i][j])

        # DeltaBest = 1 / self.best.path_length
        # # common update to all edges
        # for edge in pheromones.keys():
        #     pheromones[edge] *= (1 - rho)
        # # update to edges in best path
        # for i in range(len(self.best.path) - 1):
        #     if (self.best.path[i], self.best.path[i + 1]) in pheromones.keys():
        #         pheromones[(self.best.path[i], self.best.path[i + 1])] += rho * DeltaBest
        #     elif (self.best.path[i + 1], self.best.path[i]) in pheromones.keys():
        #         pheromones[(self.best.path[i + 1], self.best.path[i])] += rho * DeltaBest
        #     else:
        #         raise IndexError('Something went wrong.')
        # # edges that closes the cycle of the best path
        # if (self.best.path[0], self.best.path[-1]) in pheromones.keys():
        #     pheromones[(self.best.path[0], self.best.path[-1])] += rho * DeltaBest
        # elif (self.best.path[-1], self.best.path[0]) in pheromones.keys():
        #     pheromones[(self.best.path[-1], self.best.path[0])] += rho * DeltaBest
        # else:
        #     raise IndexError('Something went wrong.')


class Hormiga:
    def __init__(self, num):
        """Class that represents an ant that follows a path with a certain length."""
        self.num = num
        self.camino = []
        self.costos_camino = [[], []]

    def construir_camino_sin_feromonas(self):
        """An ant builds a new path according to pseudorandom proportional rule.
        Pheromone values are updated locally. The total distance of the path is
        computed."""
        self.camino = []  # erase previous path
        for i in range(len(self.costos_camino)):
            self.costos_camino[i] = 0

        lambd = self.num / m
        self.camino.append(0)  # start from BCN
        for _ in range(1, n):
            # select possible next stops
            pos_sig = [i for i in range(n)]
            # delete already visited cities
            for visited in self.camino:
                pos_sig.remove(visited)
            q = os.urandom(1)[0] / 255
            i = self.camino[-1]
            if q <= q0:
                max_prob = (-1, -1)
                for j in pos_sig:
                    prob_actual = (eta1[i][j] ** (lambd * beta)) * (eta2[i][j] ** ((1 - lambd) * beta))
                    if max_prob[1] < prob_actual:
                        max_prob = (j, prob_actual)
                sig_ciudad = max_prob[0]
            else:
                suma_vecinos = 0
                tau_temporal = {}
                for j in pos_sig:
                    t = (eta1[i][j] ** (lambd * beta)) * (eta2[i][j] ** ((1 - lambd) * beta))
                    tau_temporal.update({j: t})
                    suma_vecinos += t
                probabilidades = [tau_temporal[camino]/suma_vecinos for camino in tau_temporal]
                sig_ciudad = np.random.choice(list(tau_temporal.keys()), p=probabilidades)

            self.camino.append(sig_ciudad)
            self.costos_camino[0] += costo1[i][j]
            self.costos_camino[1] += costo2[i][j]
            tau[self.camino[-2]][self.camino[-1]] = tau[self.camino[-1]][self.camino[-2]] = (1 - rho) * tau[self.camino[-2]][self.camino[-1]] + rho * tau0

        tau[0][self.camino[-1]] = tau[self.camino[-1]][0] = (1 - rho) * tau[self.camino[-2]][self.camino[-1]] + rho * tau0
        self.camino.append(0)

    def construir_nuevo_camino(self):
        """An ant builds a new path according to pseudorandom proportional rule.
        Pheromone values are updated locally. The total distance of the path is
        computed."""
        self.camino = []  # erase previous path
        for i in range(len(self.costos_camino)):
            self.costos_camino[i] = 0

        lambd = self.num / m
        self.camino.append(0)  # start from BCN
        for _ in range(1, n):
            # select possible next stops
            pos_sig = [i for i in range(n)]
            # delete already visited cities
            for visited in self.camino:
                pos_sig.remove(visited)
            # probs = {}

            """ calculo de tau y eta de la hormiga """
            # Estaticos por paso
            # promedio_func_1 = promedio_pareto(1)
            # promedio_func_2 = promedio_pareto(2)
            # tau_0 = 1 / (promedio_func_1 * promedio_func_2)
            # lambda = k / m
            q = os.urandom(1)[0] / 255
            i = self.camino[-1]
            if q <= q0:
                max_prob = (-1, -1)
                for j in pos_sig:
                    prob_actual = tau[i][j] * (eta1[i][j] ** (lambd * beta)) * (eta2[i][j] ** ((1 - lambd) * beta))
                    if max_prob[1] < prob_actual:
                        max_prob = (j, prob_actual)
                sig_ciudad = max_prob[0]
            else:
                suma_vecinos = 0
                tau_temporal = {}
                for j in pos_sig:
                    t = tau[i][j] * (eta1[i][j] ** (lambd * beta)) * (eta2[i][j] ** ((1 - lambd) * beta))
                    tau_temporal.update({j: t})
                    suma_vecinos += t
                probabilidades = [tau_temporal[camino]/suma_vecinos for camino in tau_temporal]
                sig_ciudad = np.random.choice(list(tau_temporal.keys()), p=probabilidades)

            self.camino.append(sig_ciudad)
            self.costos_camino[0] += costo1[i][j]
            self.costos_camino[1] += costo2[i][j]

            # Se modifica por paso
            #


            # get tau and eta for each next possible option
            # for opt in pos_sig:
            #     if (self.camino[-1], opt) in costo1.keys():
            #         probs.update({opt: pheromones[(self.camino[-1], opt)] ** alpha /
            #                            europe[(self.camino[-1], opt)] ** beta})
            #     elif (opt, self.camino[-1]) in europe.keys():
            #         probs.update({opt: pheromones[(opt, self.camino[-1])] ** alpha /
            #                            europe[(opt, self.camino[-1])] ** beta})
            #     else:
            #         raise IndexError('Something went wrong.')
            # total = sum(list(probs.values()))
            # for k in probs:
            #     probs[k] = probs[k] / total
            # # select next destination: pseudorandom proportional rule
            # q = os.urandom(1)[0] / 255.
            # if q <= q0:
            #     next_index = np.argmax(list(probs.values()))
            #     self.camino.append(list(probs.keys())[next_index])
            # else:
            #     next_city = np.random.choice(list(probs.keys()),
            #                                  p=list(probs.values()))
            #     self.camino.append(next_city)

            tau[self.camino[-2]][self.camino[-1]] = tau[self.camino[-1]][self.camino[-2]] = (1 - rho) * tau[self.camino[-2]][self.camino[-1]] + rho * tau0

            # local pheromone update and add edge cost to total distance
            # if (self.camino[-1], self.camino[-2]) in pheromones.keys():
            #     pheromones[(self.camino[-1], self.camino[-2])] = (1 - phi) * pheromones[
            #         (self.camino[-1], self.camino[-2])] + phi * tau0
            #     self.camino_length += europe[(self.camino[-1], self.camino[-2])]
            # elif (self.camino[-2], self.camino[-1]) in pheromones.keys():
            #     pheromones[(self.camino[-2], self.camino[-1])] = (1 - phi) * pheromones[
            #         (self.camino[-2], self.camino[-1])] + phi * tau0
            #     self.camino_length += europe[(self.camino[-2], self.camino[-1])]
            # else:
            #     raise IndexError('Something went wrong.')


        tau[0][self.camino[-1]] = tau[self.camino[-1]][0] = (1 - rho) * tau[self.camino[-2]][self.camino[-1]] + rho * tau0
        self.camino.append(0)

        # after exit the loop: local pheromone and distance update on edge closing the cycle
        # if (self.camino[-1], 'Barcelona') in pheromones.keys():
        #     pheromones[(self.camino[-1], 'Barcelona')] = (1 - phi) * pheromones[
        #         (self.camino[-1], 'Barcelona')] + phi * tau0
        #     self.camino_length += europe[(self.camino[-1], 'Barcelona')]
        # elif ('Barcelona', self.camino[-1]) in pheromones.keys():
        #     pheromones[('Barcelona', self.camino[-1])] = (1 - phi) * pheromones[
        #         ('Barcelona', self.camino[-1])] + phi * tau0
        #     self.camino_length += europe[('Barcelona', self.camino[-1])]
        # else:
        #     raise IndexError('Something went wrong.')

"""Ant Colony System"""
import copy
col = Colonia()
bestSoFar = copy.deepcopy(col.pareto_set)
iters = 1
stop = 250
bestIter = 0

while iters <= stop:
    print(iters)
    col.nuevos_caminos()
    col.actualizar_frente_pareto()
    col.actualizacion_global_feromonas()
    bestSoFar = copy.deepcopy(col.pareto_set)
    bestIter = iters
    iters += 1

# output result
print('Best solution found during iteration #' + str(bestIter) + '.')
print('This clever ant visits all cities in a ' + str(bestSoFar.path_length) +
      ' km tour.')
print('Her adventorous journey is:', end = ' ')
for city in bestSoFar.camino:
  print(city, end = '-')
print(bestSoFar.camino[0])

# Inicializar
# while
#     for cada hormiga
#         contruir una solucion
#         algo = contruir_tour()
#         if x esta en conjunto pareto
#             guardar x y borrar las soluciones dominadas del conjunto pareto
#     calcular theta sub-cero prima con la formula
#     if theta sub-cero prima es mejor a theta sub-cero
#         significa que un mejor conjuto pareto se encontro
#         asignar theta sub-cero prima a theta sub-cero
#         se reinicializa theta con theta sub-cero
#     else
#         for x in conjunto pareto
#             se hace una actualizacion global
#             theta(i,j) = (1 - p) * theta(i,j) + p / algunas weas