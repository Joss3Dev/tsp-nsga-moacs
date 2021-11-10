import os
import pandas as pd
import numpy as np
import math
from nsga import lectura_archivo, Individuo
import matplotlib.pyplot as plt



# Inicializar parametros
max_iter = 20
beta = 0.10         # Beta
m = 20              # Número de Hormigas
ciudades = []       # Las ciudades
q0 = 0.6            # Probabilidad de exploracion o explotacion
tau0 = 0.1*10**-20         # Tau subcero
rho = 0.02
n = 0

# file1 = os.path.join(os.getcwd(), "instancias", "tsp_KROAB100.TSP.TXT")
# file2 = os.path.join(os.getcwd(), "instancias", "tsp_kroac100.tsp.txt")
# valores = lectura_archivo(file1)
#
# costo1 = [[float(i) for i in j] for j in valores[2][:100]]                                 # Costo 1 entre las ciudades
# costo2 = [[float(i) for i in j] for j in valores[2][100:]]                                 # Costo 2 entre las ciudades
#
# n = len(costo1)
# eta1 = [[0 for i in range(n)] for j in range(n)]    # Visibilidad de los arcos para el objetivo 1
# eta2 = [[0 for i in range(n)] for j in range(n)]    # Visibilidad de los arcos para el objetivo 2
# tau = [[0 for i in range(n)] for j in range(n)]     # Feromonas de los arcos
#
# ciudades = [i for i in range(n)]
# for ciudad in ciudades:
#     for vecino in ciudades[ciudad+1:]:
#         eta1[ciudad][vecino] = eta1[vecino][ciudad] = 1 / costo1[ciudad][vecino]
#         eta2[ciudad][vecino] = eta2[vecino][ciudad] = 1 / costo2[ciudad][vecino]
#         tau[ciudad][vecino] = 0


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

    def draw(self, pareto_set, subplot=111):
        """
        Dibuja el frente pareto.

        @param subplot: posición del gráfico.
        """
        fig = plt.figure()
        pf_ax = fig.add_subplot(subplot)
        pf_ax.set_title(u"Frente Pareto")
        for p in pareto_set:
            pf_ax.scatter(p.costos_camino[0], p.costos_camino[1], marker='o', facecolor='blue')
        plt.show()

    def distancia(self, a, b):
        dist = 0
        for i in range(len(a.costos_camino)):
            dist += (a.costos_camino[i]-b.costos_camino[i]) ** 2
        return math.sqrt(dist)

    def m1(self, y_true, frente_patero):
        suma = 0
        for p in frente_patero:
            dist = []
            for y in y_true:
                dist.append(self.distancia(p, y))
            suma = suma + min(dist)
        return suma / len(frente_patero)

    def m2(self, sigma, frente_pareto):
        suma = 0
        for h1 in frente_pareto:
            for h2 in frente_pareto:
                if h1 != h2 and self.distancia(h1,h2) > sigma:
                    suma = suma + 1
        if len(frente_pareto) - 1 > 0:
            resultado = suma / (len(frente_pareto) - 1)
        else:
            resultado = suma
        return resultado

    def m3(self, frente_pareto):
        dist_x = []
        dist_y = []
        n = len(frente_pareto)
        pareto = frente_pareto
        for i in range(n - 1):
            for j in range(i + 1, n):
                dist_x.append((pareto[i].costos_camino[0] - pareto[j].costos_camino[0]) ** 2)
                dist_y.append((pareto[i].costos_camino[1] - pareto[j].costos_camino[1]) ** 2)
        return math.sqrt(max(dist_x) + max(dist_y))

    def error(self, frente_pareto, y_true):
        interseccion = []
        for i in frente_pareto:
            agregado = False
            for j in y_true:
                if i.costos_camino[0] == j.costos_camino[0] and i.costos_camino[1] == j.costos_camino[1] and i.camino == j.camino and not agregado:
                    interseccion.append(i)
                    agregado = True
        return 1 - len(interseccion) / len(frente_pareto)


class Hormiga:
    def __init__(self, num):
        self.num = num
        self.camino = []
        self.costos_camino = [[], []]

    def construir_camino_sin_feromonas(self):
        self.camino = []  # erase previous path
        for i in range(len(self.costos_camino)):
            self.costos_camino[i] = 0

        lambd = self.num / m
        ini = int(os.urandom(1)[0] / 255 * 100)
        self.camino.append(0)
        for _ in range(1, n):
            pos_sig = [i for i in range(n)]
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
        self.camino = []
        for i in range(len(self.costos_camino)):
            self.costos_camino[i] = 0
        lambd = self.num / m
        ini = int(os.urandom(1)[0] / 255 * 100)
        self.camino.append(0)
        for _ in range(1, n):
            pos_sig = [i for i in range(n)]
            for visited in self.camino:
                pos_sig.remove(visited)
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
            tau[self.camino[-2]][self.camino[-1]] = tau[self.camino[-1]][self.camino[-2]] = (1 - rho) * tau[self.camino[-2]][self.camino[-1]] + rho * tau0

        tau[0][self.camino[-1]] = tau[self.camino[-1]][0] = (1 - rho) * tau[self.camino[-2]][self.camino[-1]] + rho * tau0
        self.camino.append(0)

print("KROAB100.TSP")
import copy


def test_metricas_moacs(file, y_true_nsga):
    valores = file

    global costo1, costo2, n, eta1, eta2, tau, ciudades

    costo1 = [[float(i) for i in j] for j in valores[2][:100]]  # Costo 1 entre las ciudades
    costo2 = [[float(i) for i in j] for j in valores[3][:100]]  # Costo 2 entre las ciudades

    n = len(costo1)
    eta1 = [[0 for i in range(n)] for j in range(n)]  # Visibilidad de los arcos para el objetivo 1
    eta2 = [[0 for i in range(n)] for j in range(n)]  # Visibilidad de los arcos para el objetivo 2
    tau = [[tau0 for i in range(n)] for j in range(n)]  # Feromonas de los arcos

    ciudades = [i for i in range(n)]
    for ciudad in ciudades:
        for vecino in ciudades[ciudad + 1:]:
            eta1[ciudad][vecino] = eta1[vecino][ciudad] = 1 / costo1[ciudad][vecino]
            eta2[ciudad][vecino] = eta2[vecino][ciudad] = 1 / costo2[ciudad][vecino]
            tau[ciudad][vecino] = tau[vecino][ciudad] = tau0

    frentes_pareto = []
    y_true_calculada = y_true_nsga
    for _ in range(5):
        col = Colonia()
        bestSoFar = copy.deepcopy(col.pareto_set)
        iters = 1
        stop = 500
        bestIter = 0

        while iters <= stop:
            # print(iters)
            col.nuevos_caminos()
            col.actualizar_frente_pareto()
            col.actualizacion_global_feromonas()
            conjunto_pareto = copy.deepcopy(col.pareto_set)
            iters += 1
        frentes_pareto.append(copy.deepcopy(conjunto_pareto))
        y_true_calculada = set(y_true_calculada).union(set(conjunto_pareto))
        y_true_calculada = col.encontrar_no_dominados(y_true_calculada)
        col.draw(col.pareto_set)
    col.draw(y_true_calculada)

    return y_true_calculada, frentes_pareto, col


def evaluar_moacs(y_true_calculada, frentes_pareto, col):
    m1 = []
    m2 = []
    m3 = []
    error = []
    for frente in frentes_pareto:
        m1.append(col.m1(y_true_calculada, frente))
        m2.append(col.m2(5000, frente))
        m3.append(col.m3(frente))
        error.append(col.error(frente, y_true_calculada))

    m1_prom = sum(m1)/len(m1)
    m2_prom = sum(m2)/len(m2)
    m3_prom = sum(m3)/len(m3)
    error_prom = sum(error)/len(error)

    print('M1 promedio: ', m1_prom)
    print('M2 promedio: ', m2_prom)
    print('M3 promedio: ', m3_prom)
    print('Error promedio: ', error_prom)

# test_metricas_moacs(file2)

# print()
# print()
# print()
# print("tsp_kroac100.TSP")
# valores = lectura_archivo(file2)
#
# costo1 = [[float(i) for i in j] for j in valores[2][:100]]                                 # Costo 1 entre las ciudades
# costo2 = [[float(i) for i in j] for j in valores[2][100:]]                                 # Costo 2 entre las ciudades
#
# n = len(costo1)
# eta1 = [[0 for i in range(n)] for j in range(n)]    # Visibilidad de los arcos para el objetivo 1
# eta2 = [[0 for i in range(n)] for j in range(n)]    # Visibilidad de los arcos para el objetivo 2
# tau = [[0 for i in range(n)] for j in range(n)]     # Feromonas de los arcos
#
# ciudades = [i for i in range(n)]
# for ciudad in ciudades:
#     for vecino in ciudades[ciudad+1:]:
#         eta1[ciudad][vecino] = eta1[vecino][ciudad] = 1 / costo1[ciudad][vecino]
#         eta2[ciudad][vecino] = eta2[vecino][ciudad] = 1 / costo2[ciudad][vecino]
#         tau[ciudad][vecino] = 0
#
# frentes_pareto = []
# y_true_calculada = []
# m1 = []
# m2 = []
# m3 = []
# error = []
# for _ in range(5):
#     col = Colonia()
#     bestSoFar = copy.deepcopy(col.pareto_set)
#     iters = 1
#     stop = 500
#     bestIter = 0
#
#     while iters <= stop:
#         # print(iters)
#         col.nuevos_caminos()
#         col.actualizar_frente_pareto()
#         col.actualizacion_global_feromonas()
#         conjunto_pareto = copy.deepcopy(col.pareto_set)
#         iters += 1
#     frentes_pareto.append(copy.deepcopy(conjunto_pareto))
#     y_true_calculada = set(y_true_calculada).union(set(conjunto_pareto))
#     y_true_calculada = col.encontrar_no_dominados(y_true_calculada)
#     col.draw(col.pareto_set)
# col.draw(y_true_calculada)
#
# for frente in frentes_pareto:
#     m1.append(col.m1(y_true_calculada, frente))
#     m2.append(col.m2(5000, frente))
#     m3.append(col.m3(frente))
#     error.append(col.error(frente, y_true_calculada))
#
# m1_prom = sum(m1)/len(m1)
# m2_prom = sum(m2)/len(m2)
# m3_prom = sum(m3)/len(m3)
# error_prom = sum(error)/len(error)
#
# print('M1 promedio: ', m1_prom)
# print('M2 promedio: ', m2_prom)
# print('M3 promedio: ', m3_prom)
# print('Error promedio: ', error_prom)


def y_true_nsga_a_moacs(y_true):
    y_true_moacs = []
    i = 1
    for ind in y_true:
        y_true_moacs.append(individuo_a_hormiga(ind, i))
        i += 1
    return y_true_moacs


def individuo_a_hormiga(individuo, i):
    h = Hormiga(i)
    h.camino = individuo.cromosoma
    h.costos_camino[0] = individuo.f1
    h.costos_camino[1] = individuo.f2
    return h


def y_true_moacs_a_nsga(y_true):
    y_true_nsga = []
    for h in y_true:
        y_true_nsga.append(hormiga_a_individuo(h))
    return y_true_nsga


def hormiga_a_individuo(h):
    i = Individuo(h.camino, costo1, costo2)
    i.cromosoma = h.camino
    i.f1 = h.costos_camino[0]
    i.f2 = h.costos_camino[1]
    return i
