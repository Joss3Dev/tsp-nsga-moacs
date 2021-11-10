import os, random, math, statistics

class Individuo():
    cromosoma = []
    f1 = 0
    f2 = 0
    ciudad_origen = None
    fitness = None

    def __init__(self, cromosoma, m1, m2):
        self.cromosoma = cromosoma
        self.ciudad_origen = self.cromosoma[0] 
        self.__evaluacion_f1(m1)
        self.__evaluacion_f2(m2)

    def __evaluacion_f1(self, m1):
        for i in range(len(self.cromosoma)):
            if i == len(self.cromosoma) -1 :
                self.f1 += float(m1[self.cromosoma[i]][self.ciudad_origen])
            else:
                self.f1 += float(m1[self.cromosoma[i]][self.cromosoma[i+1]]) 

    def __evaluacion_f2(self, m2):
        for i in range(len(self.cromosoma)):
            if i == len(self.cromosoma) -1 :
                self.f2 += float(m2[self.cromosoma[i]][self.ciudad_origen])
            else:
                self.f2 += float(m2[self.cromosoma[i]][self.cromosoma[i+1]]) 


class NSGA():
    size_poblacion = None 
    poblacion = []
    max_gen = None 
    sigma_share = None
    dummy_fitness = None
    porcentaje_elitismo = None
    porcentaje_mutacion = None 
    frentes_pareto = []

    def __init__(self, size_poblacion=50, max_gen=500, sigma_share=1.7, dummy_fitness=100, porcentaje_elitismo=0.08, porcentaje_mutacion=0.06):
        self.size_poblacion = size_poblacion
        self.max_gen = max_gen
        self.sigma_share = sigma_share
        self.dummy_fitness = dummy_fitness
        self.porcentaje_elitismo = porcentaje_elitismo
        self.porcentaje_mutacion = porcentaje_mutacion 
    
    def crear_poblacion_inicial(self, datos):
        for _ in range(self.size_poblacion):
            c = random.sample(range(int(datos[0])), int(datos[0]))
            self.poblacion.append(Individuo(c, tupla[2], tupla[3]))
    
    def iniciar_algoritmo(self):
        for _ in range(self.max_gen):
            rank_frentes = self.__ranking_de_frentes()
            nueva_pobla = self.__operadores_geneticos(rank_frentes)
            self.poblacion = nueva_pobla[:]
        return self.frentes_pareto[-1]
    
    def __ranking_de_frentes(self):
        rank = []
        clasificacion = self.poblacion[:]
        k = 1

        while len(clasificacion) > 0:
            dominancia = []
            frente = []

            for p in clasificacion:
                n = 0
                for q in clasificacion:
                    if (p.f1 <= q.f1 and p.f2 <= q.f2) and (p.f1 < q.f1 or p.f2 < q.f2) and q not in dominancia:
                        dominancia.append(q)
                    elif (q.f1 <= p.f1 and q.f2 <= p.f2) and (q.f1 < p.f1 or q.f2 < p.f2):
                        n += 1

                if n == 0 and p not in frente:
                    frente.append(p)
            
            if k == 1:
                #El mejor frente encontrado se guarda
                self.frentes_pareto.append(frente[:])

            t = (float(self.dummy_fitness / k), frente)
            f = self.__sharing_function(t)
            rank.append(f)
            k += 1 
            clasificacion = dominancia 
        return rank
    
    def __sharing_function(self,pareto):
        N = len(pareto[1])
        current_front = pareto[1]
        for i in range(N):
            niche_count = 0  
            for j in range(N):
                d = math.sqrt((i - j)**2)
                if d < self.sigma_share:
                    sh = 1 - (float((d / self.sigma_share)) ** 2) 
                else:
                    sh = 0 
                niche_count += sh 
            current_front[i].fitness = float(pareto[0] / niche_count) 
        current_front.sort(reverse=True, key=lambda x: x.fitness)
        return current_front 
    
    def __operadores_geneticos(self,rank):
        all_ind = []
        for r in rank:
            all_ind += r

        elitismo = int(self.size_poblacion * self.porcentaje_elitismo)
        nueva_poblacion = all_ind[:elitismo]
        #El resto del porcentaje restante se toma para el crossover
        crossover = all_ind[elitismo:]
        corte = random.randint(1,80)

        for i in range(0,len(crossover),2):
            nueva_poblacion.append(self.__crossover_un_punto(crossover[i].cromosoma[:], crossover[i+1].cromosoma[:], corte))
            nueva_poblacion.append(self.__crossover_un_punto(crossover[i+1].cromosoma[:], crossover[i].cromosoma[:], corte))

        cant_mut = int(self.size_poblacion * self.porcentaje_mutacion)
        rand = random.sample(range(50), cant_mut)
        for i in rand:
            nueva_poblacion[i] = self.__mutacion(nueva_poblacion[i]) 
        return nueva_poblacion 
    
    def __crossover_un_punto(self,cro1, cro2, pcorte):
        for index in range(pcorte): 
            if cro2[index] not in cro1:
                cro1[index] = cro2[index] 
            else:
                for i in range(len(cro1)):
                    if cro1[i] == cro2[index]:
                        aux = cro1[i]
                        cro1[i] = cro1[index]
                        cro1[index] = aux 
                        break 
        new_ind = Individuo(cro1, tupla[2], tupla[3])
        return new_ind
    
    def __mutacion(self,ind):
        r = random.sample(range(100),2)
        r.sort()
        cro = ind.cromosoma
        ini = cro[:r[0]]
        fin = cro[r[1]:]
        inversion = cro[r[0]:r[1]]
        inversion.reverse()
        new_cromosoma = ini + inversion + fin 
        new_ind = Individuo(new_cromosoma, tupla[2], tupla[3])
        return new_ind



def lectura_archivo(file):
    m1 = []
    m2 = []

    with open(file, "r") as f:
        line = f.readlines() 

    for x in range(1,204):
        if x == 1:
            cant_ciudades = line[x - 1].split()[0]
        elif x == 2:
            cant_objetivos = line[x -1].split()[0]
        elif x < 103:
            m1.append(line[x-1].split()) 
        elif x > 103:
            m2.append(line[x-1].split()) 
    
    return (cant_ciudades, cant_objetivos, m1, m2)


def test_metricas(tupla):
    m1 = []
    m2 = []
    m3 = []
    error = [] 
    conjunto_frentes = []

    for i in range(5):
        ngsa = NSGA()
        ngsa.crear_poblacion_inicial(tupla)
        conjunto_frentes.append(ngsa.iniciar_algoritmo()) 
    
    y_true = no_dominancia(conjunto_frentes)
    for frente in conjunto_frentes:
        m1.append(calculo_m1(y_true, frente)))
        m2.append(calculo_m2(frente))
        m3.append(calculo_m3(frente))
        error.append(calculo_error(y_true, frente)) 
    
    print("Distancia al frente Pareto: ", statistics.mean(m1))
    print("Distribución del frente Pareto: ", statistics.mean(m2))
    print("Extensión del frente Pareto: ", statistics.mean(m3))
    print("Error: ", statistics.mean(error))


def no_dominancia(fronts):
    union_frentes = []
    for i in fronts:
        union_frentes += i

    dominancia = []
    best_front = []

    for p in union_frentes:
        n = 0
        for q in union_frentes:
            if (p.f1 <= q.f1 and p.f2 <= q.f2) and (p.f1 < q.f1 or p.f2 < q.f2) and q not in dominancia:
                dominancia.append(q)
            elif (q.f1 <= p.f1 and q.f2 <= p.f2) and (q.f1 < p.f1 or q.f2 < p.f2):
                n +=    
        if n == 0 and p not in best_front:
            best_front.append(p)
    
    return best_front 


def calculo_m1(y_true, front):   
    suma = 0
    for p in front:
        dist = []
        for q in y_true:
            dist.append(distancia_euclidiana(p,q))
        suma += min(dist)
    return suma / len(front)


def calculo_m2(front):
    sigma_metric = 5000
    suma = 0 
    for p in front:
        for q in front:
            if p != q and distancia_euclidiana(p,q) > sigma_metric:
                suma += 1
    
    if len(front) - 1 > 0:
        resultado = suma / (len(front) - 1) 
    else:
        resultado = 1 
    return resultado


def calculo_m3(front):
    suma = 0
    #Dos objetivos
    for in range(2):
        dist = []
        for p in front:
            for q in front:
                if p != q:
                    dist.append(distancia_euclidiana(p,q))
        suma += max(dist)
    return math.sqrt(suma) 


def calculo_error(y_true, front):
    interseccion = set(y_true) & set(front)
    resultado = 1 -(len(interseccion) / len(front))
    return resultado


def distancia_euclidiana(p, q):
    s = ((p.f1 - q.f1)**2) + ((p.f2 - q.f2)**2)
    dist = math.sqrt(s)
    return dist 



#Primera instancia 
file1 = os.path.join(os.getcwd(), "instancias", "tsp_KROAB100.TSP.TXT")
tupla = lectura_archivo(file1);
print("\nPrimera Instancia\n")
test_metricas(tupla)


#Segunda instancia 
file2 = os.path.join(os.getcwd(), "instancias", "tsp_kroac100.tsp.txt") 
tupla = lectura_archivo(file2);
print("\nSecunda Instancia\n")
test_metricas(tupla)