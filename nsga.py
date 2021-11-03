import os

def lectura_archivo(file):
    m1 = m2 = []

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

file1 = os.path.join(os.getcwd(), "instancias", "tsp_KROAB100.TSP.TXT")
file2 = os.path.join(os.getcwd(), "instancias", "tsp_kroac100.tsp.txt") 

#valores = lectura_archivo(file1)
#print(len(valores))



