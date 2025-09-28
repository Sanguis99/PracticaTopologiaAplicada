from itertools import combinations # Para crear las caras dados los vertices

# Clase de los simplices
class Simplice:
    def __init__(self, vertices):
        self.vertices = vertices
        self.caras = self.calcular_caras()
        self.dimension = len(vertices) - 1

    def calcular_caras(self):
        caras = set()
        n = len(self.vertices)
        for k in range(1, n + 1):
            #Para calcular las caras del símplice, vemos todas las posibles combinaciones que se pueden
            #formar con los vértices, para ello se utiliza el paquete combinations
            for cara in combinations(self.vertices, k):
                caras.add(tuple(cara))
        return caras

# Clase de los complejos simpliciales
class Complejo_simplicial:
    def __init__(self, simplices): #simplices es un set de simplices
        self.simplices = simplices
        self.c = self.calcular_caras()
        #la dimensión del complejo simplicial es la dimensión máxima de los símplices
        self.d = max(s.dimension for s in simplices) if simplices else 0

    # Definimos las caras del complejo simplicial usando las caras de los símplices maximales
    def calcular_caras(self):
        caras = set()
        #Añadimos las caras de cada símplice. Las caras de cada símplice ya las calculamos en la clase Simplice
        for s in self.simplices:
            for cara in s.caras:
                caras.add(cara)
        return sorted(caras, key=lambda x: x) #lambda expression, recibe x y lo devuelve sin cambiarlo

    # Este metodo permite extraer las caras de dimensión n
    def n_caras(self, n):
        if n < 0 or n > self.d:
            print(f"No hay caras de dimensión {n} en el complejo.")
            return []
        else:
            #Miramos en nuestro atributo c (en el que se almacenan todas las caras del complejo simplicial) si tienen
            #dimension n y la añadimos
            caras_n = sorted(set([cara for cara in self.c if len(cara) == n+1]), key=lambda x: x)
            print(f"Caras de dimensión {n}: {caras_n}")
            return caras_n
    
    # Los siguientes métodos son para poder imprimir las caras y la dimensión del complejo
    def caras(self):
        print(f"Caras del complejo: {self.c}")

    def dimension(self):
        print(f"Dimensión del complejo: {self.d}")

    ################################# CLASE 2 #######################################
    # Cálculo de la característica de Euler
    def caras_por_dimension(self):
        caras_dim = [self.n_caras(i) for i in range(self.d + 1)]
        # No hace falta imprimir aquí porque n_caras ya imprime
        return caras_dim

    def Euler(self):
        chi = 0
        # Definimos el sumatorio para la característica de Euler
        for i in range(self.d + 1):
            # Calculamos el número de caras de cada dimensión, las de dimensión para se suman
            # y las de dimensión impar se restan, obteniendo así la característica de Euler
            chi += (-1) ** i * len(self.caras_por_dimension()[i])
        print(f"Característica de Euler: {chi}")
        return chi

    #La estrella de un símplice c es el conjunto de todas las cocaras de c
    def estrella(self, c):
        # Todas las caras que contienen a c
        estrella = set([cara for cara in self.c if set(c).issubset(set(cara))])
        estrella = sorted(estrella, key=lambda x: x)
        print(f"Estrella de {c}: {estrella}")
        return estrella

    #La estrella cerrada de c es el menor subcomplejo de K que contiene a la estrella de c.
    def estrella_cerrada(self, c):
        # Encuentra todas las caras que contienen al menos un vértice de c
        caras_con_v = [cara for cara in self.c if any(v in cara for v in c)]
        # Añade todas las subcaras de esas caras
        estrella_cerrada = set()
        for cara in caras_con_v:
            for k in range(1, len(cara)+1):
                for subcara in combinations(cara, k):
                    estrella_cerrada.add(tuple(sorted(subcara)))
        estrella_cerrada = sorted(estrella_cerrada, key=lambda x: x)
        print(f"Estrella cerrada de {c}: {estrella_cerrada}")
        return estrella_cerrada

    #El link de un símplice c es el conjunto de todos los símplices de la estrella cerrada de c
    #cuya intersección con la estrella de c es vacía
    def link(self, c):
        estrella_cerrada = self.estrella_cerrada(c)
        estrella = self.estrella(c)
        link = [cara for cara in estrella_cerrada if cara not in estrella]
        print(f"Link de {c}: {link}")
        return link

    def j_esqueleto(self, j):
        # Comprobamos que j es válido
        if j < 0 or j > self.d:
            print(f"No hay esqueleto de dimensión {j} en el complejo.")
            return []
        else:
            #Añadimos aquellas caras que tengan una longitud menor o igual a j+1
            esqueleto = sorted(set([cara for cara in self.c if len(cara) <= j + 1]), key=lambda x: x)
            print(f"Esqueleto de dimensión {j}: {esqueleto}")
            return esqueleto

    def componentes_conexas(self):
        visited = set()
        components = []

        # bep significa Búsqueda en Profundidad
        def bep(v, component):
            visited.add(v)
            component.append(v)
            for cara in self.c:
                if v in cara:
                    for u in cara:
                        if u not in visited:
                            bep(u, component)

        for cara in self.c:
            for v in cara:
                if v not in visited:
                    component = []
                    bep(v, component)
                    components.append(sorted(component))
        print(f"Componentes conexas: {components}")
        return components

    #Calculamos el número de componentes conexas
    def connected_components(self):
        componentes = self.componentes_conexas()
        return len(componentes)

    # El complejo será conexo si tiene una única componente conexa
    def es_conexo(self):
        if self.connected_components() == 1:
            print("El complejo es conexo.")
            return True
        else:
            print("El complejo no es conexo.")
            return False


# Ejemplo de uso
if __name__ == "__main__":
    # Definimos los símplices maximales del complejo
    # Por manejo más sencillo, los vértices son enteros
    s1 = Simplice([0,1,2])
    s2 = Simplice([2,3])
    s3 = Simplice([4])
    # Creamos el complejo simplicial
    complejo = Complejo_simplicial([s1, s2,s3])
    # Probamos los métodos
    print("####################################################")
    print("#         Ejercicio Complejos Simpliciales         #")
    print("####################################################")
    complejo.caras()
    complejo.caras_por_dimension()
    complejo.dimension()
    complejo.Euler()
    complejo.estrella((2,))
    complejo.estrella_cerrada((2,))
    complejo.link((2,))
    complejo.componentes_conexas()
    complejo.j_esqueleto(1)