from itertools import combinations # Para crear las caras dados los vertices

# Información sobre los headers de las funciones:
# Las funciones xx_aux() se usan para calcular xx y devolver el resultado.
# Las funciones xx() se usan para imprimir el resultado de xx_aux().
# A excepción de caras_por_dimension() que usa n_caras() como auxiliar.
# De esta forma, si se quiere usar el resultado de xx en otro código,
# se puede usar xx_aux() sin que imprima nada por pantalla

###################################### CLASE 1 ######################################
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
            # Para calcular las caras del símplice, vemos todas las posibles combinaciones que se pueden
            # formar con los vértices, para ello se utiliza el paquete combinations
            for cara in combinations(self.vertices, k):
                caras.add(tuple(cara))
        return caras

# Clase de los complejos simpliciales
class Complejo_simplicial:
    def __init__(self, simplices):
        self.simplices = set(simplices)
        self.c = self.calcular_caras()
        # La dimensión del complejo simplicial es la dimensión máxima de los símplices
        self.d = max(s.dimension for s in simplices) if simplices else 0

    # Definimos las caras del complejo simplicial usando las caras de los símplices maximales
    def calcular_caras(self):
        caras = set()
        # Añadimos las caras de cada símplice. Las caras de cada símplice ya las calculamos en la clase Simplice
        for s in self.simplices:
            for cara in s.caras:
                caras.add(cara)
        return sorted(caras, key=lambda x: x) # lambda expression que ordena las caras por su valor inicial

    # Este metodo permite extraer las caras de dimensión n
    def n_caras(self, n):
        if n < 0 or n > self.d:
            print(f"No hay caras de dimensión {n} en el complejo.")
            return []
        else:
            # Miramos en nuestro atributo c (Caras del complejo simplicial) si tienen dimensión n y la añadimos
            caras_n = sorted(set([cara for cara in self.c if len(cara) == n+1]), key=lambda x: x)
            return caras_n
    
    # Los siguientes métodos son para poder imprimir las caras y la dimensión del complejo
    def caras(self):
        print(f"Caras del complejo: {self.c}")

    def dimension(self):
        print(f"Dimensión del complejo: {self.d}")

###################################### CLASE 2 ######################################
    # Calculamos el número de caras por dimensión
    def caras_por_dimension(self):
        caras_dim = [self.n_caras(i) for i in range(self.d + 1)]
        for i in range(self.d + 1):
            print(f"Caras de dimensión {i}: {caras_dim[i]}")
        return caras_dim

    # Cálculo de la característica de Euler
    def Euler(self):
        chi = 0
        # Definimos el sumatorio para la característica de Euler
        for i in range(self.d + 1):
            # Calculamos el número de caras de cada dimensión, las de dimensión para se suman
            # y las de dimensión impar se restan, obteniendo así la característica de Euler
            chi += (-1) ** i * len(self.n_caras(i))
        print(f"Característica de Euler: {chi}")
        return chi

    # La estrella de un símplice c es el conjunto de todas las cocaras de c
    def estrella_aux(self, c):
        estrella = set([cara for cara in self.c if set(c).issubset(set(cara))])
        estrella = sorted(estrella, key=lambda x: x)
        return estrella
    # Usamos la función auxiliar para calcular la estrella
    # y luego la imprimimos
    def estrella(self, c):
        # Todas las caras que contienen a c
        estrella = self.estrella_aux(c)
        print(f"Estrella de {c}: {estrella}")
        return estrella

    # La estrella cerrada de c es el menor subcomplejo de K que contiene a la estrella de c.
    def estrella_cerrada_aux(self, c):
        # Encuentra todas las caras que contienen al menos un vértice de c
        caras_con_v = [cara for cara in self.c if any(v in cara for v in c)]
        # Añade todas las subcaras de esas caras
        estrella_cerrada = set()
        for cara in caras_con_v:
            for k in range(1, len(cara)+1):
                for subcara in combinations(cara, k):
                    estrella_cerrada.add(tuple(sorted(subcara)))
        estrella_cerrada = sorted(estrella_cerrada, key=lambda x: x)
        return estrella_cerrada
    # Usamos la función auxiliar para calcular la estrella cerrada
    # y luego la imprimimos
    def estrella_cerrada(self, c):
        estrella_cerrada = self.estrella_cerrada_aux(c)
        print(f"Estrella cerrada de {c}: {estrella_cerrada}")
        return estrella_cerrada

    # El link de un símplice c es el conjunto de todos los símplices de la estrella cerrada de c
    # cuya intersección con la estrella de c es vacía
    def link_aux(self, c):
        estrella_cerrada = self.estrella_cerrada_aux(c)
        estrella = self.estrella_aux(c)
        link = [cara for cara in estrella_cerrada if cara not in estrella]
        return link
    # Usamos la función auxiliar para calcular el link
    # y luego la imprimimos
    def link(self, c):
        link = self.link_aux(c)
        print(f"Link de {c}: {link}")
        return link

    def j_esqueleto_aux(self, j):
        # Comprobamos que j es válido
        if j < 0 or j > self.d:
            print(f"No hay esqueleto de dimensión {j} en el complejo.")
            return []
        else:
            # Añadimos aquellas caras que tengan una longitud menor o igual a j+1
            esqueleto = sorted(set([cara for cara in self.c if len(cara) <= j + 1]), key=lambda x: x)
            return esqueleto
    # Usamos la función para calcular el j-esqueleto
    # y luego la imprimimos
    def j_esqueleto(self, j):
        esqueleto = self.j_esqueleto_aux(j)
        print(f"{j}-esqueleto del complejo: {esqueleto}")
        return esqueleto

    # Se calculan las componentes conexas del complejo usando búsqueda en profundidad (BEP)
    def componentes_conexas_aux(self):
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
        return components
    # Usamos la función auxiliar para calcular las componentes conexas
    # y luego las imprimimos
    def componentes_conexas(self):
        components = self.componentes_conexas_aux()
        print(f"Componentes conexas del complejo: {components}")
        return components

    # Calculamos el número de componentes conexas
    def connected_components(self):
        return len(self.componentes_conexas_aux())

    # El complejo será conexo si tiene una única componente conexa
    def es_conexo(self):
        if self.connected_components() == 1:
            print("El complejo es conexo.")
            return True
        else:
            print("El complejo no es conexo.")
            return False
        
    def insert(self, simplices):
        for s in simplices:
            # Evitamos añadir símplices repetidos
            if any(set(s.vertices) == set(existing.vertices) for existing in self.simplices):
                continue
            self.simplices.add(s)
        self.c = self.calcular_caras()
        self.d = max(s.dimension for s in self.simplices) if self.simplices else 0

###################################### CLASE 3 ######################################
class Simplice_filtrado(Simplice):
    def __init__(self, vertices, index):
        super().__init__(vertices)
        self.index = float(index)

class Complejo_simplicial_filtrado(Complejo_simplicial):
    def __init__(self, simplices_filtrados):
        # Comprobamos que todos los elementos son de tipo Simplice_filtrado
        for s in simplices_filtrados:
            if not isinstance(s, Simplice_filtrado):
                raise ValueError("Todos los elementos deben ser de tipo Simplice_filtrado")
        super().__init__(simplices_filtrados)
        # Ordenamos los símplices primero por índice de filtrado y luego por dimensión
        self.update_simplices_ordenados()

    def update_simplices_ordenados(self):
        self.simplices_ordenados = sorted(self.simplices, key=lambda x: (x.index, x.dimension))

    # Insertar un conjunto de símplices con el mismo índice de filtrado
    def insert_filtrado(self, simplices, index):
        for s in simplices:
            s1 = Simplice_filtrado(s.vertices, index)
            # Si ya existe un símplice con los mismos vértices, mantenemos el de menor índice
            if any(set(s1.vertices) == set(existing.vertices) for existing in self.simplices):
                e = [existing for existing in self.simplices if set(s1.vertices) == set(existing.vertices)][0]
                if s1.index < e.index:
                    self.simplices.remove(e)
                else:
                    continue
            self.simplices.add(s1)
        self.c = self.calcular_caras()
        self.d = max(s.dimension for s in self.simplices) if self.simplices else 0
        self.update_simplices_ordenados()


    def simplices_por_filtrado_aux(self, index):
        sf = sorted([s for s in self.simplices if s.index <= index], key=lambda x: (x.index, x.dimension))
        return sf
    # Usamos la función auxiliar para calcular los símplices con índice de filtrado menor o igual a index
    # y luego los imprimimos
    def simplices_por_filtrado(self, index):
        sf = self.simplices_por_filtrado_aux(index)
        print(f"Símplices con índice de filtrado menor o igual a {index}: {[ (s.vertices, s.index) for s in sf ]}")
        return sf


###################################### Ejemplo de Uso ######################################
if __name__ == "__main__":
    # Definimos los símplices maximales del complejo
    # Por manejo más sencillo, los vértices son enteros
    s1 = Simplice([0,1,2])
    s2 = Simplice([2,3])
    s3 = Simplice([4])
    # Creamos el complejo simplicial
    complejo = Complejo_simplicial([s1, s2, s3])
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
    complejo.connected_components()
    complejo.es_conexo()
    print("####################################################")
    print("#    Ejercicio Complejos Simpliciales Filtrados    #")
    print("####################################################")
    csf = Complejo_simplicial_filtrado([])
    csf.insert_filtrado([s1, s2], 0)
    csf.insert_filtrado([s3], 1)
    csf.insert_filtrado([s1], 0.5)
    csf.caras()
    csf.caras_por_dimension()
    csf.simplices_por_filtrado(0)
    csf.simplices_por_filtrado(1)
    print(f"Simplices ordenados: {[(s.vertices, s.index) for s in csf.simplices_ordenados]}")