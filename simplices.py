from itertools import combinations # Para crear las caras dados los vertices
import numpy as np
#Bibliotecas usadas en los ejemplos de Voronoi y Delaunay
from scipy.spatial import Delaunay,Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib.colors

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
        vert = []
        for v in vertices:
            try:
                vert.append(int(v))
            except Exception:
                vert.append(v)
        self.vertices = vert
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

###################################### Clase 6 ######################################
    def matriz_borde_aux(self, p):
        if p < 0 or p > self.d:
            return None
        caras_p = self.n_caras(p)
        if p == 0:
            return np.zeros((1, len(caras_p)), dtype=int)
        caras_p_minus_1 = self.n_caras(p - 1)
        m = np.zeros((len(caras_p_minus_1), len(caras_p)), dtype=int)
        for d in caras_p:
            for c in caras_p_minus_1:
                if set(c).issubset(set(d)):
                    m[caras_p_minus_1.index(c)][caras_p.index(d)] = 1
        return m

    def matriz_borde(self, p):
        if p < 0 or p > self.d:
            m = []
        else:
            m = self.matriz_borde_aux(p)
        print(f"Matriz de borde de dimensión {p}:\n{m}")
        return m
    
    # Devuelve directamente la matriz borde p en forma normal de Smith
    def normal_Smith_aux(self, p):
        m = self.matriz_borde_aux(p)
        if m is None or m.size == 0:
            return []
        A = m.copy()
        n_rows, n_cols = A.shape
        n = min(n_rows, n_cols)
        for i in range(n):
            if A[i, i] == 0:
                # Buscar un elemento 1 en el resto de la matriz
                for j in range(i, n_rows):
                    for k in range(i, n_cols):
                        if A[j, k] == 1:
                            # Intercambiar filas i y j
                            A[i], A[j] = A[j].copy(), A[i].copy()
                            # Intercambiar columnas i y k 
                            A[:, i], A[:, k] = A[:, k].copy(), A[:, i].copy()
                            break
                    if A[i, i] == 1:
                        break
            # Si no se ha encontrado ningún 1, salimos del bucle
            if A[i, i] != 1:
                break
            # Eliminar los 1s en la fila y columna i
            for j in range(i + 1, n_cols):
                if A[i, j] == 1:
                    A[:, j] = (A[:, i] + A[:, j]) % 2
            for k in range(i + 1, n_rows):
                if A[k, i] == 1:
                    A[k, :] = (A[i, :] + A[k, :]) % 2
        return A

    def normal_Smith(self, p):
        m_smith = self.normal_Smith_aux(p)
        print(f"Matriz de borde en forma normal de Smith de dimensión {p}:\n{m_smith}")
        return m_smith

    # Podemos sacar los numeros de Betti de la matriz en forma normal de Smith
    # La dimension de Z_p será el nº de columnas de la matriz en p menos su rango
    # La dimension de B_p será el rango de la matriz en p + 1
    # El numero de Betti será Z_p - B_p
    def betti_numbers_aux(self, p):
        if p < 0 or p > self.d:
            return 0
        m_p = self.normal_Smith_aux(p)
        if type(m_p) is not type(np.ndarray([])):
            m_p = np.ndarray(m_p)
        m_p_plus_1 = self.normal_Smith_aux(p + 1)
        if type(m_p_plus_1) is not type(np.ndarray([])):
            m_p_plus_1 = np.ndarray(m_p_plus_1)
        dim_p_ciclos = m_p.shape[1] - np.linalg.matrix_rank(m_p)
        if p < self.d:
            dim_p_bordes = np.linalg.matrix_rank(m_p_plus_1)
        else:
            dim_p_bordes = 0
        return dim_p_ciclos - dim_p_bordes

    def betti(self, p):
        beta_p = self.betti_numbers_aux(p)
        print(f"Número de Betti β_{p}: {beta_p}")
        return beta_p
###################################### Fin Clase 6 ######################################

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

    def n_caras(self, n):
        if n < 0 or n > self.dimension:
            return Exception
        else:
            caras_n = sorted(set([cara for cara in self.caras if len(cara) == n+1]), key=lambda x: x)
            return caras_n

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
        self.simplices_ordenados = sorted(self.simplices, key=lambda x: (x.index, x.dimension, x.vertices))

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


    ###################################### Clase 4 ######################################

# Añadimos la clase Punto, la cual contiene los campos vértice y coordenadas.
class Punto:
    def __init__(self, nombre, coords):
        self.vertice = nombre
        self.coords = np.array(coords)

    def distancia(self, otro):
        return np.linalg.norm(self.coords - otro.coords)

    def __repr__(self):
        return f"{self.vertice}{tuple(self.coords)}"

# Clase de Complejo de Vietoris-Rips
class Complejo_Vietoris_Rips:
    def __init__(self, points):
        self.puntos = points  # points es una lista de objetos Punto
        self.space_dimension = len(points[0].coords) if points else 0 # Dimensión del espacio en el que están los puntos



    def r_complex_aux(self, r):
        complex = Complejo_simplicial_filtrado([])
        simplices = []
        for i in range(1, self.space_dimension + 2): # El simplice más grande es con n+1 puntos
            for c in combinations(self.puntos, i): # Todas las combinaciones de i puntos
                if self.verifica_radio(c, r): # dist <= 2r
                    simplices.append(Simplice_filtrado(self.vertices(c), r))
        complex.insert_filtrado(simplices, r) # Los añadimos con tiempo r
        return complex
    
    def r_complex(self, r):
        complex = self.r_complex_aux(r)
        print(f"Complejo de Vietoris-Rips para r = {r}: {[ (s.vertices, s.index) for s in complex.simplices_ordenados ]}")
        return complex
    
    # Comprueba que no haya ninguna distancia entre puntos mayor a 2r
    def verifica_radio(self, puntos, r):
        for p1, p2 in combinations(puntos, 2):
            if p1.distancia(p2) > 2 * r: # dist <= 2r
                return False
        return True
    
    def vertices(self, puntos):
        return [p.vertice for p in puntos]

###################################### Clase 5 ######################################

# Debemos crearnos una función que calcule la filtración de alfa-complejos asociada a un conjunto de puntos en el plano
class AlfaComplejo:
    def __init__(self, points, radius):
        self.puntos = points  # points es una lista de objetos Punto
        self.coords_puntos = np.array([p.coords for p in points])
        self.complex = self.alfa_complejo(radius)
        self.radius = radius

    def r_circuncirculo(self, s, puntos):
        A = puntos[s[0]].coords
        B = puntos[s[1]].coords
        C = puntos[s[2]].coords
        # Calcular el circuncentro y el radio
        M1 = (A + B) / 2
        m1 = (B[1] - A[1])/(B[0] - A[0]) if (B[0] - A[0]) != 0 else None
        M2 = (B + C) / 2
        m2 = (C[1] - B[1])/(C[0] - B[0]) if (C[0] - B[0]) != 0 else None
        if m1 is not None and m2 is not None:
            x_circ = (m1 * M1[0] - m2 * M2[0] + M2[1] - M1[1]) / (m1 - m2)
            y_circ = m1 * (x_circ - M1[0]) + M1[1]
            r = np.linalg.norm([x_circ - A[0], y_circ - A[1]])
            return r
        return None

    def alfa_complejo(self, r):
        Del = Delaunay(self.coords_puntos)  
        simplices = []
        for s in Del.simplices:  # lista de triángulos de Delaunay
            radio_circuncirculo = self.r_circuncirculo(s, self.puntos)
            if radio_circuncirculo is not None and radio_circuncirculo <= r:
                simplices.append(Simplice_filtrado([s[0], s[1], s[2]], r))
            # Comprobamos las aristas
            else:
                arr_dist_aristas = [self.puntos[s[i]].distancia(self.puntos[s[(i+1)%3]]) for i in range(3)] # [6 3 9] -> [63 39 96] -> [dist(6,3), dist(3,9), dist(9,6)]
                for i in range(3):
                    d = arr_dist_aristas[i]
                    if d <= 2*r: # Primer caso de aristas
                        simplices.append(Simplice_filtrado(sorted([s[i], s[(i+1)%3]]), r))
                        # No hay segundo caso ya que la arista se añadira con el triangulo, y puesto que hemos comprobado antes si se añade el triángulo o no
                        # no hace falta volver a comprobarlo.
                    else: # Si no están ni la arista ni el triángulo, añadimos los vértices
                        simplices.append(Simplice_filtrado([s[i]], r))
        # Comprobamos que no haya ninguna cara que ya sea añadida por otra del complejo
        s_aux = []
        # Recorremos los símplices ordenados por dimensión decreciente
        for s in sorted(simplices, key=lambda x: len(x.vertices), reverse=True):
            # Una cara ya está añadida si sus vértices son un subconjunto de los vértices de algún símplice ya añadido
            if any(set(s.vertices).issubset(set(existing.vertices)) for existing in s_aux):
                continue
            # Si no está añadida, la añadimos
            s_aux.append(s)
        complejo = Complejo_simplicial_filtrado([])
        complejo.insert_filtrado(s_aux, r)
        return complejo

    def show_voronoi_delaunay(self):
        vor = Voronoi(self.coords_puntos)
        Del = Delaunay(self.coords_puntos)
        fig = voronoi_plot_2d(vor,show_vertices=False,line_width=2, line_colors='blue' )
        c=np.ones(len(self.coords_puntos))
        cmap = matplotlib.colors.ListedColormap("limegreen")
        plt.tripcolor(self.coords_puntos[:,0],self.coords_puntos[:,1],Del.simplices, c, edgecolor="k", lw=2,
        cmap=cmap)
        plt.plot(self.coords_puntos[:,0], self.coords_puntos[:,1], 'ko')
        plt.show()
    
    # Dibuja el diagrama de Voronoi junto con el alfa-complejo
    def show_voronoi_alfa(self):
        vor = Voronoi(self.coords_puntos)
        fig = voronoi_plot_2d(vor, show_vertices=False, line_width=2, line_colors='blue')
        # Ya dibujado el diagrama de Voronoi, dibujamos el alfa-complejo
        p = self.coords_puntos
        aristas = set()
        triangulos = set()
        for s in self.complex.simplices_ordenados:
            # Intentamos extraer las aristas y triángulos
            try:
                for e in s.n_caras(1):
                    aristas.add(tuple(int(v) for v in e))
                for t in s.n_caras(2):
                    triangulos.add(tuple(int(v) for v in t))
            except Exception:
                pass
        # Dibujamos los triángulos si existen
        if triangulos:
            for tri in triangulos:
                coords = p[list(tri)]
                plt.fill(coords[:, 0], coords[:, 1], facecolor='limegreen', edgecolor='k', alpha=0.3)
        # Dibujamos las aristas, incluidas las que ya dibujo el triangulo para que todas tengan la misma forma
        if aristas:
            for edge in aristas:
                xs = p[list(edge), 0]
                ys = p[list(edge), 1]
                plt.plot(xs, ys, color='k', linewidth=2)
        # Por ultimo dibujamos los puntos
        plt.plot(p[:, 0], p[:, 1], 'ko')
        # Obtenemos el sistema de coordenadas haciendo que una unidad en x valga lo mismo que una en y
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    # Imprime el complejo alfa
    def print_complex(self):
        print(f"Alfa-complejo con radio {self.radius}: {[ (s.vertices, s.index) for s in self.complex.simplices_ordenados ]}")
        return self.complex
    
    def print_points(self):
        print("Puntos del alfa-complejo:")
        for pt in self.puntos:
            x, y = pt.coords
            print(f"{pt.vertice}: ({float(x):.4f}, {float(y):.4f})")
        return self.puntos

###################################### Ejemplo de Uso ######################################
if __name__ == "__main__":
    # # Definimos los símplices maximales del complejo
    # # Por manejo más sencillo, los vértices son enteros
    # s1 = Simplice([0,1,2])
    # s2 = Simplice([2,3])
    # s3 = Simplice([4])
    # # Creamos el complejo simplicial
    # complejo = Complejo_simplicial([s1, s2, s3])
    # # Probamos los métodos
    # print("####################################################")
    # print("#         Ejercicio Complejos Simpliciales         #")
    # print("####################################################")
    # complejo.caras()
    # complejo.caras_por_dimension()
    # complejo.dimension()
    # complejo.Euler()
    # complejo.estrella((2,))
    # complejo.estrella_cerrada((2,))
    # complejo.link((2,))
    # complejo.componentes_conexas()
    # complejo.j_esqueleto(1)
    # complejo.connected_components()
    # complejo.es_conexo()
    # print("####################################################")
    # print("#    Ejercicio Complejos Simpliciales Filtrados    #")
    # print("####################################################")
    # csf = Complejo_simplicial_filtrado([])
    # csf.insert_filtrado([s1, s2], 0)
    # csf.insert_filtrado([s3], 1)
    # csf.insert_filtrado([s1], 0.5)
    # csf.caras()
    # csf.caras_por_dimension()
    # csf.simplices_por_filtrado(0)
    # csf.simplices_por_filtrado(1)
    # print(f"Simplices ordenados: {[(s.vertices, s.index) for s in csf.simplices_ordenados]}")
    # print("####################################################")
    # print("#             Ejercicios Vietoris-Rips             #")
    # print("####################################################")
    # vr = Complejo_Vietoris_Rips([Punto(0, (0,0)), Punto(1, (1,0)), Punto(2, (0,1)), Punto(3, (1,1))])
    # vr.r_complex(0)
    # vr.r_complex(0.25)
    # vr.r_complex(0.5)
    # vr.r_complex(1)
    print("####################################################")
    print("#             Ejercicios Alfa-Complejos            #")
    print("####################################################")
    points = np.random.rand(10,2)
    p = [Punto(i, points[i]) for i in range(len(points))]
    ac = AlfaComplejo(p, 0.25)
    ac.print_complex()
    ac.print_points()
    # ac.show_voronoi_delaunay()
    ac.show_voronoi_alfa()
    print("####################################################")
    print("#              Ejercicios Matriz Borde             #")
    print("####################################################")
    s1 = Simplice([0,1])
    s2 = Simplice([1,2,3,4])
    s3 = Simplice([4,5])
    s4 = Simplice([4,6])
    s5 = Simplice([5,6])
    s6 = Simplice([6,7,8])
    s7 = Simplice([8,9])
    diapositiva_4 = Complejo_simplicial([s1, s2, s3, s4, s5, s6, s7])
    diapositiva_4.matriz_borde(1)
    print("####################################################")
    print("#              Ejercicios Normal Smith             #")
    print("####################################################")
    s1 = Simplice([0,1,2,3])
    tetraedro = Complejo_simplicial([s1])
    tetraedro.normal_Smith(1)
    print("####################################################")
    print("#                Ejercicios Betti                  #")
    print("####################################################")
    print("----------------------------------------------------")
    print("Numeros de Betti del tetraedro")
    print("----------------------------------------------------")
    for i in range(tetraedro.d + 1):
        tetraedro.betti(i)
    print("----------------------------------------------------")
    print("Numeros de Betti del ejemplo de la diapositiva 4")
    print("----------------------------------------------------")
    for i in range(diapositiva_4.d + 1):
        diapositiva_4.betti(i)