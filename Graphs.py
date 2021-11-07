class Vertex:

    def __init__(self, elmt, lat=None, long=None):
        self.elmt = elmt
        self.lat = lat
        self.long = long

    def __str__(self):
        return str(self.elmt)

    def __lt__(self, vertex):
        return self.elmt < vertex.element()

    def element(self):
        return self.elmt


class Edge:

    def __init__(self, v1, v2, element):
        self._v1 = v1
        self._v2 = v2
        self._vertices = (v1, v2)
        self.elmt = element

    def __str__(self):
        return ('(' + str(self._vertices[0]) + '--'
                + str(self._vertices[1]) + ' : '
                + str(self.elmt) + ')')

    def vertices(self):
        return self._vertices

    def start(self):
        """ Return the first vertex in the ordered pair. """
        return self._vertices[0]

    def end(self):
        """ Return the second vertex in the ordered pair. """
        return self._vertices[1]

    def opposite(self, v1):
        if self._vertices[0] == v1:
            return self._vertices[1]
        elif self._vertices[1] == v1:
            return self._vertices[0]
        else:
            return None

    def element(self):
        return self.elmt


class Graph:

    def __init__(self):
        self._structure = {}

    def vertices(self):
        return list(self._structure)

    def edges(self):
        """ Return a list of all edges in the graph. """
        edgelist = []
        for v in self._structure:
            for w in self._structure[v]:
                # to avoid duplicates, only return if v is the first vertex
                if self._structure[v][w].start() == v:
                    edgelist.append(self._structure[v][w])
        return edgelist

    def num_vertices(self):
        return len(self._structure)

    def num_edges(self):
        num_edges = 0
        for v1 in self._structure:
            num_edges += len(self._structure[v1])
        return num_edges // 2

    def get_vertex_by_label(self, element):
        """ Return the first vertex that matches element. """
        for v in self._structure:
            if v.element() == element:
                return v
        return None

    def get_edge(self, v1, v2):
        if (self._structure is not None
                and v1 in self._structure
                and v2 in self._structure[v1]):
            return self._structure[v1][v2]

    def degree(self, v1):
        return len(self._structure[v1])

    def get_edges(self, v1):
        if v1 in self._structure:
            edges = []
            for v2 in self._structure[v1]:
                edges.append(self._structure[v1][v2])
            return edges
        return None

    def add_vertex(self, elmt, lat, long):
        new_vertex = Vertex(elmt)
        self._structure[new_vertex] = {}
        return new_vertex

    def add_vertex_if_new(self, element):
        for v in self._structure:
            if v.element() == element:
                return v
        return self.add_vertex(element, None, None)

    def add_edge(self, v1, v2, elmt):
        if v1 not in self._structure or v2 not in self._structure:
            return None
        new_edge = Edge(v1, v2, elmt)
        self._structure[v1][v2] = new_edge
        self._structure[v2][v1] = new_edge
        return new_edge

    def add_edge_pairs(self, elist):
        """ add all vertex pairs in elist as edges with empty elements.

        Args:
            elist - a list of pairs of vertex objects
        """
        for (v, w) in elist:
            self.add_edge(v, w, None)

    def remove_vertex(self, v1):
        if v1 in self._structure:
            self._structure.pop(v1)

    def remove_edge(self, name):
        v1 = name.start()
        v2 = name.end()
        self._structure[v1].pop(v2)
        self._structure[v2].pop(v1)

    def highestdegreevertex(self):
        """ Return the vertex with highest degree. """
        hd = -1
        hdv = None
        for v in self._structure:
            if self.degree(v) > hd:
                hd = self.degree(v)
                hdv = v
        return hdv

    def __str__(self):
        """ Return a string representation of the graph. """
        hstr = ('|V| = ' + str(self.num_vertices())
                + '; |E| = ' + str(self.num_edges()))
        vstr = '\nVertices: '
        for v in self._structure:
            vstr += str(v) + '-'
        vstr = vstr.strip('-')
        edges = self.edges()
        estr = '\nEdges: '
        for e in edges:
            estr += str(e) + ' '
        return hstr + vstr + estr

    def depthFirstSearch(self, vertex):
        results = {str(vertex): "Start"}
        self._depthFirstSearch(vertex, results)
        return results

    def _depthFirstSearch(self, vertex, marked):
        for v2 in self._structure[vertex]:
            if str(v2) not in marked:
                marked[str(v2)] = str(self.get_edge(vertex, v2))
                self._depthFirstSearch(v2, marked)

    def breadthFirstSearch(self, vertex):
        results = []
        queue = [vertex]
        while queue:
            vertex = queue.pop(0)
            if str(vertex) not in results:
                results.append(str(vertex))
                neighbours = self._structure[vertex]

                for neighbour in neighbours:
                    queue.append(neighbour)
        return results

    def Dijkstra(self, s):
        open = APQ()
        locs = dict()
        closed = dict()
        preds = {s: None}
        start_point = open.add(0, s)
        locs[s] = start_point
        while open.length() > 0:
            v = open.remove_min()
            locs.pop(v.value)
            predecessor = preds.pop(v.value, None)
            closed[v.value] = (v.cost, predecessor)
            for edge in self.get_edges(v.value):
                w = edge.opposite(v.value)
                if w not in closed.keys():
                    newcost = v.cost + edge.elmt
                    if w not in locs.keys():
                        preds[w] = v.value
                        next_vertex = open.add(newcost, w)
                        locs[w] = next_vertex
                    elif newcost < locs[w].cost:
                        preds[w] = v.value
                        open.update_key(locs[w], newcost)
        return closed


class Element:

    def __init__(self, k, v, i):
        self.cost = k
        self.value = v
        self.index = i

    def __eq__(self, other):
        return self.cost == other.cost

    def __lt__(self, other):
        return self.cost < other.cost

    def __str__(self):
        return str(self.value)


class APQ:

    def __init__(self):
        self.queue = []

    def bubbleup(self, i):
        while i > 0:
            parent = (i - 1) // 2
            if self.queue[i] < self.queue[parent]:
                self.queue[i], self.queue[parent] = self.queue[parent], self.queue[i]
                self.queue[i].index = i
                self.queue[parent].index = parent
                i = parent
            else:
                i = 0

    def bubbledown(self, i, last):
        while last > (i * 2):
            lc = i * 2 + 1
            rc = i * 2 + 2
            minc = lc
            if last > lc and self.queue[rc] < self.queue[lc]:
                minc = rc
            if self.queue[i] > self.queue[minc]:
                self.queue[i], self.queue[minc] = self.queue[minc], self.queue[i]
                self.queue[i].index = i
                self.queue[minc].index = minc
                i = minc
            else:
                i = last

    def add(self, cost, item):
        item_added = Element(cost, item, len(self.queue))
        self.queue.append(item_added)
        if len(self.queue) > 1:
            self.bubbleup(item_added.index)
        return item_added

    def min(self):
        return self.queue[0]

    def remove_min(self):
        ret = self.queue[0]
        if len(self.queue) == 1:
            self.queue.pop()
        elif len(self.queue) > 1:
            self.queue[0] = self.queue.pop()
            self.queue[0].index = 0
            self.bubbledown(0, len(self.queue) - 1)
        return ret

    def length(self):
        return len(self.queue)

    def update_key(self, element, newcost):
        index = element.index
        self.queue[index].cost = newcost
        # compare keys & bubbleup
        parent = (index - 1) // 2
        if self.queue[index].cost < self.queue[parent].cost:
            self.bubbleup(index)

    def get_key(self, element):
        return element.cost


class RouteMap(Graph):

    def __init__(self):
        super().__init__()
        self.keyElements = {}

    def __str__(self):
        if self.num_edges() + self.num_vertices() > 100:
            return "Graph too large to represent"
        else:
            return super().__str__()

    def sp(self, v, w):
        output = []
        vertices = self.Dijkstra(v)
        current = w
        while current != v:
            output.append(current)
            current = vertices[current][1]
        output.append(v)
        output.reverse()
        return output

    def get_vertex_by_label(self, element):
        vertex = self.keyElements[str(element)]
        return vertex

    def add_vertex(self, elmt, lat, long):
        new_vertex = Vertex(elmt, lat, long)
        self._structure[new_vertex] = {}
        self.keyElements[str(elmt)] = new_vertex
        return new_vertex

    def printSP(self, v, w):
        shortPath = self.sp(v, w)
        allShortest = self.Dijkstra(v)
        print("type" + "\t" + "latitude" + "\t" + "longitude" + "\t" + "element" + "\t" + "cost")
        for label in shortPath:
            vertex = self.get_vertex_by_label(label)
            print("W \t %s \t %s \t %s \t %s" % (vertex.lat, vertex.long, vertex, allShortest[vertex][0]))


def graphreader(filename):
    """ Read and return the route map in filename. """
    graph = RouteMap()
    file = open(filename, 'r')
    entry = file.readline()  # either 'Node' or 'Edge'
    num = 0
    while entry == 'Node\n':
        num += 1
        nodeid = int(file.readline().split()[1])
        coords = file.readline().split()[1:]
        lat, long = float(coords[0]), float(coords[1])
        vertex = graph.add_vertex(nodeid, lat, long)
        entry = file.readline()  # either 'Node' or 'Edge'
    print('Read', num, 'vertices and added into the graph')
    num = 0
    while entry == 'Edge\n':
        num += 1
        source = int(file.readline().split()[1])
        sv = graph.get_vertex_by_label(source)
        target = int(file.readline().split()[1])
        tv = graph.get_vertex_by_label(target)
        length = float(file.readline().split()[1])
        time = float(file.readline().split()[1])
        edge = graph.add_edge(sv, tv, time)
        file.readline()  # read the one-way data
        entry = file.readline()  # either 'Node' or 'Edge'
    print('Read', num, 'edges and added into the graph')
    print(graph)
    return graph


if __name__ == "__main__":
    def test_graph():
        graph = graphreader("corkCityData.txt")
        vertex1 = graph.get_vertex_by_label(3777201945)
        vertex2 = graph.get_vertex_by_label(330068634)
        graph.printSP(vertex1, vertex2)



    test_graph()
