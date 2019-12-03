__author__ = '@gavruskin, Lena Collienne, Kieran Elmes'

from collections import deque
import time
import numpy as np

# Graph methods.

def Seidel_recursive(A, count):
    print("Going down: {} at {}".format(count, time.strftime("%D %H:%M:%S")))
    n = np.shape(A)[0]
    if all(A[i][j] for i in range(n) for j in range(n) if i != j):
        return A
    print("starting dot")
    Z = np.dot(A, A)
    print("done dot, starting B")
    B = np.array([
        [1 if i != j and (A[i][j] == 1 or Z[i][j] > 0) else 0 for j in range(n)]
    for i in range(n)])
    print("done B")
    T = Seidel_recursive(B, count + 1)
    X = np.dot(T,A)
    degree = [sum(A[i][j] for j in range(n)) for i in range(n)]
    D = np.array([
        [2 * T[i][j] if X[i][j] >= T[i][j] * degree[j] else 2 * T[i][j] - 1 for j in range(n)]
    for i in range(n)])
    print("Coming up: {} at {}".format(count, time.strftime("%D %H:%M:%S")))
    return D

class Graph(object):
    def __init__(self, vertices, edges):
        # graph is a set of vertices and a list of edges, which a unordered pairs of vertices.
        self.vertices = vertices
        self.edges = edges

    def print_graph(self):
        print("Vertices:")
        for i in self.vertices:
            print(i)
        print("\n")
        print("Edges:")
        for i in self.edges:
            print(i)
        print("\n")

    def print_graph_to_file(self,filename):
        # Prints graph to file called filename (creates one if it doesn't exist).
        # Unlike print_graph, prints edges as ordered pairs.
        f = open(filename, 'w')
        for i in self.edges:
            f.write(str(i) + "\n")
        f.close()

    def adjacency_list(self):
        # returns the adjacency list of the given graph
        adj = dict()
        for e in self.edges:
            # Edges are strings '(x,y)'
            e = e.replace('(','')
            e = e.replace(')','')
            e = e.split(',')
            if int(e[0]) in adj.keys():
                adj[int(e[0])].add(e[1])
            else:
                adj[int(e[0])] = set()
                adj[int(e[0])].add(e[1])

            if int(e[1]) in adj.keys():
                adj[int(e[1])].add(e[0])
            else:
                adj[int(e[1])] = set()
                adj[int(e[1])].add(e[0])
        return adj

    def shortest_path(self, start, end):
            # Returns the distance between two vertices in the graph using breadth first search and a shortest path between the given nodes
            distance = dict()
            pred = dict()
            visited = dict()
            queue = deque([])
            # Initialization of the dicts distance, pred and visited
            for i in self.vertices:
                distance[i] = float("inf")
                pred[i] = -1
                visited[i] = False
            distance[start] = 0
            queue.append(start)
            visited[start] = True
            k = 0 #variable to print the distance of the currently visiting node (to gain information on progress of algorithm)
            currentdist = list()
            while len(queue) > 0:
                node = queue.popleft()
                # Stop when end node is reached, return distance between start and end node
                if node == end:
                    # Return the pair (distance, shortest path):
                    current_node = end
                    path = [current_node]
                    while current_node is not start:
                        current_node = pred[current_node]
                        path.append(current_node)
                    path.reverse()
                    return (distance[end], path)
                adj = self.adjacency_list()
                for next in adj[node]:
                    next = int(next)
                    if visited[next] == False:
                        pred[next] = node
                        visited[next] = True
                        distance[next] = distance[node] + 1
                        queue.append(next)

    def Floyd_Warshall(self):
        # PROBLEM: indices here are not the indices of the trees in the uRNNI graph !!!!!
        #FLOYD-WARSHALL algorithm: returns all distances between all pairs of taxa as a nxn distance matrix
        n = len(self.vertices)
        m = len(self.edges)
        index = dict() #Save the vertices in a dictionary where they are assigned their positions in the distance matrix (needed if vertices don't have names 1,..,n)
        counter = 0
        #Initalise the distance matrix with values inf
        distance = np.array(np.ones((n,n)) * np.inf)
        for i in self.edges:
            edge = i.replace('(','').replace(')','').split(',') #edges are saved as strings > convert them to list of integers
            edge[0] = int(edge[0])
            edge[1] = int(edge[1])

            # a = index of edge[0]
            if edge[0] in index:
                a = index[edge[0]]
            else:
                index[edge[0]] = counter
                a = counter
                counter = counter + 1
            # b = index of edge[0]
            if edge[1] in index:
                b = index[edge[1]]
            else:
                index[edge[1]] = counter
                b = counter
                counter = counter + 1
            distance[a][b] = 1 #Can/Should we reduce the distance matrix to its half?
            distance[b][a] = 1
        for i in self.vertices:
            distance[index[i]][index[i]] = 0
        for k in range(0,n):
            for i in range(0,n):
                for j in range(0,n):
                    if distance[i][j] > distance[i][k] + distance[k][j]:
                        distance[i][j] = distance[i][k] + distance[k][j]
        # Save the distances in output/distance_matrix_n_taxa.txt, indices are given in dict that is returned when computing RNNI graoh
        d_file = open('output/distance_matrix.txt', 'w')
        d_file.write(str(distance))
        d_file.close()
        return distance, index


    def Seidel(self):
        # Convert the graph into a numpy adjacency matrix
        n = len(self.vertices)
        m = len(self.edges)
        index = dict() #Save the vertices in a dictionary where they are assigned their positions in the distance matrix (needed if vertices don't have names 1,..,n)
        counter = 0

        adjacency = np.zeros((n,n), dtype=bool)
        for i in self.edges:
            edge = i.replace('(','').replace(')','').split(',') #edges are saved as strings > convert them to list of integers
            edge[0] = int(edge[0])
            edge[1] = int(edge[1])

            # a = index of edge[0]
            if edge[0] in index:
                a = index[edge[0]]
            else:
                index[edge[0]] = counter
                a = counter
                counter = counter + 1
            # b = index of edge[0]
            if edge[1] in index:
                b = index[edge[1]]
            else:
                index[edge[1]] = counter
                b = counter
                counter = counter + 1
            adjacency[a][b] = 1 #Can/Should we reduce the distance matrix to its half?
            adjacency[b][a] = 1
        for i in self.vertices:
            adjacency[index[i]][index[i]] = 0

        # Find all-pairs distance
        A = Seidel_recursive(adjacency, 0)
        return A, index


    def Dijkstra(self, start):
        #Dijkstra algortihm for getting shortest path from one node (start) to all others
        Q = dict()
        dist = dict()
        prev = dict()
        adj = self.adjacency_list()
        for v in self.vertices:
            dist[v] = float('inf')
            prev[v] = None
            Q[v] = dist[v]
        dist[start] = 0
        Q[start] = 0
        while len(Q)!=0:
            u = min(Q, key=lambda x: Q.get(x))
            del Q[u]
            for neighbour in adj[int(u)]:
                neighbour = int(neighbour)
                alt = dist[u] + 1
                if alt < dist[int(neighbour)]:
                    dist[neighbour] = alt
                    Q[neighbour] = alt
                    prev[neighbour] = u
        return dist, prev

    def subgraph(self, S):
        #Compute a subgraph defined by a set of vertices
        H = Graph(S,set())
        for e in self.edges:
            f = e
            e = e.replace('(','')
            e = e.replace(')','')
            e = e.split(',')
            if (int(e[0]) in S) and (int(e[1]) in S):
                H.edges.add(f)
        return H
