import heapq
def clarkson(g, Vertices, Edges, w, d):#Clarkson's Greedy Agortihm

    #initialization
    V = set([])
    for v in Vertices:
        V.add(v)
    E = set([])
    for e in Edges:
        E.add(e)
#    Y = {}
#    for e in E:
#        Y[e] = 0
    C = set([])
    W = {}
    for v in w:
        W[v] = w[v]
    D = {}
    for v in d:
        D[v] = d[v]
    #main loop
    while E != set([]):
        #pick v from V for which W[v]/D[v] is minimized
        WD = {}
        for v in V:
            WD[v] = float(W[v])/D[v]
        v = min(WD.keys(), key=(lambda k:WD[k]))
        u_set = g[v]
        for u in u_set:
            if v < u:
                e = (v, u)
            else:
                e = (u, v)
            if e in E:
                E.remove(e)
            W[u] = W[u] - float(W[v]/D[v])
            D[u] = D[u] - 1
            if D[u]==0 and u in V:
                V.remove(u)
#            Y[e] = float(W[v]/D[v])
        C.add(v)
        V.remove(v)
        W[v] = 0
    return C

def VC(g, w):
    Vertices = set([])
    Edges = set([])
    d = {}
    for v in g:
        Vertices.add(v)
        u_set = g[v]
        for u in u_set:
            if v < u:
                Edges.add((v, u))
            else:
                Edges.add((u, v))
        d[v] = len(u_set)
    print g
    print d
    C = clarkson(g, Vertices, Edges, w, d)
    return C


