# load("GitHub/vertex-minors/LMQC.sage")
load("LMQC.sage")
#Finding all graphs T LMQC to Bell 2
V = 4
Bell2 = SimpleGraphLMQC({0:[3],1:[2]})
T = [[0,1],[2],[3]]
Bell2.set_partition(T)

T_LMQC_to_Bell2 = []

for i in powerset(combinations(range(V),2)):
    G = Graph(V)
    G.add_edges(i)
    if Bell2.is_LMQC_eq(G)[0]:
        T_LMQC_to_Bell2.append(G)
print "There are ", len(T_LMQC_to_Bell2), "T-LMQC equivalent to Bell2"

#Lemma 4.7
lemma_true = True
for G in T_LMQC_to_Bell2:
    for i in [2,3]:
        if not (G.has_edge(i,0) or G.has_edge(i,1)):
            lemma_true = False

print "Lemma 4.7 is ", lemma_true
