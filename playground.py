from pyvis.network import Network

#from graph import *
from transitions import *
'''
G = nx.complete_bipartite_graph(4, 2)

#################

G.add_node(6)
G.add_node(7)
G.add_node(8)

G.add_edge(4,5)
G.add_edge(6,4)
G.add_edge(6,5)
G.add_edge(7,4)
G.add_edge(7,5)
G.add_edge(6,8)
G.add_edge(7,8)

#G = nx.complete_graph(3)
#R = 2
#G = assign_random_labels(G, [0, 1])
PVG2 = Network()
PVG2.from_nx(G)
PVG2.show(f'prima_es_un_muerto.html', notebook = False)

GStates = apply_transitions(G, mukwembi_transition, parameters = {"R":5}, gens = 10, fullList = True)

PVG = Network()
PVG.from_nx(G)
prepare_drawing(PVG)
PVG.show(f'test0.html', notebook = False)
for i, CG in enumerate(GStates):
    PVG = Network()
    PVG.from_nx(CG)
    prepare_drawing(PVG)
    PVG.show(f'test{i+1}.html', notebook = False)


for R in range(2, 8):
    TG = TransitionGraph(mukwembi_transition,  parameters = {"R":R} )
    TG.expand_all_states(G, valid_states=[0, 1])
    TGX = TG.to_networkx()
    connectedComponents = list(nx.weakly_connected_components(TGX))
    initialStatesCompsCounter = 0
    for cc in connectedComponents:
        hasInitialStates = False
        for n in cc:
            if not '2' in n:
                hasInitialStates = True
                break
        if hasInitialStates: initialStatesCompsCounter += 1
    print(f'R = {R}, connected_components_with_initial_states = {initialStatesCompsCounter}, connectedComponents = {len(connectedComponents)}')

    #build_interactive_view(G, TG, TGX, filename=f"combined{R}.html", node_color_function=mukwembi_transition_coloring, graph_name=f"K{3}, R={R}")
'''

G = nx.complete_graph(3)
R = 1
G = assign_random_labels(G, [0, 1])
TG = TransitionGraph(mukwembi_transition,  parameters = {"R":R} )
TG.expand_all_states(G, valid_states=[0, 1, 2])
TGX = TG.to_networkx()
build_interactive_view(G, TG, TGX, filename=f'K{3}R{R}.html', node_color_function=mukwembi_transition_coloring, graph_name=f"K{3}, R={R}")