import copy
import networkx as nx

def mukwembi_transition(_G, R = 2):
    """
    Apply one step of Mukwembi's epidemic spread transition to a graph.

    Computes the next state of every node based on its current 'group'
    and the groups of its neighbors in `_G`, using a healthy -> infected
    -> dead -> (infected or healthy) cycle:

        0: healthy  — becomes infected (1) if any neighbor is infected.
        1: infected — always becomes dead (2) on the next step.
        2: dead     — becomes infected (1) again if at least R neighbors
                      are infected; otherwise reverts to healthy (0).

    All next-state values are computed from the *original* input graph
    `_G` (not from intermediate updates within the same call), so the
    transition is applied simultaneously across all nodes — order of
    iteration does not affect the result.

    Parameters
    ----------
    _G : nx.Graph
        The current graph state. Each node must have a 'group' attribute
        with value 0 (healthy), 1 (infected), or 2 (dead). Not mutated —
        a deep copy is made internally and returned.
    R : int, default=2
        Resistance threshold: the minimum number of infected neighbors
        required for a dead node to become infected again instead of
        reverting to healthy.

    Returns
    -------
    nx.Graph
        A new graph (deep copy of `_G`) with each node's 'group' updated
        according to the transition rules above.
    """
    G = copy.deepcopy(_G)
    for i in range(len(G.nodes)):
        state = _G.nodes[i]["group"]
        neighborsStates = [_G.nodes[j]["group"] for j in _G.neighbors(i)]
        if state == 0:
            if 1 in neighborsStates: G.nodes[i]["group"] = 1
        elif state == 1: G.nodes[i]["group"] = 2
        else:
            if neighborsStates.count(1) >= R: G.nodes[i]["group"] = 1
            else: G.nodes[i]["group"] = 0
    return G

def mukwembi_transition_coloring(label):
    """
    Choose a display color for a transition-graph state node based on its
    label (a tuple/sequence of per-node 'group' values from the original
    base graph at that state).

    Rules
    -----
    - If any node in the state is dead (group 2), color green (#62D484).
    - Else if any node is infected (group 1), color purple/blue (#9797F1).
    - Else (every node is healthy), color white (#FFFFFF).

    Intended for use as the `node_color_function` passed to
    `build_interactive_view`, to visually distinguish transition-graph
    states by severity at a glance.

    Parameters
    ----------
    label : tuple or sequence
        The per-node group values making up a single state (e.g. as
        stored in `TGX.nodes[n]["state"]`).

    Returns
    -------
    str
        A hex color code representing the state's severity.
    """
    if 2 in label:
        return "#62D484"
    elif 1 in label:
        return "#9797F1"
    else:
        return "#FFFFFF"

def is_not_vanishing(TG):
    """
    Check whether a digraph has more than one weakly connected component
    containing initial states (a node whose id does not include '2').

    Parameters
    ----------
    TGX : nx.DiGraph
        A transition digraph, complete or partially expanded.

    Returns
    -------
    bool
        True if more than one weakly connected component contains an initial state ('2' in its node id). False otherwise.
    """
    TGX = TG.to_networkx()
    components = nx.weakly_connected_components(TGX)
    initialStatesCounter = 0
    for cc in components:
        for n in cc:
            if not '2' in n: #It is an initial state
                initialStatesCounter += 1
                break
        if initialStatesCounter > 1: break
    return initialStatesCounter > 1