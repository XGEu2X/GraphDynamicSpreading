import copy
import networkx as nx
import random
import json
import itertools

from mukwembi_utils import *

def assign_random_labels(G, valid_states):
    """
    Randomly assign a label from valid_states to each node in G.

    Each node in G is independently assigned a value drawn from
    `valid_states` via random.choice, stored under both a 'label' and a
    'group' attribute (identical values — 'group' is used for numeric/
    logic purposes such as transition rules, 'label' for display). Node
    identities and graph structure are unchanged; only attributes are set.

    Parameters
    ----------
    G : nx.Graph
        The graph whose nodes will be labeled. Modified in place.
    valid_states : list
        Pool of possible values to randomly assign to each node.

    Returns
    -------
    nx.Graph
        The same graph G, with 'label' and 'group' node attributes set.
    """
    labels = {node: random.choice(valid_states) for node in G.nodes()}
    nx.set_node_attributes(G, labels, name="label")
    nx.set_node_attributes(G, labels, name="group")
    return G

#Prepares PVG settings for drawing.
def prepare_drawing(PVG, Pallete = ['#ffffff', '#ff0000', '#00ff00']):
    """
    Prepares a pyvis Network's node settings for drawing the base graph.

    Colors each node's background according to its 'group' attribute
    (used as an index into `Pallete`), sets a 1-indexed numeric label
    (position + 1) for display, and applies a large font size for
    readability.

    Assumes PVG's nodes are ordered such that enumeration index `i`
    matches each node's actual id (e.g. sequential integer node ids
    starting at 0), since `PVG.nodes[i]` is looked up directly by that
    index.

    Parameters
    ----------
    PVG : pyvis.network.Network
        The pyvis network (already built from the base graph, e.g. via
        `PVG.from_nx(G)`) whose node display settings will be set.
    Pallete : list of str, default=['#ffffff', '#ff0000', '#00ff00']
        Hex colors indexed by each node's 'group' value (0, 1, 2, ...).

    Returns
    -------
    None
        Modifies PVG's nodes in place.
    """
    for i, node in enumerate(PVG.nodes): 
        PVG.nodes[i]['color'] = { "background": Pallete[node["group"]], "border": "#000000"}
        PVG.nodes[i]['label'] = f'{i + 1}'
        PVG.nodes[i]['font'] = {"size": 50}
    return

def prepare_transition_drawing(PVG, fontsize=20):
    """
    Prepares a pyvis Network's node settings for drawing a transition
    digraph built by TransitionGraph.to_networkx().

    Labels each node with its state tuple (joined into a single string,
    e.g. (0, 1, 2) -> "012"), and colors it by severity based on that
    label's characters: red if any node in the state is dead ('2'
    appears), blue if any node is infected ('1' appears, and no '2'),
    otherwise white (all healthy).

    Parameters
    ----------
    PVG : pyvis.network.Network
        The pyvis network (already built from a transition digraph, e.g.
        via `PVG.from_nx(TGX)`) whose node display settings will be set.
    fontsize : int, default=20
        Font size (in px) used for node labels.

    Returns
    -------
    None
        Modifies PVG's nodes in place.
    """
    for node in PVG.nodes:
        state = node.get("state", node["id"])
        node["label"] = "".join(str(s) for s in state)
        node["font"] = {"size": fontsize}
        if '2' in node["label"]: node["color"] = { "background": '#ff0000', "border": "#000000"}
        elif '1' in node["label"]: node["color"] = { "background": '#0000ff', "border": "#000000"}
        else: node["color"] = { "background": '#ffffff', "border": "#000000"}
    return


def apply_transitions(G, transition_function, parameters, gens = 10, fullList = False):
    """
        Repeatedly apply a transition function to a graph over multiple generations.

        Starting from G, this runs `transition_function` for `gens` iterations, 
        feeding the output of each iteration back in as the input graph for the 
        next (via parameters['G']).

        Parameters
        ----------
        G : nx.Graph
            The initial graph to evolve. It is deep-copied before any transitions
            are applied, so the original G is never mutated.
        transition_function : callable
            A function that takes the current graph (and any other parameters)
            and returns the next graph state. Must accept 'G' as a keyword argument.
        parameters : dict
            Keyword arguments to pass to `transition_function` on each call.
            Note: this dict is mutated in place — 'G' is overwritten with the
            current graph state before each call.
        gens : int, default=10
            Number of generations (iterations) to run the transition function.
        fullList : bool, default=False
            If True, return a list containing the graph state after every
            generation (length == gens). If False, return only the final
            graph state as a single-element list.
    """
    Result = []
    CurrentG = copy.deepcopy(G)
    for i in range(gens):
        FullParameters = parameters | {"G": CurrentG}
        CurrentG = transition_function(**FullParameters)
        if fullList: Result.append(copy.deepcopy(CurrentG))
    if fullList: return Result
    return CurrentG

class TransitionGraph:
    """
    Represents the digraph of transitions between distinct graph states.

    Each node in the internal digraph corresponds to a unique *state*
    (a canonical, hashable representation of a graph's node groups —
    by default the tuple of `group` values ordered by node id). Applying
    `transition_function` to a state produces the next state; an edge is
    added from the current state to the next one. States are deduplicated:
    if a transition leads back to a previously seen state, no new node is
    created and expansion along that path stops (cycle/fixed point reached).

    This lets you build the transition digraph lazily and partially —
    starting from one or more initial graphs — instead of enumerating
    every possible state up front.

    Parameters
    ----------
    transition_function : callable
        Function that maps a graph (via a 'G' keyword argument) to its
        next state, e.g. `mukwembi_transition`.
    parameters : dict, optional
        Extra keyword arguments to pass to `transition_function` on every
        call (e.g. {'R': 2}). 'G' is injected automatically and should not
        be included here.
    state_key : callable, optional
        Function mapping a graph -> hashable state key. Defaults to the
        tuple of each node's 'group' attribute, ordered by node id.
    """

    def __init__(self, transition_function, parameters=None, state_key=None):
        self.transition_function = transition_function
        self.parameters = parameters or {}
        self.state_key = state_key or self._default_state_key
        self.digraph = nx.DiGraph()   # nodes = state keys, edges = transitions
        self.graphs = {}              # state key -> representative graph

    @staticmethod
    def _default_state_key(G):
        return tuple(G.nodes[n]["group"] for n in sorted(G.nodes))

    def add_state(self, G):
        """
        Register a graph's state if not already present.

        Returns the state key, creating a new digraph node (and storing a
        deep copy of G as its representative graph) only if this exact
        state hasn't been seen before.
        """
        key = self.state_key(G)
        if key not in self.graphs:
            self.graphs[key] = copy.deepcopy(G)
            self.digraph.add_node(key)
        return key

    def expand_from(self, G, max_steps=None):
        """
        Grow the transition digraph starting from graph G.

        Repeatedly applies `transition_function` to walk forward from G's
        state, adding a state node and a transition edge at each step.
        Stops when a previously visited state is reached (the resulting
        edge closes the cycle) or when `max_steps` is hit — whichever
        comes first. Already-known states are never recomputed or
        duplicated, so calling this multiple times from different
        starting graphs safely shares/merges the same digraph.

        Parameters
        ----------
        G : nx.Graph
            Starting graph for this expansion.
        max_steps : int, optional
            Maximum number of transitions to apply. None means expand
            until a state repeats.

        Returns
        -------
        list
            The sequence of state keys visited, in order, starting with
            G's own state.
        """
        path = [self.add_state(G)]
        visited = {path[0]}
        steps = 0

        while max_steps is None or steps < max_steps:
            current_key = path[-1]
            current_G = self.graphs[current_key]
            call_params = self.parameters | {"_G": copy.deepcopy(current_G)}
            next_G = self.transition_function(**call_params)
            next_key = self.add_state(next_G)

            self.digraph.add_edge(current_key, next_key)
            path.append(next_key)
            steps += 1

            if next_key in visited:
                break  # cycle (or fixed point) reached — stop expanding
            visited.add(next_key)

        return path

    def get_graph(self, key):
        """Return the representative graph stored for a given state key."""
        return self.graphs[key]

    def states(self):
        """Return all distinct state keys currently known."""
        return list(self.digraph.nodes)

    def transitions(self):
        """Return all transitions (edges) currently known, as (from, to) state key pairs."""
        return list(self.digraph.edges)

    def to_networkx(self, include_graphs=False):
        """
        Build a standalone networkx DiGraph from the recorded states and transitions.

        Node IDs are stringified versions of each state key (pyvis/vis.js
        require string or int node IDs), with the original state tuple kept
        as the 'state' attribute. Edges correspond to observed transitions.

        Parameters
        ----------
        include_graphs : bool, default=False
            If True, attach each state's representative graph as a node
            attribute 'graph' (a deep copy).

        Returns
        -------
        nx.DiGraph
            A digraph with one (string-ID) node per distinct state and one
            edge per observed transition between states.
        """
        def to_id(key):
            return "".join(str(s) for s in key)

        D = nx.DiGraph()
        for key in self.states():
            D.add_node(to_id(key), state=key)
        for src, dst in self.transitions():
            D.add_edge(to_id(src), to_id(dst))

        if include_graphs:
            for key in self.states():
                D.nodes[to_id(key)]["graph"] = copy.deepcopy(self.graphs[key])

        return D

    def expand_all_states(self, base_graph, valid_states, stop_method = None):
        """
        Fully populate the transition digraph for a given base graph by
        expanding from every possible state (i.e. every combination of
        valid_states assigned to the base graph's nodes).

        This guarantees the resulting digraph is complete: every state
        reachable from any possible starting configuration is included,
        not just the ones reachable from a single initial graph.

        Warning: this is exponential in the number of nodes — there are
        len(valid_states) ** n_nodes possible starting states. Only use
        this for small graphs / small state spaces.

        Parameters
        ----------
        base_graph : nx.Graph
            The graph structure (edges) to use for every state. Node
            'group'/'label' attributes on this graph are ignored — every
            possible combination is generated and tried instead.
        valid_states : list
            The set of possible values for each node's 'group' attribute.

        Returns
        -------
        set
            The set of all state keys discovered (equivalent to
            set(self.states()) after this call).
        """
        nodes = sorted(base_graph.nodes)
        n = len(nodes)

        for combo in itertools.product(valid_states, repeat=n):
            key = tuple(combo)
            if key in self.graphs:
                continue  # already covered by a previous expansion

            candidate = copy.deepcopy(base_graph)
            nx.set_node_attributes(
                candidate, dict(zip(nodes, combo)), name="group"
            )
            self.expand_from(candidate)
            if stop_method is not None:
                if stop_method(self): return None

        return set(self.states())

def build_interactive_view(G, TG, TGX, graph_name="Base Graph",
                            pallete=('#ffffff', "#bd5454", "#8cd48c"),
                            node_color_function=None,
                            filename="interactive.html", height="600px"):
    """
    Builds a single HTML page with two linked vis.js networks:
    - Left: the base graph G (fixed structure, colors change on click)
    - Right: the transition digraph TGX (click a state node to preview it)

    Clicking a state node in the transition digraph updates the base
    graph's node colors to match that state.

    Parameters
    ----------
    G : nx.Graph
        The base graph (structure is fixed; only colors change).
    TG : TransitionGraph
        The TransitionGraph instance, used to look up each state's group
        values (TG.graphs[state_key]).
    TGX : nx.DiGraph
        The transition digraph produced by TG.to_networkx() (nodes carry
        a 'state' attribute with the original tuple key).
    graph_name : str
        Name of the base graph, used in the page title:
        "Spread Transition Graph of {graph_name}".
    pallete : tuple of str
        Hex colors indexed by group value (0, 1, 2, ...), used for the
        base graph's nodes.
    node_color_function : callable, optional
        Function that takes a state's label (the tuple stored in
        TGX.nodes[n]["state"]) and returns a hex color string, used to
        color the transition digraph's nodes. Defaults to a flat
        light-gray color for every node if not provided.
    filename : str
        Output HTML file path.
    """
    title = f"Spread Transition Graph of {graph_name}"

    if node_color_function is None:
        node_color_function = lambda label: "#f0f0f0"

    # --- Base graph data (fixed node ids/edges, initial colors from current G) ---
    base_nodes = [
        {"id": n, "label": str(n), "color": pallete[G.nodes[n]["group"]]}
        for n in G.nodes
    ]
    base_edges = [{"from": u, "to": v} for u, v in G.edges]

    # --- Transition digraph data ---
    trans_nodes = [
        {
            "id": n,
            "label": "".join(str(s) for s in TGX.nodes[n]["state"]),
            "color": node_color_function(TGX.nodes[n]["state"]),
        }
        for n in TGX.nodes
    ]
    trans_edges = [{"from": u, "to": v} for u, v in TGX.edges]

    # --- Precompute, for each state node id, the color list for the base graph ---
    base_node_order = list(G.nodes)  # must match the order used in state_key (sorted)
    state_colors = {}
    for n in TGX.nodes:
        key = TGX.nodes[n]["state"]
        state_G = TG.graphs[key]
        state_colors[n] = [
            {"id": node_id, "color": pallete[state_G.nodes[node_id]["group"]]}
            for node_id in base_node_order
        ]

    html_template = f"""
<!DOCTYPE html>
<html>
<head>
  <title>{title}</title>
  <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
  <style>
    body {{ font-family: sans-serif; }}
    h1 {{ text-align: center; margin: 10px 0; }}
    #container {{ display: flex; }}
    #base, #trans {{ width: 50%; height: {height}; border: 1px solid #ccc; }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <div id="container">
    <div id="base"></div>
    <div id="trans"></div>
  </div>

  <script>
    const baseNodes = new vis.DataSet({json.dumps(base_nodes)});
    const baseEdges = new vis.DataSet({json.dumps(base_edges)});
    const baseNetwork = new vis.Network(
        document.getElementById("base"),
        {{ nodes: baseNodes, edges: baseEdges }},
        {{ physics: {{ enabled: true }} }}
    );

    const transNodes = new vis.DataSet({json.dumps(trans_nodes)});
    const transEdges = new vis.DataSet({json.dumps(trans_edges)});
    const transNetwork = new vis.Network(
        document.getElementById("trans"),
        {{ nodes: transNodes, edges: transEdges }},
        {{ physics: {{ enabled: true }}, edges: {{ arrows: "to" }} }}
    );

    const stateColors = {json.dumps(state_colors)};

    transNetwork.on("click", function (params) {{
        if (params.nodes.length > 0) {{
            const clickedId = params.nodes[0];
            const colorUpdate = stateColors[clickedId];
            if (colorUpdate) {{
                baseNodes.update(colorUpdate);
            }}
        }}
    }});
  </script>
</body>
</html>
"""
    with open(filename, "w") as f:
        f.write(html_template)
    return filename