import sys
import graph_tool as gt
import graph_tool.topology as gtt
import graph_tool.stats as gts
import graph_tool.util as gtu
import itertools as it
from find_bridges import find_bridges

# =============================================================================
def print_usage():
    print(' ')
    print('        usage: python3 must.py network_file seed_file hub_penalty num_trees tolerance outfile_name(optional)')
    print('        -----------------------------------------------------------------')
    print('        network_file       : The network file must be provided as graphml or gt format.')
    print('                             Nodes attribute "primaryDomainId" is used in graphml as node id.')
    print('        seed_file          : Table containing the seed proteins (if table contains')
    print('                             more than one column they must be tab-separated;')
    print('                             the first column will be used only)')
    print('        hub_penalty        : A float value between 0.0-1.0 to be used as Penalty parameter for hubs.')
    print('                             Set edge weight to 1 + (hub_penalty / 2) (e.source.degree + e.target.degree).')
    print('        num_trees          : An integer value between 0-25. Number of Steiner trees to return')
    print('                             Notice that the algorithm might not be able to reach this number of Steiner trees.')
    print('        tolerance          : A float value between 0-20. The error tolerance of the subsequent Steiner trees')
    print('                             w.r.t. the first one in percent.')
    print('        outfile_name       : Results will be saved under this file name')
    print('                             by default the outfile_name is set to "must_result.txt"')
    print(' ')


# =============================================================================
def check_input_style(input_list):
    try:
        network_file = input_list[1]
        seeds_file = input_list[2]
        hub_penalty = float(input_list[3])
        num_trees = int(input_list[4])
        tolerance = float(input_list[5])
        outfile_name = 'must_hp_%(hp).2f_nt_%(nt)d_tlr_%(tl).0f.txt' % {"hp": hub_penalty, "nt": num_trees, "tl": tolerance}
    # if no input is given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)

    if len(input_list) == 7:
        try:
            outfile_name = input_list[6]
        except:
            print_usage()
            sys.exit(0)

    return network_file, seeds_file, hub_penalty, num_trees, tolerance, outfile_name


# =============================================================================
def read_input(network_file, seed_file):
    """
    Reads the network and the list of seed proteins from external files.
    * The network must be provided as a graphml or gt format file.
    Nodes should have "primaryDomainId" attribute as their unique identifier.
    * The seed proteins mus be provided as a table. If the table has more
    than one column, they must be tab-separated. The first column will
    be used only.
    * Lines that start with '#' will be ignored in both cases
    """

    # p_type = "Protein"
    d_type = "Drug"
    # d_type = "drug-approved"
    # pp_type = "ProteinInteractsWithProtein"
    dp_type = "DrugHasTarget"
    node_name_attribute = "name"
    # node_name_attribute = "name"

    # read the seed proteins:
    seeds = list()
    for line in open(seed_file, 'r'):
        # lines starting with '#' will be ignored
        if line[0] == '#':
            continue
        # the first column in the line will be interpreted as a seed
        # protein:
        line_data = line.strip().split('\t')
        seed = line_data[0]
        seeds.append(seed)

    seeds.sort()
    # read the network:
    G = gt.load_graph_from_csv(network_file, csv_options={"delimiter":"\t"})
    G_lcc = gtt.extract_largest_component(G, directed=False, prune=True)

    seed_ids = []
    is_matched = {protein: False for protein in seeds}
    print("Number of nodes in G {} and in G_lcc {}".format(G.num_vertices(), G_lcc.num_vertices()))
    print(G_lcc)
    # obtain seeds and drugs ids
    for node in range(G_lcc.num_vertices()):
        if G_lcc.vertex_properties[node_name_attribute][node] in seeds:
            seed_ids.append(node)
            is_matched[G_lcc.vertex_properties[node_name_attribute][node]] = True

    # Check that all seed seeds have been matched and print the ones not found.
    for protein, found in is_matched.items():
        if not found:
            print("No node named {} found in the network read from {}".format(protein, network_file))
            # raise ValueError("Invalid seed protein {}. No node named {} in {}.".format(protein, protein, network_file))

    return G_lcc, seeds, seed_ids

# =============================================================================
def edge_weights(g, hub_penalty, inverse=False):
    avdeg = gts.vertex_average(g, "total")[0]
    weights = g.new_edge_property("double", val=avdeg)
    if hub_penalty <= 0:
        return weights
    if hub_penalty > 1:
        raise ValueError("Invalid hub penalty {}.".format(hub_penalty))
    for e in g.edges():
        edge_avdeg = float(e.source().out_degree() + e.target().out_degree()) / 2.0
        penalized_weight = (1.0 - hub_penalty) * avdeg + hub_penalty * edge_avdeg
        if inverse:
            weights[e] = 1.0 / penalized_weight
        else:
            weights[e] = penalized_weight
    return weights

# =============================================================================
def steiner_tree(g, seeds, seed_map, weights, non_zero_hub_penalty):

    node_name_attribute = "name" # nodes in the input network which is created from RepoTrialDB have primaryDomainId as name attribute
    mc = gt.Graph(directed=False)
    eprop_dist = mc.new_edge_property("int")
    mc.ep['dist'] = eprop_dist
    vprop_name = mc.new_vertex_property("string")
    mc.vp[node_name_attribute] = vprop_name

    eprop_path = mc.new_edge_property("object")
    mc.ep['path'] = eprop_path

    mc_vertex_map = dict()
    mc_id_map = dict()
    for i in range(len(seeds)):
        vert = mc.add_vertex()
        vprop_name[i] = seeds[i]
        mc_vertex_map[seeds[i]] = vert
        mc_id_map[vert] = i

    for u, v in it.combinations(seeds, 2):
        _, elist = gtt.shortest_path(g, g.vertex(seed_map[u]), g.vertex(seed_map[v]), weights=weights, negative_weights=False, pred_map=None, dag=False)
        e = mc.add_edge(mc_vertex_map[u], mc_vertex_map[v])
        eprop_dist[e] = len(elist)
        mc.ep.path[e] = list(elist)

    mst = gtt.min_spanning_tree(mc, weights=eprop_dist, root=None, tree_map=None)
    mc.set_edge_filter(mst)

    g2 = gt.Graph(directed=False)
    vprop_name = g2.new_vertex_property("string")
    g2.vp[node_name_attribute] = vprop_name

    g2_vertex_map = dict()
    g2_id_map = dict()
    addedNodes = set()
    for i in range(len(seeds)):
        vert = g2.add_vertex()
        vprop_name[i] = seeds[i]
        g2_vertex_map[seeds[i]] = vert
        g2_id_map[vert] = i
        addedNodes.add(seeds[i])

    allmcedges = []

    for mc_edges in mc.edges():
        path = mc.ep.path[mc_edges]
        allmcedges.extend(path)

    j = len(seeds)
    allmcedges_g2 = []
    for e in allmcedges:
        # sourceName = g.vertex_properties["name"][e.source()]
        # targetName = g.vertex_properties["name"][e.target()]
        sourceName = g.vertex_properties[node_name_attribute][e.source()]
        targetName = g.vertex_properties[node_name_attribute][e.target()]
        if sourceName not in addedNodes:
            vert = g2.add_vertex()
            vprop_name[j] = sourceName
            g2_vertex_map[sourceName] = vert
            g2_id_map[vert] = j
            addedNodes.add(sourceName)
            j += 1
        if targetName not in addedNodes:
            vert = g2.add_vertex()
            vprop_name[j] = targetName
            g2_vertex_map[targetName] = vert
            g2_id_map[vert] = j
            addedNodes.add(targetName)
            j += 1
        allmcedges_g2.append(g2.add_edge(g2_vertex_map[sourceName], g2_vertex_map[targetName]))
    weights_g2 = g2.new_edge_property("double", val=1.0)
    if non_zero_hub_penalty:
        for e, e_g2 in zip(allmcedges, allmcedges_g2):
            weights_g2[e_g2] = weights[e]
    mst2 = gtt.min_spanning_tree(g2, root=None, tree_map=None, weights=weights_g2)
    g2.set_edge_filter(mst2)
    # vw = gt.GraphView(g2, efilt=mst2)
    # g3 = Graph(vw, prune=True)

    # g3 = Graph(g22)

    while True:
        noneSteinerLeaves = []
        for i in range(g2.num_vertices()):
            if g2.vertex(i).out_degree() == 1 and g2.vertex_properties[node_name_attribute][i] not in seeds:
                noneSteinerLeaves.append(i)
        if len(noneSteinerLeaves) == 0:
            break
        noneSteinerLeaves = reversed(sorted(noneSteinerLeaves))
        for node in noneSteinerLeaves:
            # outarray = g3.get_out_edges(node)
            g2.remove_edge(g2.edge(g2.vertex(node), g2.get_all_neighbors(node)[0]))
            g2.remove_vertex(node)

    return g2

# =============================================================================
def must(G_lcc, seed_ids, hub_penalty, num_trees, tolerance, outfile=None):
    """
    Runs the TrustRank w.r.t. seeds
    Input:
    ------
     - G:
        The network
     - seed_ids:
            a set of seed proteins
     - outfile:
              filename for the output generates by the algorithm,
              if not given the program will assign a name based on the input parameters
     - hub_penalty:
                a float value between 0.0-1.0 to be used as Penalty parameter for hubs.
                Set edge weight to 1 + (hub_penalty / 2) (e.source.degree + e.target.degree
    - num_trees:
            an integer value between 0-25. Number of Steiner trees to return.
            Notice that the algorithm might not be able to reach this number of Steiner trees.
    - tolerance:
            a float value between 0-20. The error tolerance of the subsequent Steiner trees w.r.t. the first one in
            percent.
     Returns:
     --------
      - score_map: A map of scores for all drugs in the network
        A sorted map based on the score values is also written in the outfile

      Notes
    -----
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.

    """

    # Set number of threads if OpenMP support is enabled.
    # if gt.openmp_enabled():
    #     gt.openmp_set_num_threads(num_threads)

    node_name_attribute = "name"
    # node_name_attribute = "name"

    seed_map = {G_lcc.vertex_properties[node_name_attribute][node]: node for node in seed_ids}
    weights = edge_weights(G_lcc, hub_penalty)

    # Find first steiner trees
    first_tree = steiner_tree(G_lcc, seeds, seed_map, weights, hub_penalty > 0)
    num_found_trees = 1
    tree_edges = []
    for tree_edge in first_tree.edges():
        source_name = first_tree.vertex_properties[node_name_attribute][first_tree.vertex_index[tree_edge.source()]]
        target_name = first_tree.vertex_properties[node_name_attribute][first_tree.vertex_index[tree_edge.target()]]
        tree_edges.append((gtu.find_vertex(G_lcc, prop=G_lcc.vertex_properties[node_name_attribute], match=source_name)[0], gtu.find_vertex(G_lcc, prop=G_lcc.vertex_properties[node_name_attribute], match=target_name)[0]))
    cost_first_tree = sum([weights[G_lcc.edge(source, target)] for source, target in tree_edges])
    returned_nodes = set(int(gtu.find_vertex(G_lcc, prop=G_lcc.vertex_properties[node_name_attribute], match=first_tree.vertex_properties[node_name_attribute][node])[0]) for node in range(first_tree.num_vertices()))

    if num_trees > 1:
        is_bridge = find_bridges(G_lcc)
        edge_filter = G_lcc.new_edge_property("boolean", True)
        found_new_tree = True
        while len(tree_edges) > 0:
            found_new_tree = False
            tree_edge = tree_edges.pop()
            g_edge = G_lcc.edge(tree_edge[0], tree_edge[1])
            if not is_bridge[g_edge]:
                edge_filter[g_edge] = False
                G_lcc.set_edge_filter(edge_filter)
                next_tree = steiner_tree(G_lcc, seeds, seed_map, weights, hub_penalty > 0)
                next_tree_edges = set()
                for next_tree_edge in next_tree.edges():
                    source_name = next_tree.vertex_properties[node_name_attribute][next_tree.vertex_index[next_tree_edge.source()]]
                    target_name = next_tree.vertex_properties[node_name_attribute][next_tree.vertex_index[next_tree_edge.target()]]
                    next_tree_edges.add((gtu.find_vertex(G_lcc, prop=G_lcc.vertex_properties[node_name_attribute], match=source_name)[0],gtu.find_vertex(G_lcc, prop=G_lcc.vertex_properties[node_name_attribute], match=target_name)[0]))
                cost_next_tree = sum([weights[G_lcc.edge(source, target)] for source, target in next_tree_edges])
                if cost_next_tree <= cost_first_tree * ((100.0 + tolerance) / 100.0):
                    found_new_tree = True
                    num_found_trees += 1
                    for node in range(next_tree.num_vertices()):
                        returned_nodes.add(int(gtu.find_vertex(G_lcc, prop=G_lcc.vertex_properties[node_name_attribute],match=next_tree.vertex_properties[node_name_attribute][node])[0]))
                    removed_edges = []
                    for source, target in tree_edges:
                        if not ((source, target) in set(next_tree_edges)) or ((target, source) in set(next_tree_edges)):
                            removed_edges.append((source, target))
                    for edge in removed_edges:
                        tree_edges.remove(edge)
                G_lcc.clear_filters()
                edge_filter[g_edge] = True
            if num_found_trees >= num_trees:
                break

    returned_edges = []
    for node in returned_nodes:
        for neighbor in G_lcc.get_all_neighbors(node):
            if int(neighbor) > node and int(neighbor) in returned_nodes:
                returned_edges.append((node, int(neighbor)))

    subgraph = {"nodes": [G_lcc.vertex_properties[node_name_attribute][node] for node in returned_nodes],
                "edges": [{"from": G_lcc.vertex_properties[node_name_attribute][source], "to": G_lcc.vertex_properties[node_name_attribute][target]} for
                          source, target in returned_edges]}
    is_seed = {G_lcc.vertex_properties[node_name_attribute][node]: node in set(seed_ids) for node in returned_nodes}

    with open(outfile, 'w') as fout:
        fout.write('\t'.join(['source_node', 'target_node']) + '\n')
        for s, t in returned_edges:
            fout.write('\t'.join([G_lcc.vertex_properties[node_name_attribute][s], G_lcc.vertex_properties[node_name_attribute][t]]) + '\n')

    return subgraph

# ===========================================================================

if __name__ == '__main__':
    # -----------------------------------------------------
    # Checking for input from the command line:
    # -----------------------------------------------------
    #
    # [1] file providing the network in graphml or gt format
    # [2] file with the seed proteins (if table contains more than one
    #     column they must be tab-separated; the first column will be
    #     used only)
    # [3] float value between 0.0-1.0 to be used as hub penalty parameter
    # [4] integer value between 0-25 to be used as number of trees parameter
    # [5] float value between 0-20 to be used as tolerance parameter
    # [6] (optional) name for the results file

    # check if input style is correct
    input_list = sys.argv
    network_file, seeds_file, hub_penalty, num_trees, tolerance, outfile_name = check_input_style(input_list)
    print('Network file: ' + network_file)
    print('Seed file:' + seeds_file)
    print('Hub penalty: ' + str(hub_penalty))
    print('Number of trees: ' + str(num_trees))
    print('Tolerance: ' + str(tolerance))
    print('Output file: ' + outfile_name)

    # read the network and the seed proteins:
    G_lcc, seeds, seed_ids = read_input(network_file, seeds_file)
    # kick out terminals not in graph
    seeds = list(set(seeds).intersection(set(G_lcc.vertex_properties['name'])))
    # run must
    result_subgraph = must(G_lcc, seed_ids, hub_penalty, num_trees, tolerance, outfile=outfile_name)

    print("\n results have been saved to '%s' \n" % outfile_name)
    print(result_subgraph)