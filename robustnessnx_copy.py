# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 21:44:53 2020

@author: BALAMLAPTOP2
"""
def nodestoremove(p, q, nodes_0_c, nodes_1_c, df_atc, df_icd, type_method):
    if type_method == 0:
        atcNodesToRemove = random.sample(nodes_1_c, int(round(p * len(nodes_1_c))))
        icdNodesToRemove = random.sample(nodes_0_c, int(round(q * len(nodes_0_c))))
    elif type_method == 1:
        atcNodesToRemove = list(df_atc.head(int(round(p * len(nodes_1_c))))['node'])
        icdNodesToRemove = list(df_icd.head(int(round(p * len(nodes_0_c))))['node'])
    
    return icdNodesToRemove, atcNodesToRemove
        
for i in range(1,1000):
    index_p = 0
    index_q = 0
    for p in p_vector:
        for q in p_vector:
            unfrozen_graph = nx.Graph(C)
            icdNodesToRemove, atcNodesToRemove = nodestoremove(p, q, nodes_0_c, nodes_1_c, df_atc, df_icd, type_method)
            unfrozen_graph.remove_nodes_from(icdNodesToRemove)
            unfrozen_graph.remove_nodes_from(atcNodesToRemove)
            if type_proj == 0:
                nodes_lst = [x for x in nodes_0_c if x not in icdNodesToRemove]
            elif type_proj == 1:
                nodes_lst = [x for x in nodes_1_c if x not in atcNodesToRemove]
            GP = bipartite.projected_graph(unfrozen_graph, nodes_lst)
            components = sorted(nx.connected_components(GP), key=len, reverse=True)
            nodes_connected = sum(list(map(lambda c: len(c), components)))
            mean_size_components = nodes_connected / len(components)
            nodes_unconnected = len(nodes_0_c) - nodes_connected
            degrees = GP.degree()
            sum_of_edges = sum(list(dict(degrees).values()))
            avg_degree = sum_of_edges / GP.number_of_nodes()
            avg_clustering = nx.average_clustering(GP)
            if i == 1:
                np_avg_degree[index_p][index_q] = avg_degree
                np_conn_nodes[index_p][index_q] = nodes_connected
                np_unconn_nodes[index_p][index_q] = nodes_unconnected
                np_conn_components[index_p][index_q] = len(components)
                np_mean_size[index_p][index_q] = mean_size_components
                np_clustering[index_p][index_q] = avg_clustering
            else:
                np_avg_degree[index_p][index_q] = np_avg_degree[index_p][index_q] + avg_degree
                np_conn_nodes[index_p][index_q] = np_conn_nodes[index_p][index_q] + nodes_connected
                np_unconn_nodes[index_p][index_q] = np_unconn_nodes[index_p][index_q] + nodes_unconnected
                np_conn_components[index_p][index_q] = np_conn_components[index_p][index_q] + len(components)
                np_mean_size[index_p][index_q] = np_mean_size[index_p][index_q] + mean_size_components
                np_clustering[index_p][index_q] = np_clustering[index_p][index_q] + avg_clustering
            index_q += 1
        index_q = 0
        index_p += 1

np_avg_degree = np.true_divide(np_avg_degree, 1000)
np_conn_nodes = np.true_divide(np_conn_nodes, 1000)
np_unconn_nodes = np.true_divide(np_unconn_nodes, 1000)
np_conn_components = np.true_divide(np_conn_components, 1000)
np_mean_size = np.true_divide(np_mean_size, 1000)
np_clustering = np.true_divide(np_clustering, 1000)
            
for i in range(1,1000):
    index_p = 0
    index_q = 0
    if type_method == 0:
        if type_proj == 0: # Disease projection
            print("Calculating randomly robustness disease projection ...")
            for p in p_vector:
                for q in p_vector:
                    unfrozen_graph = nx.Graph(C)
                    atcNodesToRemove = random.sample(nodes_1_c, int(round(p * len(nodes_1_c))))
                    icdNodesToRemove = random.sample(nodes_0_c, int(round(q * len(nodes_0_c))))
                    unfrozen_graph.remove_nodes_from(icdNodesToRemove)
                    unfrozen_graph.remove_nodes_from(atcNodesToRemove)
                    icd_lst = [x for x in nodes_0_c if x not in icdNodesToRemove]
                    GPCIE = bipartite.projected_graph(unfrozen_graph, icd_lst)
                    components = sorted(nx.connected_components(GPCIE), key=len, reverse=True)
                    nodes_connected = sum(list(map(lambda c: len(c), components)))
                    mean_size_components = nodes_connected / len(components)
                    nodes_unconnected = len(nodes_0_c) - nodes_connected
                    degrees = GPCIE.degree()
                    sum_of_edges = sum(list(dict(degrees).values()))
                    avg_degree = sum_of_edges / GPCIE.number_of_nodes()
                    avg_clustering = nx.average_clustering(GPCIE)
                    if i == 1:
                        np_avg_degree[index_p][index_q] = avg_degree
                        np_conn_nodes[index_p][index_q] = nodes_connected
                        np_unconn_nodes[index_p][index_q] = nodes_unconnected
                        np_conn_components[index_p][index_q] = len(components)
                        np_mean_size[index_p][index_q] = mean_size_components
                        np_clustering[index_p][index_q] = avg_clustering
                    else:
                        np_avg_degree[index_p][index_q] = np_avg_degree[index_p][index_q] + avg_degree
                        np_conn_nodes[index_p][index_q] = np_conn_nodes[index_p][index_q] + nodes_connected
                        np_unconn_nodes[index_p][index_q] = np_unconn_nodes[index_p][index_q] + nodes_unconnected
                        np_conn_components[index_p][index_q] = np_conn_components[index_p][index_q] + len(components)
                        np_mean_size[index_p][index_q] = np_mean_size[index_p][index_q] + mean_size_components
                        np_clustering[index_p][index_q] = np_clustering[index_p][index_q] + avg_clustering
                    index_q += 1
                index_q = 0
                index_p += 1
        elif type_proj == 1: # Active Ingredient projection
            print("Calculating randomly robustness active ingredients projection ...")
            for p in p_vector:
                for q in p_vector:
                    unfrozen_graph = nx.Graph(C)
                    icdNodesToRemove = random.sample(nodes_0_c, int(round(p * len(nodes_0_c))))
                    atcNodesToRemove = random.sample(nodes_1_c, int(round(q * len(nodes_1_c))))
                    unfrozen_graph.remove_nodes_from(icdNodesToRemove)
                    unfrozen_graph.remove_nodes_from(atcNodesToRemove)
                    atc_lst = [x for x in nodes_1_c if x not in atcNodesToRemove]
                    GPATC = bipartite.projected_graph(unfrozen_graph, atc_lst)
                    components = sorted(nx.connected_components(GPATC), key=len, reverse=True)
                    nodes_connected = sum(list(map(lambda c: len(c), components)))
                    mean_size_components = nodes_connected / len(components)
                    nodes_unconnected = len(nodes_1_c) - nodes_connected
                    degrees = GPATC.degree()
                    sum_of_edges = sum(list(dict(degrees).values()))
                    avg_degree = sum_of_edges / GPATC.number_of_nodes()
                    avg_clustering = nx.average_clustering(GPATC)
                    np_avg_degree[index_p][index_q] = avg_degree
                    np_conn_nodes[index_p][index_q] = nodes_connected
                    np_unconn_nodes[index_p][index_q] = nodes_unconnected
                    np_conn_components[index_p][index_q] = len(components)
                    np_mean_size[index_p][index_q] = mean_size_components
                    np_clustering[index_p][index_q] = avg_clustering
                    index_q += 1
                index_q = 0
                index_p += 1
    elif type_method == 1:    
        # Dirigido
        if type_proj == 0: # Disease projection
            print("Calculating directed robustness disease projection ...")
            for p in p_vector:
    #            atcNodesToRemove = random.sample(nodes_1_c, int(round(p * len(nodes_1_c))))
                for q in p_vector:
                    unfrozen_graph = nx.Graph(C)
    #                icdNodesToRemove = random.sample(nodes_0_c, int(round(q * len(nodes_0_c))))
                    atcNodesToRemove = list(df_atc.head(int(round(p * len(nodes_1_c))))['node'])
                    icdNodesToRemove = list(df_icd.head(int(round(p * len(nodes_0_c))))['node'])
                    unfrozen_graph.remove_nodes_from(icdNodesToRemove)
                    unfrozen_graph.remove_nodes_from(atcNodesToRemove)
                    icd_lst = [x for x in nodes_0_c if x not in icdNodesToRemove]
                    GPCIE = bipartite.projected_graph(unfrozen_graph, icd_lst)
                    components = sorted(nx.connected_components(GPCIE), key=len, reverse=True)
                    nodes_connected = sum(list(map(lambda c: len(c), components)))
                    mean_size_components = nodes_connected / len(components)
                    nodes_unconnected = len(nodes_0_c) - nodes_connected
                    degrees = GPCIE.degree()
                    sum_of_edges = sum(list(dict(degrees).values()))
                    avg_degree = sum_of_edges / GPCIE.number_of_nodes()
                    avg_clustering = nx.average_clustering(GPCIE)
                    np_avg_degree[index_p][index_q] = avg_degree
                    np_conn_nodes[index_p][index_q] = nodes_connected
                    np_unconn_nodes[index_p][index_q] = nodes_unconnected
                    np_conn_components[index_p][index_q] = len(components)
                    np_mean_size[index_p][index_q] = mean_size_components
                    np_clustering[index_p][index_q] = avg_clustering
                    index_q += 1
                index_q = 0
                index_p += 1
        elif type_proj == 1: # Active Ingredient projection
            print("Calculating directed robustness active ingredients projection ...")
            for p in p_vector:
                for q in p_vector:
                    unfrozen_graph = nx.Graph(C)
    #                atcNodesToRemove = random.sample(nodes_1_c, int(round(q * len(nodes_1_c))))
                    icdNodesToRemove = list(df_icd.head(int(round(p * len(nodes_0_c))))['node'])
                    atcNodesToRemove = list(df_atc.head(int(round(p * len(nodes_1_c))))['node'])
                    unfrozen_graph.remove_nodes_from(icdNodesToRemove)
                    unfrozen_graph.remove_nodes_from(atcNodesToRemove)
                    atc_lst = [x for x in nodes_1_c if x not in atcNodesToRemove]
                    GPATC = bipartite.projected_graph(unfrozen_graph, atc_lst)
                    components = sorted(nx.connected_components(GPATC), key=len, reverse=True)
                    nodes_connected = sum(list(map(lambda c: len(c), components)))
                    mean_size_components = nodes_connected / len(components)
                    nodes_unconnected = len(nodes_1_c) - nodes_connected
                    degrees = GPATC.degree()
                    sum_of_edges = sum(list(dict(degrees).values()))
                    avg_degree = sum_of_edges / GPATC.number_of_nodes()
                    avg_clustering = nx.average_clustering(GPATC)
                    
                    np_avg_degree[index_p][index_q] = avg_degree
                    np_conn_nodes[index_p][index_q] = nodes_connected
                    np_unconn_nodes[index_p][index_q] = nodes_unconnected
                    np_conn_components[index_p][index_q] = len(components)
                    np_mean_size[index_p][index_q] = mean_size_components
                    np_clustering[index_p][index_q] = avg_clustering
                    index_q += 1
                index_q = 0
                index_p += 1