# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 01:35:18 2020

@author: Irene López-Rodríguez

Input variables: 
    type_method = [0 : random, 1 : degree]
    projection = [0 : ICD, 1 : ATC]
"""

import sys
import networkx as nx
from networkx.algorithms import bipartite
import collections
import pandas as pd
import random
import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt

def nodestoremove(p, q, nodes_0_c, nodes_1_c, df_atc, df_icd, type_method):
    if type_method == 0:
        atcNodesToRemove = random.sample(nodes_1_c, int(round(p * len(nodes_1_c))))
        icdNodesToRemove = random.sample(nodes_0_c, int(round(q * len(nodes_0_c))))
    elif type_method == 1:
        atcNodesToRemove = list(df_atc.head(int(round(p * len(nodes_1_c))))['node'])
        icdNodesToRemove = list(df_icd.head(int(round(p * len(nodes_0_c))))['node'])
    
    return icdNodesToRemove, atcNodesToRemove


if __name__ == '__main__':
    # Remover al azar 0 o dirigido 1
    # Proyección ICD 0 o ATC 1 
    type_method = int(sys.argv[1])
    type_proj = int(sys.argv[2])
    
    print("Reading file ...")

    vdmdata = pd.read_csv('vdmdata_reduce.csv', encoding = 'utf-8-sig')
    
    nodes_0 = []
    nodes_1 = []
    for m in vdmdata.iterrows():
        nodes_0.append(m[1][0]) #ICD
        nodes_1.append(m[1][1]) #ATC
        
    nodes_0 = list(dict.fromkeys(nodes_0))
    nodes_1 = list(dict.fromkeys(nodes_1))
    
    print("Building a bipartite graph ...")
    # Build a bipartite graph:
    G = nx.Graph()
    # Add nodes ATC - ICD
    G.add_nodes_from(nodes_0, bipartite=0) # Add the node attribute “bipartite” disease
    G.add_nodes_from(nodes_1, bipartite=1) # active substance
    
    # Add edges without weight
    for m in vdmdata.iterrows():
        enfermedad = m[1][0];
        #peso = m[1][3];
        sustancia = m[1][1];
        G.add_edge(enfermedad, sustancia)
        
    print("Getting largest component ...")
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    largest_component = components[0]
    C = G.subgraph(largest_component)
        
    degX,degY=bipartite.degrees(C,nodes_0)
    degATC = dict(degX).values()
    degCIE = dict(degY).values()
    counterATC = collections.Counter(degATC)
    counterCIE = collections.Counter(degCIE)
    
    
    df_icd = pd.DataFrame(dict(degY).items(), columns=['node', 'degree'])
    df_atc = pd.DataFrame(dict(degX).items(), columns=['node', 'degree'])
    
    df_icd = df_icd.sort_values(by=['degree', 'node'], ascending=False)
    df_atc = df_atc.sort_values(by=['degree', 'node'], ascending=False)
    
    nodes_0_c = []
    nodes_1_c = []
    for n in C.nodes(data=True):
        if n[1]['bipartite'] == 0:
            nodes_0_c.append(n[0])
        if n[1]['bipartite'] == 1:
            nodes_1_c.append(n[0])
            
#    unfrozen_graph = nx.Graph(C)
    
    print("Defining variables ...")
#    p_vector = np.arange(0,1,0.1)
#    p_vector = [ round(x * 0.001, 3) for x in range(0, 1000)]
    p_vector = [ round(x * 0.01, 2) for x in range(0, 100)]
#    [ round(x * 0.01, 2) for x in range(0, 100)]
    
    np_avg_degree = np.zeros((len(p_vector), len(p_vector)))
    np_conn_nodes = np.zeros((len(p_vector), len(p_vector)))
    np_unconn_nodes = np.zeros((len(p_vector), len(p_vector)))
    np_conn_components = np.zeros((len(p_vector), len(p_vector)))
    np_mean_size = np.zeros((len(p_vector), len(p_vector)))
#    np_clustering = np.zeros((len(p_vector), len(p_vector)))
#    type_method = 0
#    type_proj = 1
    for i in range(1,11):
        print("Runing iteration ... "+str(i))
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
#                avg_clustering = nx.average_clustering(GP)
                if i == 1:
                    np_avg_degree[index_p][index_q] = avg_degree
                    np_conn_nodes[index_p][index_q] = nodes_connected
                    np_unconn_nodes[index_p][index_q] = nodes_unconnected
                    np_conn_components[index_p][index_q] = len(components)
                    np_mean_size[index_p][index_q] = mean_size_components
#                    np_clustering[index_p][index_q] = avg_clustering
                else:
                    np_avg_degree[index_p][index_q] = np_avg_degree[index_p][index_q] + avg_degree
                    np_conn_nodes[index_p][index_q] = np_conn_nodes[index_p][index_q] + nodes_connected
                    np_unconn_nodes[index_p][index_q] = np_unconn_nodes[index_p][index_q] + nodes_unconnected
                    np_conn_components[index_p][index_q] = np_conn_components[index_p][index_q] + len(components)
                    np_mean_size[index_p][index_q] = np_mean_size[index_p][index_q] + mean_size_components
#                    np_clustering[index_p][index_q] = np_clustering[index_p][index_q] + avg_clustering
                index_q += 1
            index_q = 0
            index_p += 1

    np_avg_degree = np.true_divide(np_avg_degree, i)
    np_conn_nodes = np.true_divide(np_conn_nodes, i)
    np_unconn_nodes = np.true_divide(np_unconn_nodes, i)
    np_conn_components = np.true_divide(np_conn_components, i)
    np_mean_size = np.true_divide(np_mean_size, i)
#    np_clustering = np.true_divide(np_clustering, i)
    
    print("Creating files ... ")
    df = pd.DataFrame(np_avg_degree, index=p_vector, columns=p_vector)
    df.to_csv('np_avg_degree_'+str(type_proj)+'_'+str(type_method)+'.csv', index=True, header=True, sep=',', encoding = 'utf-8-sig')
    
    df = pd.DataFrame(np_conn_nodes, index=p_vector, columns=p_vector)
    df.to_csv('np_conn_nodes_'+str(type_proj)+'_'+str(type_method)+'.csv', index=True, header=True, sep=',', encoding = 'utf-8-sig')
    
    df = pd.DataFrame(np_unconn_nodes, index=p_vector, columns=p_vector)
    df.to_csv('np_unconn_nodes_'+str(type_proj)+'_'+str(type_method)+'.csv', index=True, header=True, sep=',', encoding = 'utf-8-sig')
    
    df = pd.DataFrame(np_conn_components, index=p_vector, columns=p_vector)
    df.to_csv('np_conn_components_'+str(type_proj)+'_'+str(type_method)+'.csv', index=True, header=True, sep=',', encoding = 'utf-8-sig')
    
    df = pd.DataFrame(np_mean_size, index=p_vector, columns=p_vector)
    df.to_csv('np_mean_size_'+str(type_proj)+'_'+str(type_method)+'.csv', index=True, header=True, sep=',', encoding = 'utf-8-sig')
    
#    df = pd.DataFrame(np_clustering, index=p_vector, columns=p_vector)
#    df.to_csv('np_clustering_'+str(type_proj)+'_'+str(type_method)+'.csv', index=True, header=True, sep=',', encoding = 'utf-8-sig')


#np.savetxt('distancia_centroides.out', geodis_centroids, delimiter=',')


    # Randomly
    
#    for i in range(1,1000):
#        index_p = 0
#        index_q = 0
#        if type_method == 0:
#            if type_proj == 0: # Disease projection
#                print("Calculating randomly robustness disease projection ...")
#                for p in p_vector:
#                    for q in p_vector:
#                        unfrozen_graph = nx.Graph(C)
#                        atcNodesToRemove = random.sample(nodes_1_c, int(round(p * len(nodes_1_c))))
#                        icdNodesToRemove = random.sample(nodes_0_c, int(round(q * len(nodes_0_c))))
#                        unfrozen_graph.remove_nodes_from(icdNodesToRemove)
#                        unfrozen_graph.remove_nodes_from(atcNodesToRemove)
#                        icd_lst = [x for x in nodes_0_c if x not in icdNodesToRemove]
#                        GPCIE = bipartite.projected_graph(unfrozen_graph, icd_lst)
#                        components = sorted(nx.connected_components(GPCIE), key=len, reverse=True)
#                        nodes_connected = sum(list(map(lambda c: len(c), components)))
#                        mean_size_components = nodes_connected / len(components)
#                        nodes_unconnected = len(nodes_0_c) - nodes_connected
#                        degrees = GPCIE.degree()
#                        sum_of_edges = sum(list(dict(degrees).values()))
#                        avg_degree = sum_of_edges / GPCIE.number_of_nodes()
#                        avg_clustering = nx.average_clustering(GPCIE)
#                        if i == 1:
#                            np_avg_degree[index_p][index_q] = avg_degree
#                            np_conn_nodes[index_p][index_q] = nodes_connected
#                            np_unconn_nodes[index_p][index_q] = nodes_unconnected
#                            np_conn_components[index_p][index_q] = len(components)
#                            np_mean_size[index_p][index_q] = mean_size_components
#                            np_clustering[index_p][index_q] = avg_clustering
#                        else:
#                            np_avg_degree[index_p][index_q] = np_avg_degree[index_p][index_q] + avg_degree
#                            np_conn_nodes[index_p][index_q] = np_conn_nodes[index_p][index_q] + nodes_connected
#                            np_unconn_nodes[index_p][index_q] = np_unconn_nodes[index_p][index_q] + nodes_unconnected
#                            np_conn_components[index_p][index_q] = np_conn_components[index_p][index_q] + len(components)
#                            np_mean_size[index_p][index_q] = np_mean_size[index_p][index_q] + mean_size_components
#                            np_clustering[index_p][index_q] = np_clustering[index_p][index_q] + avg_clustering
#                        index_q += 1
#                    index_q = 0
#                    index_p += 1
#            elif type_proj == 1: # Active Ingredient projection
#                print("Calculating randomly robustness active ingredients projection ...")
#                for p in p_vector:
#                    for q in p_vector:
#                        unfrozen_graph = nx.Graph(C)
#                        icdNodesToRemove = random.sample(nodes_0_c, int(round(p * len(nodes_0_c))))
#                        atcNodesToRemove = random.sample(nodes_1_c, int(round(q * len(nodes_1_c))))
#                        unfrozen_graph.remove_nodes_from(icdNodesToRemove)
#                        unfrozen_graph.remove_nodes_from(atcNodesToRemove)
#                        atc_lst = [x for x in nodes_1_c if x not in atcNodesToRemove]
#                        GPATC = bipartite.projected_graph(unfrozen_graph, atc_lst)
#                        components = sorted(nx.connected_components(GPATC), key=len, reverse=True)
#                        nodes_connected = sum(list(map(lambda c: len(c), components)))
#                        mean_size_components = nodes_connected / len(components)
#                        nodes_unconnected = len(nodes_1_c) - nodes_connected
#                        degrees = GPATC.degree()
#                        sum_of_edges = sum(list(dict(degrees).values()))
#                        avg_degree = sum_of_edges / GPATC.number_of_nodes()
#                        avg_clustering = nx.average_clustering(GPATC)
#                        np_avg_degree[index_p][index_q] = avg_degree
#                        np_conn_nodes[index_p][index_q] = nodes_connected
#                        np_unconn_nodes[index_p][index_q] = nodes_unconnected
#                        np_conn_components[index_p][index_q] = len(components)
#                        np_mean_size[index_p][index_q] = mean_size_components
#                        np_clustering[index_p][index_q] = avg_clustering
#                        index_q += 1
#                    index_q = 0
#                    index_p += 1
#        elif type_method == 1:    
#            # Dirigido
#            if type_proj == 0: # Disease projection
#                print("Calculating directed robustness disease projection ...")
#                for p in p_vector:
#        #            atcNodesToRemove = random.sample(nodes_1_c, int(round(p * len(nodes_1_c))))
#                    for q in p_vector:
#                        unfrozen_graph = nx.Graph(C)
#        #                icdNodesToRemove = random.sample(nodes_0_c, int(round(q * len(nodes_0_c))))
#                        atcNodesToRemove = list(df_atc.head(int(round(p * len(nodes_1_c))))['node'])
#                        icdNodesToRemove = list(df_icd.head(int(round(p * len(nodes_0_c))))['node'])
#                        unfrozen_graph.remove_nodes_from(icdNodesToRemove)
#                        unfrozen_graph.remove_nodes_from(atcNodesToRemove)
#                        icd_lst = [x for x in nodes_0_c if x not in icdNodesToRemove]
#                        GPCIE = bipartite.projected_graph(unfrozen_graph, icd_lst)
#                        components = sorted(nx.connected_components(GPCIE), key=len, reverse=True)
#                        nodes_connected = sum(list(map(lambda c: len(c), components)))
#                        mean_size_components = nodes_connected / len(components)
#                        nodes_unconnected = len(nodes_0_c) - nodes_connected
#                        degrees = GPCIE.degree()
#                        sum_of_edges = sum(list(dict(degrees).values()))
#                        avg_degree = sum_of_edges / GPCIE.number_of_nodes()
#                        avg_clustering = nx.average_clustering(GPCIE)
#                        np_avg_degree[index_p][index_q] = avg_degree
#                        np_conn_nodes[index_p][index_q] = nodes_connected
#                        np_unconn_nodes[index_p][index_q] = nodes_unconnected
#                        np_conn_components[index_p][index_q] = len(components)
#                        np_mean_size[index_p][index_q] = mean_size_components
#                        np_clustering[index_p][index_q] = avg_clustering
#                        index_q += 1
#                    index_q = 0
#                    index_p += 1
#            elif type_proj == 1: # Active Ingredient projection
#                print("Calculating directed robustness active ingredients projection ...")
#                for p in p_vector:
#                    for q in p_vector:
#                        unfrozen_graph = nx.Graph(C)
#        #                atcNodesToRemove = random.sample(nodes_1_c, int(round(q * len(nodes_1_c))))
#                        icdNodesToRemove = list(df_icd.head(int(round(p * len(nodes_0_c))))['node'])
#                        atcNodesToRemove = list(df_atc.head(int(round(p * len(nodes_1_c))))['node'])
#                        unfrozen_graph.remove_nodes_from(icdNodesToRemove)
#                        unfrozen_graph.remove_nodes_from(atcNodesToRemove)
#                        atc_lst = [x for x in nodes_1_c if x not in atcNodesToRemove]
#                        GPATC = bipartite.projected_graph(unfrozen_graph, atc_lst)
#                        components = sorted(nx.connected_components(GPATC), key=len, reverse=True)
#                        nodes_connected = sum(list(map(lambda c: len(c), components)))
#                        mean_size_components = nodes_connected / len(components)
#                        nodes_unconnected = len(nodes_1_c) - nodes_connected
#                        degrees = GPATC.degree()
#                        sum_of_edges = sum(list(dict(degrees).values()))
#                        avg_degree = sum_of_edges / GPATC.number_of_nodes()
#                        avg_clustering = nx.average_clustering(GPATC)
#                        
#                        np_avg_degree[index_p][index_q] = avg_degree
#                        np_conn_nodes[index_p][index_q] = nodes_connected
#                        np_unconn_nodes[index_p][index_q] = nodes_unconnected
#                        np_conn_components[index_p][index_q] = len(components)
#                        np_mean_size[index_p][index_q] = mean_size_components
#                        np_clustering[index_p][index_q] = avg_clustering
#                        index_q += 1
#                    index_q = 0
#                    index_p += 1
                
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    cax = ax.matshow(np_clustering)
#    fig.colorbar(cax)
#    
#    ax.set_xticklabels(p_vector)
#    ax.set_yticklabels(p_vector)
#    
#    plt.show()
    