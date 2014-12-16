'''

A class for finding a vertex cut. A vertex cut is a a set of nodes, ideally the smallest number, that when removed yields a disconnected graph.
Apparently this is np-hard, so we settle for the following heuristic: iteratively remove the node that has the most edges until the graph is disconnected.


Examples:

    #>>> sparse_input_sequence = [("a","b",3),("b","a",3),("a","a",3),("b","b",3),("c","c",3)]
    #>>> VertexCut().work(sparse_input_sequence)
    #['a']

    >>> matrix = np.array([[3, 3, 1],[3, 3, 1],[1, 1, 3]])
    >>> VertexCut().work(matrix,2)
    [0]

'''

import numpy as np
import scipy as sp
from collections import defaultdict
import logging
import itertools

class VertexCut(object):


    def work(self, matrix, minvalue):
        assert(len(matrix.shape) == 2 and matrix.shape[0] == matrix.shape[1])

        graph = self._load_graph_from_matrix(matrix, minvalue)

        node_list = []
        while self._piece_count(graph) < len(graph):
            logging.info("len(graph)={0}".format(len(graph)))
            aMostConnectedNode = self._find_a_most_connected_node(graph)
            logging.info("aMostConnectedNode={0}".format(aMostConnectedNode))
            self._remove_node(graph, aMostConnectedNode)
            logging.info("removed")
            node_list.append(aMostConnectedNode)
        return node_list

    def _load_graph_from_matrix(self, matrix, minvalue):
        where = np.where(matrix >= minvalue)
        sparse_input_sequence = itertools.izip(where[0],where[1])
        graph = self._load_graph(sparse_input_sequence,len(where[0]))
        return graph

    def _load_graph(self, sparse_input_sequence,len_input_sequence):
        graph = defaultdict(list)
        i = -1
        for node1, node2 in sparse_input_sequence:
            i += 1
            #don't need graph[node1] because singleton's don't need to be included in the graph
            if node1 != node2 :
                if i % 100000 == 0:
                    logging.info("added {0} of {1} ({2}%)".format(i, len_input_sequence, float(i)/len_input_sequence*100.0))
                graph[node1].append(node2)
        #!!cmk self._check_that_symmetric(graph)
        return graph

    def _check_that_symmetric(self, graph):
        for node1, list in graph.iteritems():
            for node2 in list:
                if not node1 in graph[node2]:
                    raise Exception("expect symmetric graph {0}, {1}".format(node1, node2))

    def _remove_node(self, graph, node1):
        node2List = graph[node1]
        del graph[node1]
        for node2 in node2List:
            graph[node2].remove(node1)

    def _find_a_most_connected_node(self, graph):
        best_node, best_list = max(graph.iteritems(), key=lambda pair : len(pair[1])) # find the node connected to the most other nodes
        logging.info("Removing a node with {0} connections".format(len(best_list)))
        return best_node

    def _piece_count(self, graph):
        unassigned = set(graph.iterkeys())
        pieceCount = 0
        while len(unassigned) > 0:
            seed = unassigned.pop()
            pieceCount += 1

            workList = [seed]
            nextWorkList = []

            while len(workList) > 0:
                for node1 in workList:
                    for node2 in graph[node1]:
                        if node2 in unassigned:
                            unassigned.remove(node2)
                            nextWorkList.append(node2)
                workList = nextWorkList
                nextWorkList = []

        logging.info("piece_count={0}".format(pieceCount))
        return pieceCount

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import doctest
    doctest.testmod()

