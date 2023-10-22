#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
from msilib import sequence
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
from statistics import stdev
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")
from collections import defaultdict

__author__ = "Desvilles Aurélien"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Desvilles Aurélien"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Desvilles Aurélien"
__email__ = "desvillesaurelien@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, "r") as handle:
        seq_id = None
        sequence = None
        quality = None

        for line in handle:
            line = line.strip()

            if not seq_id:
                seq_id = line
            elif not sequence:
                sequence = line
            elif not quality:
                quality = line
            else:
                # Ensure that all fields are non-empty
                if seq_id and sequence and quality:
                    yield (sequence)

                # Reset the variables
                seq_id = None
                sequence = None
                quality = None


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(0, len(read) - kmer_size + 1):
        yield read[i:i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_counts = {}

    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, kmer_size):
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

    return kmer_counts


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    G = nx.DiGraph()
    
    for kmer, weight in kmer_dict.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        if not G.has_node(prefix):
            G.add_node(prefix)
        if not G.has_node(suffix):
            G.add_node(suffix)
        
        if G.has_edge(prefix, suffix):
            G[prefix][suffix]['weight'] += weight
        else:
            G.add_edge(prefix, suffix, weight=weight)

    return G


def remove_paths(graph, path_list, delete_entry_node = True, delete_sink_node = True):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    
    modified_graph = graph.copy()

    for path in path_list:
        if delete_entry_node:
            modified_graph.remove_node_from(path[0])
        if delete_sink_node:
            modified_graph.remove_node_from(path[-1])

        for i in range(len(path) - 1):
            source_node = path[i]
            target_node = path[i + 1]
            modified_graph.remove_edge(source_node, target_node)

    return modified_graph


def select_best_path(graph, path_list, path_lengths, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    paths_to_remove = []

    stdev_weights = stdev(weight_avg_list)
    stdev_lengths = stdev(path_lengths)

    for i in range(len(path_list)):
        path = path_list[i]
        if stdev_weights > 0:
            
            if weight_avg_list[i] < max(weight_avg_list):
                paths_to_remove.append(path)
        elif stdev_lengths > 0:
            
            if path_lengths[i] < max(path_lengths):
                paths_to_remove.append(path)
    
    if not paths_to_remove:
        
        random_index = random.randint(0, len(path_list) - 1)
        paths_to_remove.pop(random_index)
    
    for path in paths_to_remove:
        for i in range(len(path) - 1):
            source_node = path[i]
            target_node = path[i + 1]
            graph.remove_edge(source_node, target_node)

        if delete_entry_node:
            graph.remove_node(path[0])

        if delete_sink_node:
            graph.remove_node(path[-1])

    return graph


def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    
    paths_between_nodes = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))

    if not paths_between_nodes:
        return graph  

    
    path_lengths = [len(path) for path in paths_between_nodes]
    path_weights = [sum(graph[source][target]['weight'] for source, target in zip(path, path[1:])) for path in paths_between_nodes]

    
    cleaned_graph = select_best_path(graph, paths_between_nodes, path_lengths, path_weights, delete_entry_node=True, delete_sink_node=True)

    return cleaned_graph


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False
    
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for i in range(len(predecessors)):
                for j in range(i+1, len(predecessors)):
                    node_ancestor = nx.lowest_common_ancestor(graph, predecessors[i], predecessors[j])
                    if node_ancestor is not None:
                        bubble = True
                        break
                if bubble:
                    break
            if bubble:
                break

    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, node_ancestor, node))
    
    return graph


def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    cleaned_graph = graph.copy()

    for entry_node in starting_nodes:
        predecessors = list(cleaned_graph.predecessors(entry_node))
        if len(predecessors) > 1:
            
            for predecessor in predecessors:
                if predecessor not in starting_nodes:
                    cleaned_graph.remove_edge(predecessor, entry_node)
        elif len(predecessors) == 1:
            
            cleaned_graph.remove_node(entry_node)

    return cleaned_graph


def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    cleaned_graph = graph.copy()

    for out_node in ending_nodes:
        successors = list(cleaned_graph.successors(out_node))
        if len(successors) > 1:
            
            for successor in successors:
                if successor not in ending_nodes:
                    cleaned_graph.remove_edge(out_node, successor)
        elif len(successors) == 1:
            
            cleaned_graph.remove_node(out_node)

    return cleaned_graph

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    starting_nodes = [node for node in graph.nodes() if len(list(graph.predecessors(node))) == 0]
    return starting_nodes

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    sink_nodes = [node for node in graph.nodes() if len(list(graph.successors(node))) == 0]
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []

    for start_node in starting_nodes:
        for sink_node in ending_nodes:
            paths = nx.all_simple_paths(graph, source=start_node, target=sink_node)
            sequence = []
            for path in paths:
                sequence = path[0]
                sequence += sum(node for node in path[1:])
                contigs.append((sequence, len(sequence)))
    return contigs

def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, "w") as file:
        for i, (contig, length) in enumerate(contigs_list, start=1):  
            header = f">Contig_{i} Length={length}\n"
            file.write(header)
            contig_wrapped = textwrap.fill(contig, width=80)
            file.write(contig_wrapped + "\n")


def draw_graph(graph, graphimg_file): 
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    
    
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # 1. Lecture du fichier et construction du graphe
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)

    # 2. Résolution des bulles
    graph = simplify_bubbles(graph)

    # 3. Résolution des pointes d'entrée et de sortie
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)

    # 4. Écriture du/des contigs
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs, args.output_file)

    # Sauvegarder le graphe en image si un nom de fichier est fourni
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)

if __name__ == '__main__': # pragma: no cover
    main()
