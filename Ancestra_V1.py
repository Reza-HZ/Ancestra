import os
import sys
import argparse
import logging
import json
from datetime import datetime
from pathlib import Path
import csv
import re
import numpy as np
import random
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio import BiopythonWarning
import warnings
from ete3 import Tree
from matplotlib import cm
import mplcursors
from functools import lru_cache
warnings.filterwarnings("ignore", message="Attempting to set identical low and high xlims*")
warnings.simplefilter('ignore', BiopythonWarning)


IUPAC_CODES = {
    'A': {'a'}, 'C': {'c'}, 'G': {'g'}, 'T': {'t'},
    'R': {'a', 'g'}, 'Y': {'c', 't'}, 'W': {'a', 't'},
    'S': {'g', 'c'}, 'K': {'g', 't'}, 'M': {'a', 'c'},
    'B': {'c', 'g', 't'}, 'D': {'a', 'g', 't'}, 'H': {'a', 'c', 't'},
    'V': {'a', 'c', 'g'}, 'N': {'a', 'c', 'g', 't'}
}


class Pathogen:
    def __init__(self, name, a_epitopes):
        self.name = name
        self.a_epitopes = a_epitopes
        self.a_aligner = self._init_a_aligner()
        self.a_align_scores = [self.a_aligner.score(ep, ep) for ep in self.a_epitopes]

    def _init_a_aligner(self):
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
        aligner.mode = 'local'
        return aligner
    
    @lru_cache(maxsize=None)
    def _cached_score(self, s1, s2):
        return self.a_aligner.score(s1, s2)
    
    def calculate_affinity_a(self, bcr_asequence):
        max_affinity = 0
        for epitope in self.a_epitopes:
            alignment_score = self._local_a_alignment_score_normalized(bcr_asequence, epitope)
            max_affinity = max(max_affinity, alignment_score)
        return max_affinity

    
    def _local_a_alignment_score_normalized(self, a_seq1, a_seq2):
        raw_score = self._cached_score(a_seq1, a_seq2)
        if raw_score == 0:
            return 0.0
        max_self_score = self.a_align_scores[self.a_epitopes.index(a_seq2)]
        if max_self_score == 0:
            max_self_score = 1.0
        normalized_score = raw_score / max_self_score
        return max(0.0, min(1.0, normalized_score))

@lru_cache(maxsize=None)
def translate_nucleotide_to_protein_min_stops_cached(nuc_seq):
    return translate_nucleotide_to_protein_min_stops(nuc_seq)

class BCRSimulator:
    def __init__(self, v_path, d_path, j_path):
        self.vdj_segments = self._initialize_vdj_segments(v_path, d_path, j_path)
        self.mutation_patterns = self._initialize_mutation_patterns()
        self.bcr_counter = 1

    def _initialize_vdj_segments(self, v_path, d_path, j_path):
        segments = {
            "V": read_fasta(v_path),
            "D": read_fasta(d_path),
            "J": read_fasta(j_path)
        }
        return segments

    def _initialize_mutation_patterns(self):
        hotspots = {
            "WRC": 4,
            "GYW": 8,
            "DGYW": 10,
            "WRCH": 10
        }
        return hotspots

    def _generate_bcr_id(self):
        bcr_id = f"seq{self.bcr_counter}"
        self.bcr_counter += 1
        return bcr_id
    
    def _generate_initial_bcr(self):
        while True:
            v_segment = random.choice(self.vdj_segments["V"])
            d_segment = random.choice(self.vdj_segments["D"])
            j_segment = random.choice(self.vdj_segments["J"])
            flag = False
            for _ in range(10):
                n1 = ''.join(random.choices("acgt", k=random.randint(0, 10)))
                n2 = ''.join(random.choices("acgt", k=random.randint(0, 10)))
                trimmed_v = segment_trim(v_segment, seg_type='V', max_trim=12)
                trimmed_d = segment_trim(d_segment, seg_type='D', max_trim=12)
                trimmed_j = segment_trim(j_segment, seg_type='J', max_trim=12)
                bcr_sequence = trimmed_v + n1 + trimmed_d + n2 + trimmed_j
                bcr_a_sequence, min_stop, best_frame = translate_nucleotide_to_protein_min_stops(bcr_sequence)
                if min_stop == 0:
                    v_junc_end = (len(trimmed_v + n1) + 2)//3
                    v_junc_start = v_junc_end - 10
                    j_junc = len(n2 + trimmed_j)//3
                    C_pos = bcr_a_sequence[v_junc_start:v_junc_end].rfind("C")
                    W_pos = bcr_a_sequence[-j_junc:].find('W')
                    F_pos = bcr_a_sequence[-j_junc:].find('F')
                    WF_pos = F_pos if F_pos >= 0 else W_pos
                    if C_pos >= 0 and WF_pos >=0:
                        cdr3_start = C_pos + v_junc_start
                        cdr3_end = len(bcr_a_sequence) - (j_junc - WF_pos)
                        flag = True
                        break
            if flag:
                break
        return {
            "sequence": bcr_sequence,
            "a_sequence": bcr_a_sequence,
            "frame": best_frame,
            "v_gene": self.vdj_segments["V"].index(v_segment),
            "d_gene": self.vdj_segments["D"].index(d_segment),
            "j_gene": self.vdj_segments["J"].index(j_segment),
            "trimmed_v": trimmed_v,
            "trimmed_d": trimmed_d,
            "trimmed_j": trimmed_j,
            "v_len": len(v_segment),
            "d_len": len(d_segment),
            "j_len": len(j_segment),
            "n1_len": len(n1),
            "n2_len": len(n2),
            "cdr3_start": cdr3_start,
            "cdr3_end": cdr3_end,
            "cdr3_len": cdr3_end - cdr3_start + 1,
            "generation": 0,
            "mutations": 0,
            "parent": None,
            "affinity": 0.0,
            "id": "naive"
        }
    
    def _mutate_sequence(self, sequence, mutation_rate, p_trans):
        nucleotides = list(sequence)
        mutations = 0

        HOTSPOT_TARGETS = {
        "WRC": 2,
        "GYW": 0,
        "DGYW": 1,
        "WRCH": 2
        }

        coldspots = {
            "SYC": 0.6,
            "GRS": 0.5
        }
        seq_len = len(nucleotides)
        for i in range(seq_len):
            raw_score = 0.0
            target_idx = i
            for motif, weight in self.mutation_patterns.items():
                sub_seq = sequence[i:i+len(motif)]
                if len(sub_seq) == len(motif) and self._match_hotspot(sub_seq, motif):
                    target_idx = i + HOTSPOT_TARGETS[motif]
                    raw_score += weight

            mutation_prob = raw_score * mutation_rate if raw_score > 0 else mutation_rate

            for cold, penalty in coldspots.items():
                sub_seq = sequence[i:i+len(cold)]
                if len(sub_seq) == len(cold) and self._match_hotspot(sub_seq, cold):
                    mutation_prob *= penalty
            
            for cold_motif, penalty in coldspots.items():
                motif_len = len(cold_motif)
                for offset in range(-motif_len + 1, 1):
                    start = target_idx + offset
                    end = start + motif_len
                    if start >= 0 and end <= len(sequence):
                        window = sequence[start:end]
                        if self._match_hotspot(window, cold_motif):
                            mutation_prob *= penalty
                            break

            if random.random() < mutation_prob and target_idx < seq_len:
                original = nucleotides[target_idx]
                if original in ['a','c','g','t']:
                    mutated = biased_mutation(original, p_trans)
                    nucleotides[target_idx] = mutated
                    mutations += 1
        mutated_sequence = ''.join(nucleotides)
        if mutated_sequence == sequence:
            mutations = 0
        return mutated_sequence, mutations

    
    def _match_hotspot(self, sequence, hotspot_pattern):
        if len(sequence) != len(hotspot_pattern):
            return False

        for base, code in zip(sequence, hotspot_pattern):
            valid_bases = IUPAC_CODES.get(code, {code})
            if base not in valid_bases:
                return False
        return True

    def generate_repertoire_dynamically(self, pathogen, affinity_max_threshold, affinity_min_threshold, 
                                        max_generations, mutation_rate, p_sub, p_stop, p_min, p_trans):
        root = self._generate_initial_bcr()
        root["affinity"] = pathogen.calculate_affinity_a(root["a_sequence"][root['cdr3_start']:root['cdr3_end']+1])
        root["abundance"] = 1
        all_bcrs = [root]
        current_generation = [root]
        generation = 1
        thresholds = np.linspace(affinity_min_threshold, affinity_max_threshold, max_generations+1)
        while generation <= max_generations:
            next_generation = []
            for parent in current_generation:
                for _ in range(2):
                    mutated_seq, num_mutations = self._mutate_sequence(parent["sequence"], mutation_rate, p_trans)
                    if num_mutations > 0:
                        mutated_a_seq, min_stop, frame = translate_nucleotide_to_protein_min_stops_cached(mutated_seq)
                        if mutated_a_seq == parent["a_sequence"]:
                            cdr3_start = parent['cdr3_start']
                            cdr3_end = parent['cdr3_end']
                            affinity = parent['affinity']
                        else:
                            cdr3_start = mutated_a_seq.find('C', parent['cdr3_start']-2, parent['cdr3_start']+2)
                            cdr3_end = min((mutated_a_seq.find(ch, parent['cdr3_end']-3, parent['cdr3_end']+3) for ch in ['F', 'W'] if mutated_a_seq.find(ch, parent['cdr3_end']-3, parent['cdr3_end']+3) != -1), default=-1)
                            if cdr3_start == -1 or cdr3_end == -1:
                                continue
                            affinity = pathogen.calculate_affinity_a(mutated_a_seq[cdr3_start:cdr3_end+1])
                    else:
                        mutated_a_seq = parent['a_sequence']
                        min_stop = mutated_a_seq.count('*')
                        cdr3_start = parent['cdr3_start']
                        cdr3_end = parent['cdr3_end']
                        affinity = parent['affinity']
                        frame = parent['frame']
                    temp_aff_thresh = thresholds[generation]
                    abundance = 1
                    
                    if affinity >= temp_aff_thresh and min_stop == 0:
                        selection_p = 1
                    elif affinity < temp_aff_thresh and min_stop == 0:
                        selection_p = p_sub
                    elif affinity >= temp_aff_thresh and not min_stop == 0:
                        selection_p = p_stop / min_stop
                    else:
                        selection_p = p_min / min_stop
                    
                    if random.random() > selection_p:
                        continue

                    child = {
                        "id": self._generate_bcr_id(),
                        "sequence": mutated_seq,
                        "a_sequence":mutated_a_seq,
                        "frame": frame,
                        "generation": generation,
                        "parent": parent["id"],
                        "affinity": affinity,
                        "abundance": abundance,
                        "mutations": num_mutations,
                        "cdr3_start": cdr3_start,
                        "cdr3_end": cdr3_end,
                        "cdr3_len": cdr3_end - cdr3_start + 1,
                    }
                    all_bcrs.append(child)
                    next_generation.append(child)

            if not next_generation:
                break

            current_generation = next_generation
            generation += 1
        
        all_bcrs = clean_sequences(all_bcrs)
        merged_all_bcrs = merge_sequences(all_bcrs)
        return all_bcrs, merged_all_bcrs
    
    def generate_repertoire_dynamically_fullBCR(self, pathogen, affinity_max_threshold, affinity_min_threshold, 
                                        max_generations, mutation_rate, p_sub, p_stop, p_min, p_trans):
        root = self._generate_initial_bcr()
        root["affinity"] = pathogen.calculate_affinity_a(root["a_sequence"])
        root["abundance"] = 1
        all_bcrs = [root]
        current_generation = [root]
        generation = 1
        thresholds = np.linspace(affinity_min_threshold, affinity_max_threshold, max_generations+1)
        while generation <= max_generations:
            next_generation = []
            for parent in current_generation:
                for _ in range(2):
                    mutated_seq, num_mutations = self._mutate_sequence(parent["sequence"], mutation_rate, p_trans)
                    if num_mutations > 0:
                        mutated_a_seq, min_stop, frame = translate_nucleotide_to_protein_min_stops_cached(mutated_seq)
                        if mutated_a_seq == parent["a_sequence"]:
                            cdr3_start = parent['cdr3_start']
                            cdr3_end = parent['cdr3_end']
                            affinity = parent['affinity']
                        else:
                            cdr3_start = mutated_a_seq.find('C', parent['cdr3_start']-2, parent['cdr3_start']+2)
                            cdr3_end = min((mutated_a_seq.find(ch, parent['cdr3_end']-3, parent['cdr3_end']+3) for ch in ['F', 'W'] if mutated_a_seq.find(ch, parent['cdr3_end']-3, parent['cdr3_end']+3) != -1), default=-1)
                            if cdr3_start == -1 or cdr3_end == -1:
                                continue
                            affinity = pathogen.calculate_affinity_a(mutated_a_seq)
                    else:
                        mutated_a_seq = parent['a_sequence']
                        min_stop = mutated_a_seq.count('*')
                        cdr3_start = parent['cdr3_start']
                        cdr3_end = parent['cdr3_end']
                        affinity = parent['affinity']
                        frame = parent['frame']
                    temp_aff_thresh = thresholds[generation]
                    abundance = 1
                    
                    if affinity >= temp_aff_thresh and min_stop == 0:
                        selection_p = 1
                    elif affinity < temp_aff_thresh and min_stop == 0:
                        selection_p = p_sub
                    elif affinity >= temp_aff_thresh and not min_stop == 0:
                        selection_p = p_stop / min_stop
                    else:
                        selection_p = p_min / min_stop
                    
                    if random.random() > selection_p:
                        continue

                    child = {
                        "id": self._generate_bcr_id(),
                        "sequence": mutated_seq,
                        "a_sequence":mutated_a_seq,
                        "frame": frame,
                        "generation": generation,
                        "parent": parent["id"],
                        "affinity": affinity,
                        "abundance": abundance,
                        "mutations": num_mutations,
                        "cdr3_start": cdr3_start,
                        "cdr3_end": cdr3_end,
                        "cdr3_len": cdr3_end - cdr3_start + 1,
                    }
                    all_bcrs.append(child)
                    next_generation.append(child)

            if not next_generation:
                break

            current_generation = next_generation
            generation += 1
        
        all_bcrs = clean_sequences(all_bcrs)
        merged_all_bcrs = merge_sequences(all_bcrs)
        return all_bcrs, merged_all_bcrs

def clean_sequences(seq_list):
    id_map = {seq.get("id"): seq for seq in seq_list if "id" in seq}

    to_remove = []

    for seq in seq_list:
        if seq['id']=='naive':
            continue
        if seq.get("mutations", None) == 0:
            parent_id = seq.get("parent")
            abundance = seq.get("abundance", 0)

            # Add abundance to parent if parent exists
            if parent_id in id_map:
                id_map[parent_id]["abundance"] = id_map[parent_id].get("abundance", 0) + abundance

            # Redirect children to this sequence's parent
            for child in seq_list:
                if child.get("parent") == seq.get("id"):
                    child["parent"] = parent_id

            to_remove.append(seq)

    # Remove marked sequences
    for seq in to_remove:
        seq_list.remove(seq)

    return seq_list



def at_content(seq):
    return (seq.count('A') + seq.count('T')) / len(seq)

def calc_trim_probability(seq, max_trim=15):
    trim_size = min(len(seq), max_trim)
    at_tail = at_content(seq[-max_trim:]) if len(seq) > max_trim else at_content(seq)
    base_weights = [0.6 if i == 0 else
                    0.2 if 1 <= i <= 4 else
                    0.01 for i in range(trim_size + 1)]
    adjusted_weights = [
        w * (1.0 + at_tail * (1 - i / trim_size)) for i, w in enumerate(base_weights)
    ]
    total = sum(adjusted_weights)
    return [w / total for w in adjusted_weights]

def segment_trim(seq, seg_type='V', max_trim=15):
    if seg_type == 'V':
        probs = calc_trim_probability(seq, max_trim)
        trim_len = random.choices(range(len(probs)), weights=probs, k=1)[0]
        if trim_len >= len(seq):
            return ''
        return seq[:-trim_len] if trim_len > 0 else seq
    elif seg_type == 'D':
        probs_5 = calc_trim_probability(seq[:max_trim], max_trim)
        probs_3 = calc_trim_probability(seq[-max_trim:], max_trim)
        trim_start = random.choices(range(len(probs_5)), weights=probs_5, k=1)[0]
        trim_end = random.choices(range(len(probs_3)), weights=probs_3, k=1)[0]
        return seq[trim_start: len(seq)-trim_end] if (trim_start + trim_end) < len(seq) else ''

    elif seg_type == 'J':
        probs = calc_trim_probability(seq[:max_trim], max_trim)
        trim_len = random.choices(range(len(probs)), weights=probs, k=1)[0]
        if trim_len >= len(seq):
            return ''
        return seq[trim_len:] if trim_len > 0 else seq
    else:
        raise ValueError("Segment type must be one of 'V', 'D', or 'J'.")



def calculate_stats(bcrs):
    affinities = [b['affinity'] for b in bcrs]
    generations = [b['generation'] for b in bcrs]
    abundances = [b['abundance'] for b in bcrs]
    unique_seqs = set(b['sequence'] for b in bcrs)
    cdr_3s = [b['cdr3_len'] for b in bcrs]
    return {
        "naive_affinity": bcrs[0]['affinity'],
        "mean_affinity": np.mean(affinities),
        "median_affinity": np.median(affinities),
        "max_affinity": max(affinities),
        "min_affinity": min(affinities),
        "mean_abundance": np.mean(abundances),
        "max_abundance": max(abundances),
        "min_abundance": min(abundances),
        "total_abundance": sum(abundances),
        "max_generation": max(generations),
        "num_bcrs": len(bcrs),
        "num_unique_sequences": len(unique_seqs),
        "cdr3s_len": cdr_3s
    }


def plot_newick_bcellTree(newick, node_weights, output_file, title_name = 'Hi' ,title_num = 0):
    def rectangular_layout(tree):
        positions = {}
        x_offset = 0
        y_offset = 0
        level_spacing = 50  # Vertical spacing between levels
        sibling_spacing = 100  # Horizontal spacing between siblings

        def assign_positions(node, x, y):
            nonlocal x_offset
            if node.is_leaf():
                positions[node.name] = (x_offset, y)
                x_offset += sibling_spacing
            else:
                child_positions = []
                for child in node.children:
                    assign_positions(child, x, y - level_spacing)
                    child_positions.append(positions[child.name][0])
                positions[node.name] = (sum(child_positions) / len(child_positions), y)

        assign_positions(tree, x_offset, y_offset)
        return positions
    node_info = {}
    for key, val in node_weights.items():
        node_info[key] = f'Node_name: {key}\n\nAbundancy: {int(val)}'
    t = Tree(newick, format=1)
    tree_nodes = {node.name for node in t.traverse() if node.name}
    node_info.update({key: 'Node_name: ''\n\nAbundancy: 1' for key in tree_nodes if key not in node_info})
    node_weights.update({key: 1 for key in tree_nodes if key not in node_weights})
    max_weight = max(node_weights.values())
    node_sizes = {node: (weight / max_weight) * 1000 for node, weight in node_weights.items()}
    weights = np.array(list(node_weights.values()))
    if weights.max() == weights.min():
        normalized_weights = np.ones_like(weights) 
    else:
        normalized_weights = (weights - weights.min()) / (weights.max() - weights.min())  # Scale between 0 and 1
    colors = cm.viridis(normalized_weights)
    
    # Create a mapping of nodes to colors
    node_colors = {node: colors[idx] for idx, node in enumerate(node_weights.keys())}
    pos = rectangular_layout(t)
    
    # Draw the tree with rectangular edges
    fig = plt.figure(figsize=(8, 5))
    fig.canvas.manager.set_window_title(f'{title_name} T{title_num+1}')

    x_coords = [x for x, y in pos.values()]
    y_coords = [y for x, y in pos.values()]
    x_range = max(x_coords) - min(x_coords)
    y_range = max(y_coords) - min(y_coords)
    x_margin = x_range * 0.3 if x_range > 0 else 1.0
    y_margin = y_range * 0.3 if y_range > 0 else 1.0
    plt.gca().set_xlim(min(x_coords) - x_margin, max(x_coords) + x_margin)
    plt.gca().set_ylim(min(y_coords) - y_margin, max(y_coords) + y_margin)
    
    # Draw nodes
    scatter = plt.scatter([], [], s=[], alpha=0.9, color=[], edgecolor="black", zorder=2)
    scatter.set_offsets([pos[node] for node in tree_nodes])
    scatter.set_sizes([node_sizes[node] for node in tree_nodes])
    scatter.set_color([node_colors[node] for node in tree_nodes])

    for i, node in enumerate(tree_nodes):
        if node == "naive":
            plt.scatter(
                pos[node][0], pos[node][1],
                s=20,
                alpha=0.9,
                color="black",
                edgecolor="black",
                zorder=2, 
                marker="^" 
            )
    
    for node in t.traverse("postorder"):
        if not node.is_root():
            parent = node.up
            x_start, y_start = pos[parent.name]
            x_end, y_end = pos[node.name]
            plt.plot([x_start, x_end], [y_start, y_start], color="black", lw=1, zorder=1)
            plt.plot([x_end, x_end], [y_start, y_end], color="black", lw=1, zorder=1)


    plt.axis("off")
    
    cursor = mplcursors.cursor(scatter, hover=True)

    @cursor.connect("add")
    def on_add(sel):
        node_name = list(tree_nodes)[sel.index]
        info_text = node_info.get(node_name, "No information available")
        # width = 20  # Adjust the width as per your preference
        # centered_text = "\n".join(line.center(width) for line in info_text.splitlines())
        # sel.annotation.set_text(centered_text)
        sel.annotation.set_text(info_text)
        sel.annotation.set_multialignment('left') 
        bbox_properties = dict(
        boxstyle="round,pad=0.7",
        edgecolor="black",
        facecolor="yellow",
        linewidth=1,
        alpha=0.7
        )
        sel.annotation.set_bbox(bbox_properties)
        sel.annotation.arrowprops = None
        sel.annotation.get_bbox_patch()
        
        @sel.annotation.axes.figure.canvas.mpl_disconnect
        def remove_annotation(event):
            sel.annotation.set_visible(False)
            sel.annotation.axes.figure.canvas.draw_idle()
    plt.savefig(output_file, format='png')
    plt.close()

def export_to_fasta(bcrs, filename):
    with open(filename, 'w') as f:
        for b in bcrs:
            f.write(f">{b['id']}@{b['abundance']}\n{b['sequence']}\n")

def export_to_newick(bcrs, filename):
    from collections import defaultdict
    tree = defaultdict(list)
    seq_map = {}
    parent_map = {}
    
    for b in bcrs:
        b_id = b['id']
        parent = b['parent']
        if parent:
            tree[parent].append(b_id)
        seq_map[b_id] = b['sequence']
        parent_map[b_id] = parent

    roots = [b['id'] for b in bcrs if b['parent'] is None]
    if not roots:
        return

    def mutation_distance(seq1, seq2):
        return sum(a != b for a, b in zip(seq1, seq2))

    def recurse(node_id):
        children = tree[node_id]
        label = f"{node_id}"
        if parent_map[node_id] is None:
            label = "naive"
            dist = 0
        else:
            parent_seq = seq_map[parent_map[node_id]]
            child_seq = seq_map[node_id]
            dist = mutation_distance(parent_seq, child_seq)
        if not children:
            return f"{label}:{dist}"
        else:
            subtree = ",".join(recurse(c) for c in children)
            return f"({subtree}){label}:{dist}"

    newick = '('+recurse(roots[0])+')'+';'
    with open(filename, 'w') as f:
        f.write(newick)
    return newick

def translate_nucleotide_to_protein(nuc_seq):
    seq_obj = Seq(nuc_seq.upper().replace('U', 'T'))
    longest_protein = ""
    for frame in range(3):
        protein_seq = seq_obj[frame:].translate(to_stop=True)
        protein_str = str(protein_seq)
        if len(protein_str) > len(longest_protein):
            longest_protein = protein_str
    return longest_protein

def translate_nucleotide_to_protein_min_stops(nuc_seq):
    seq_obj = Seq(nuc_seq.upper().replace('U', 'T'))
    best_protein = None
    min_stops = None
    for frame in range(3):
        protein_seq = seq_obj[frame:].translate(to_stop=False)
        protein_str = str(protein_seq)
        total_stops = protein_seq.count('*')
        if min_stops is None or total_stops < min_stops:
            min_stops = total_stops
            best_protein = protein_str
            best_frame = frame
    return best_protein, min_stops, best_frame

def group_by_sequence(seq_list):

    def extract_seq_num(seq_id):
        match = re.search(r'(\d+)', seq_id)
        return int(match.group(1)) if match else float('inf')
    
    seq_to_ids = defaultdict(list)
    for entry in seq_list:
        seq_to_ids[entry['sequence']].append(entry['id'])

    result = {}
    for seq, ids in seq_to_ids.items():
        sorted_ids = sorted(ids, key=extract_seq_num)
        key_id = sorted_ids[0]
        result[key_id] = sorted_ids
    return result


def merge_sequences(seq_list):
    grouped = defaultdict(list)
    for entry in seq_list:
        key = (entry['sequence'], entry['parent'])
        grouped[key].append(entry)

    merged_map = {} 
    merged_list = []
    for entries in grouped.values():
        if len(entries) == 1:
            merged_list.append(entries[0])
        else:
            entries.sort(key=lambda x: int(x['id'][3:]))
            merged_id = entries[0]['id']
            total_abundance = sum(e['abundance'] for e in entries)
            merged_entry = entries[0].copy()
            merged_entry['id'] = merged_id
            merged_entry['abundance'] = total_abundance
            merged_list.append(merged_entry)
            for e in entries:
                merged_map[e['id']] = merged_id

    for entry in merged_list:
        parent = entry['parent']
        if parent in merged_map:
            entry['parent'] = merged_map[parent]
    for entry in merged_list:
        if entry['id'] not in merged_map:
            merged_map[entry['id']] = entry['id'] 

    for entry in merged_list:
        if entry['parent'] in merged_map:
            entry['parent'] = merged_map[entry['parent']]

    return merged_list


def biased_mutation(base, p_trans):
    transitions = {'a': 'g', 'g': 'a', 'c': 't', 't': 'c'}
    transversions = {
        'a': ['c', 't'],
        'g': ['c', 't'],
        'c': ['a', 'g'],
        't': ['a', 'g']
    }
    base = base.lower()
    if random.random() < p_trans:
        return transitions.get(base)
    else:
        return random.choice(transversions.get(base))

def read_fasta(filepath):
            return [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]






"""
BCR Simulator - Generate B-cell receptor repertoires

Usage:
  python bcr_simulator.py [--config CONFIG.yaml] [--help]
  
For quick start with defaults:
  python bcr_simulator.py
"""



# Configure logging early
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="BCR Simulator: Generate affinity-matured B-cell repertoires",
        epilog="Example: python bcr_simulator.py --clones 5 --max-gen 20 --min-seq 50"
    )
    
    # Input paths
    parser.add_argument("--input-dir", type=str, 
                       default=os.path.join(os.path.dirname(__file__), "IGH genes"),
                       help="Directory containing V/D/J FASTA and epitope.txt files")
    
    # Simulation parameters
    parser.add_argument("--clones", type=int, default=1, help="Number of independent clones (in a repertoire) to simulate")
    parser.add_argument("--max-gen", type=int, default=15, help="Maximum depth of lineage expansion")
    parser.add_argument("--min-seq", type=int, default=10, help="Minimum unique sequences required for acceptance")
    parser.add_argument("--t-max", type=float, default=0.8, help="Upper bound for selecting high-affinity BCRs during simulation")
    parser.add_argument("--t-min", type=float, default=0.3, help="Minimum affinity required for BCRs to survive early generations")
    parser.add_argument("--mu", type=float, default=0.001, help="Probability of baseline mutation per nucleotide per sequence")
    parser.add_argument("--p-sub", type=float, default=0.3, help="Survival probability for in-frame low-affinity sequences")
    parser.add_argument("--p-stop", type=float, default=0.005, help="Base survival probability for high-affinity sequences with stop codons")
    parser.add_argument("--p-min", type=float, default=0.001, help="Base survival probability for low-affinity sequences with stop codons")
    parser.add_argument("--p-trans", type=float, default=0.7, help="Bias in nucleotide substitution favoring transitions over transversions")
    parser.add_argument("--use-full-BCR", type=bool, default=False, help="Affinity calculations should use the full BCR sequence (T) or only the CDR3 region (F)")
    
    # Output control
    parser.add_argument("--output-dir", type=str, default=os.path.dirname(__file__),
                       help="Base directory for simulation outputs")
    parser.add_argument("--plot-tree", action="store_true", help="Generate lineage tree visualizations")
    parser.add_argument("--verbose", action="store_true", help="Show debug-level messages")

    return parser.parse_args()


def validate_inputs(args):
    required_files = {
        "epitope.txt": os.path.join(args.input_dir, "epitope.txt"),
        "V.fasta": os.path.join(args.input_dir, "V.fasta"),
        "D.fasta": os.path.join(args.input_dir, "D.fasta"),
        "J.fasta": os.path.join(args.input_dir, "J.fasta")
    }
    
    missing = [name for name, path in required_files.items() if not os.path.exists(path)]
    if missing:
        logger.error(f"Missing required files in {args.input_dir}: {', '.join(missing)}")
        sys.exit(1)
    
    # Validate parameter ranges
    if not (0 < args.t_min <= args.t_max <= 1.0):
        logger.error("Affinity thresholds must satisfy: 0 < t_min <= t_max <= 1.0")
        sys.exit(1)


def save_run_metadata(run_dir, args, stats, success):
    """Save run configuration and results for reproducibility."""
    metadata = {
        "timestamp": datetime.now().isoformat(),
        "parameters": vars(args),
        "stats": stats,
        "success": success
    }
    with open(os.path.join(run_dir, "run_metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)


def main():
    args = parse_arguments()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info("🚀 Starting Ancestra Simulator (guaranteeing %d successful clones)", args.clones)
    validate_inputs(args)
    
    # Setup simulator components (same as before)
    epitopes_path = os.path.join(args.input_dir, "epitope.txt")
    v_path = os.path.join(args.input_dir, "V.fasta")
    d_path = os.path.join(args.input_dir, "D.fasta")
    j_path = os.path.join(args.input_dir, "J.fasta")
    
    with open(epitopes_path, 'r') as f:
        epitope_list = [line.strip() for line in f if line.strip()]
    
    sim = BCRSimulator(v_path, d_path, j_path)
    pathogen = Pathogen("patho", a_epitopes=epitope_list)
    
    successful_runs = 0
    attempts = 0
    base_output_dir = Path(args.output_dir)
    base_output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"🎯 Target: {args.clones} successful clones (min {args.min_seq} seqs, affinity ≥ {args.t_max})\n")
    
    while successful_runs < args.clones:
        attempts += 1
        logger.info(f"\n{'='*60}")
        logger.info(f"Attempt #{attempts} | Successful: {successful_runs}/{args.clones}")
        logger.info(f"{'='*60}")
        
        try:
            # Run simulation
            if args.use_full_BCR:
                repertoire, merged_repertoire = sim.generate_repertoire_dynamically_fullBCR(
                    pathogen,
                    args.t_max,
                    args.t_min,
                    args.max_gen,
                    args.mu,
                    args.p_sub,
                    args.p_stop,
                    args.p_min,
                    args.p_trans
                )
            else:
                repertoire, merged_repertoire = sim.generate_repertoire_dynamically(
                    pathogen,
                    args.t_max,
                    args.t_min,
                    args.max_gen,
                    args.mu,
                    args.p_sub,
                    args.p_stop,
                    args.p_min,
                    args.p_trans
                )
            
            stats = calculate_stats(merged_repertoire)
            success = (stats['num_unique_sequences'] >= args.min_seq and 
                      stats['max_affinity'] >= args.t_max)
            
            # Determine output folder: success gets numbered run, failure goes to "rejected"
            if success:
                successful_runs += 1
                run_folder = base_output_dir / "successful_clones" / f"run_{successful_runs}_gen{args.max_gen}_minSeq{args.min_seq}"
                status_emoji = "✅"
                status_msg = "ACCEPTED"
            else:
                run_folder = base_output_dir / "rejected_clones" / f"rejected_attempt_{attempts}"
                status_emoji = "❌"
                status_msg = "REJECTED"
            
            run_folder.mkdir(parents=True, exist_ok=True)
            
            # with open(run_folder / "Statistics.txt", 'w') as f:
            #     for k, v in stats.items():
            #         f.write(f"{k}: {v}\n")
            
            export_to_fasta(merged_repertoire, run_folder / "repertoire.fasta")
            newick = export_to_newick(merged_repertoire, run_folder / "repertoire.nk")
            
            if args.plot_tree and success:  # Only plot trees for successful runs
                node_weights = {item['id']: item['abundance'] for item in merged_repertoire}
                plot_newick_bcellTree(
                    newick, node_weights, 
                    run_folder / "lineage_tree.png",
                    title_name=f"RUN {successful_runs}",
                    title_num=successful_runs
                )
            
            # Save sequence info
            keys_to_include = ['id', 'a_sequence', 'frame', 'cdr3_start', 'cdr3_end', 
                             'cdr3_len', 'generation', 'parent', 'mutations', 
                             'affinity', 'abundance']
            with open(run_folder / "repertoire_info.csv", 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=keys_to_include)
                writer.writeheader()
                for entry in merged_repertoire:
                    writer.writerow({k: entry[k] for k in keys_to_include if k in entry})
            
            # Save metadata
            save_run_metadata(run_folder, args, stats, success)
            
            # Log outcome
            logger.info(f"{status_emoji} {status_msg} | Unique seqs: {stats['num_unique_sequences']}/{args.min_seq} | "
                       f"Max affinity: {stats['max_affinity']:.3f}/{args.t_max}")
            logger.info(f"📁 Saved to: {run_folder.name}")
            
            # Early exit if we hit target
            if successful_runs == args.clones:
                logger.info(f"\n🎉 Target reached! {args.clones} successful clones generated in {attempts} attempts.")
                break
                
        except Exception as e:
            logger.exception(f"💥 Simulation crashed on attempt #{attempts}: {e}")
            # Continue to next attempt rather than aborting
    
    # Final summary
    logger.info("\n" + "="*60)
    if successful_runs == args.clones:
        logger.info(f"✅ SUCCESS: Generated {successful_runs}/{args.clones} successful clones")
        logger.info(f"   Total attempts: {attempts}")
        logger.info(f"   Success rate: {successful_runs/attempts*100:.1f}%")
        logger.info(f"   Outputs in: {base_output_dir}")
    else:
        logger.error(f"❌ FAILED: Only {successful_runs}/{args.clones} successful clones after {attempts} attempts")
        logger.error(f"   Possible causes:")
        logger.error(f"   • Thresholds too strict (min_seq={args.min_seq}, t_max={args.t_max})")
        logger.error(f"   • max_gen={args.max_gen} too low for affinity maturation")
        logger.error(f"   • Check rejected_attempt_* folders for diagnostics")
        sys.exit(1)
    logger.info("="*60)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("\n⚠️ Simulation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.exception(f"💥 Unexpected error: {e}")
        sys.exit(1)
