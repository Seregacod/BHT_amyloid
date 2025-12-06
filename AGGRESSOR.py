#!/usr/bin/env python3
"""
In Silico Mutagenesis Script with Rule-Based Mutations
Performs point mutations and insertions at specified positions in protein sequences
with rule-based mutagenesis in specified regions

Rules apply only when amino acids from the rule are clustered together
Multiple rules can apply simultaneously to the same motif

Usage: python mutagenesis.py <input_file> [options]
"""

import argparse
import re
import sys
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Set, Any, Union, Optional

# Defining constants
VALID_AAS = set('ACDEFGHIKLMNPQRSTVWY')
FASTA_LINE_LENGTH = 60
MAX_GAP_FOR_MERGING = 2

# Default mutation list - customizable here
DEFAULT_MUTATIONS = ['P', 'G', 'D', 'K']

# Gatekeeping amino acids (only applied to edge positions)
GATEKEEPING_AAS = ['Y']

# Data classes for cluster structures
@dataclass
class Cluster:
    """Represents a cluster of amino acid positions matching a rule"""
    positions: List[int]
    residues: List[str]
    rule_name: str
    aggregation_score: int
    size: int = field(init=False)
    span: int = field(init=False)
    
    def __post_init__(self):
        """Calculate size and span after initialization"""
        self.size = len(self.positions)
        if self.positions:
            self.span = max(self.positions) - min(self.positions) + 1
        else:
            self.span = 0
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for backward compatibility"""
        return {
            'positions': self.positions,
            'residues': self.residues,
            'size': self.size,
            'span': self.span,
            'rule_name': self.rule_name,
            'aggregation_score': self.aggregation_score
        }

@dataclass
class HydrophobicAromaticCluster(Cluster):
    """Special cluster for hydrophobic-aromatic interactions"""
    pair_count: int = 0
    condition: str = ""
    hydrophobic_count: int = 0
    aromatic_count: int = 0
    pairs: List[Dict] = field(default_factory=list)
    nearby_hydrophobics: List[Dict] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary with special fields"""
        result = super().to_dict()
        result.update({
            'pair_count': self.pair_count,
            'condition': self.condition,
            'hydrophobic_count': self.hydrophobic_count,
            'aromatic_count': self.aromatic_count,
            'pairs': self.pairs,
            'nearby_hydrophobics': self.nearby_hydrophobics
        })
        return result

@dataclass
class MultiRuleCluster:
    """Represents a cluster matching multiple rules"""
    positions: List[int]
    residues: List[str]
    rules: List[str]
    combined_aggregation_score: int
    size: int = field(init=False)
    span: int = field(init=False)
    is_multi_rule: bool = True
    
    def __post_init__(self):
        """Calculate size and span after initialization"""
        self.size = len(self.positions)
        if self.positions:
            self.span = max(self.positions) - min(self.positions) + 1
        else:
            self.span = 0
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for backward compatibility"""
        return {
            'positions': self.positions,
            'residues': self.residues,
            'size': self.size,
            'span': self.span,
            'rules': self.rules,
            'combined_aggregation_score': self.combined_aggregation_score,
            'is_multi_rule': self.is_multi_rule
        }

# Define rules for mutagenesis with clustering requirement
RULES = {
    'hydrophobic_aliphatic': {
        'description': 'V, I, L, A, M present and clustered (≥3 residues within 4 positions)',
        'residues': {'V', 'I', 'L', 'A', 'M'},
        'min_cluster_size': 3,
        'max_gap': 4,
        'mutations': DEFAULT_MUTATIONS,
        'priority': 1,
        'aggregation_score': 3
    },
    'aromatic': {
        'description': 'F, Y, W present and clustered (any combination within 3 positions)',
        'residues': {'F', 'Y', 'W'},
        'min_cluster_size': 2,
        'max_gap': 3,
        'mutations': DEFAULT_MUTATIONS,
        'priority': 2,
        'aggregation_score': 2
    },
    'amide': {
        'description': 'Q, N present and clustered (any combination within 3 positions)',
        'residues': {'Q', 'N'},
        'min_cluster_size': 2,
        'max_gap': 3,
        'mutations': DEFAULT_MUTATIONS,
        'priority': 3,
        'aggregation_score': 1
    },
    'hydrophobic_and_aromatic': {
        'description': 'Hydrophobic (V,I,L,A,M) adjacent to aromatic (F,Y,W) (≥2 such pairs OR 1 pair + hydrophobic)',
        'residues': {'V', 'I', 'L', 'A', 'M', 'F', 'Y', 'W'},
        'min_cluster_size': 2,
        'max_gap': 1,
        'mutations': DEFAULT_MUTATIONS,
        'priority': 4,
        'aggregation_score': 2,
        'special_rule': True
    }
}

def print_usage_example():
    """Print detailed usage example"""
    example = """
USAGE EXAMPLES:
=========================================================

1. Rule-based mutagenesis in specific regions:
   python mutagenesis.py protein.fasta --regions 10:20 30:40 50:60

2. Rule-based with custom mutations:
   python mutagenesis.py protein.fasta --regions 5:15 -m A D E

3. Rule-based with specific rules only:
   python mutagenesis.py protein.fasta --regions 10:30 --rules hydrophobic_aliphatic aromatic

4. Combined approach (rules + specific positions):
   python mutagenesis.py protein.fasta --regions 10:20 --positions 15 25 --mutations P G

5. With insertions and rule-based:
   python mutagenesis.py protein.fasta --regions 5:15 --insert-positions 10 --insert-aas K

6. With gatekeeping amino acids (only for edge positions):
   python mutagenesis.py protein.fasta --regions 10:20 --gatekeeping Y K

7. Detailed verbose output:
   python mutagenesis.py protein.fasta --regions 10:20 -v

AVAILABLE RULES:
=========================================================
• hydrophobic_aliphatic    : Triggers if ≥3 V, I, L, A, M residues within 4 positions of each other
• aromatic                 : Triggers if ≥2 F, Y, W residues within 3 positions of each other
• amide                    : Triggers if ≥2 Q, N residues within 3 positions of each other
• hydrophobic_and_aromatic : Triggers if ≥2 hydrophobic-aromatic adjacent pairs 
                             OR 1 pair + at least 1 hydrophobic residue within 3 positions

AGGREGATION SCORE RANKING (highest to lowest):
1. hydrophobic_aliphatic: 3
2. hydrophobic_and_aromatic: 2
3. aromatic: 2
4. amide: 1

GATEKEEPING AMINO ACIDS:
=========================================================
• Gatekeeping amino acids (default: Y) are only applied to positions at the edge of motifs
  or directly adjacent to them (within 1 position of motif boundary)
• Regular mutations are applied to all positions in motifs
• Use --gatekeeping option to specify custom gatekeeping amino acids

REQUIRED PARAMETERS:
=========================================================
• input_file    : Input FASTA file containing protein sequence (the multifasta format is not supported!)
• Either --positions OR --regions must be specified
"""
    print(example)

def read_fasta(filepath: str) -> Tuple[str, str]:
    """
    Read a single sequence FASTA file
    Returns: (header, sequence)
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        header = ""
        sequence = ""
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if header:  # If we already have a header, we're done (single sequence)
                    break
                header = line
            elif line and not header:
                raise ValueError("FASTA file must start with header line (>)")
            elif header:
                sequence += line.upper()
        
        if not header or not sequence:
            raise ValueError("Invalid FASTA file or empty sequence")
        
        return header, sequence
    
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")
    except Exception as e:
        raise ValueError(f"Error reading FASTA file: {e}")

def validate_amino_acids(input_data: Union[str, List[str]], name: str = "amino acids", strict: bool = False) -> bool:
    """
    Unified function to validate amino acid codes.
    
    Args:
        input_data: Can be a string (sequence) or list of strings (amino acid codes)
        name: Description for error messages (used only in strict mode)
        strict: If True, raises ValueError on invalid input. If False, returns bool.
    
    Returns:
        bool: True if all amino acids are valid (only when strict=False)
    
    Raises:
        ValueError: If strict=True and invalid amino acids are found
    """
    if isinstance(input_data, str):
        # Validate sequence string: check if all characters are valid AAs
        invalid_chars = [char for char in set(input_data.upper()) if char not in VALID_AAS]
        if invalid_chars:
            if strict:
                raise ValueError(f"Invalid {name}: contains invalid characters {invalid_chars}. "
                               f"Valid amino acids: {', '.join(sorted(VALID_AAS))}")
            return False
        return True
    
    elif isinstance(input_data, list):
        # Validate list of amino acid codes: check if all items are single valid AAs
        invalid = [aa for aa in input_data if len(aa) != 1 or aa.upper() not in VALID_AAS]
        if invalid:
            if strict:
                raise ValueError(f"Invalid {name}: {invalid}. Valid amino acids: {', '.join(sorted(VALID_AAS))}")
            return False
        return True
    
    else:
        raise TypeError(f"input_data must be str or list, got {type(input_data).__name__}")

def validate_sequence(sequence: str) -> bool:
    """Validate that the sequence contains only valid amino acid codes (non-strict)"""
    return validate_amino_acids(sequence, "sequence", strict=False)

def parse_region(region_str: str, seq_length: int) -> Tuple[int, int]:
    """Parse region string in format start:stop (1-indexed)"""
    try:
        if ':' not in region_str:
            raise ValueError("Region must be in format start:stop")
        
        start_str, stop_str = region_str.split(':')
        start = int(start_str.strip())
        stop = int(stop_str.strip())
        
        # Validate bounds
        if start < 1 or stop > seq_length:
            raise ValueError(f"Region {start}:{stop} out of bounds (1-{seq_length})")
        if start > stop:
            raise ValueError(f"Start position {start} cannot be greater than stop position {stop}")
        
        return start, stop
    except ValueError as e:
        raise ValueError(f"Invalid region format '{region_str}': {e}")

def find_clusters(positions: List[int], max_gap: int = 3) -> List[List[int]]:
    """
    Find clusters of positions where positions are within max_gap of each other
    
    Args:
        positions: Sorted list of positions
        max_gap: Maximum allowed gap between consecutive positions in a cluster
    
    Returns:
        List of clusters, each cluster is a list of positions
    """
    if not positions:
        return []
    
    clusters = []
    current_cluster = [positions[0]]
    
    for i in range(1, len(positions)):
        if positions[i] - current_cluster[-1] <= max_gap:
            current_cluster.append(positions[i])
        else:
            if len(current_cluster) >= 2:  # Only keep clusters with at least 2 positions
                clusters.append(current_cluster)
            current_cluster = [positions[i]]
    
    # Don't forget the last cluster
    if len(current_cluster) >= 2:
        clusters.append(current_cluster)
    
    return clusters

def find_hydrophobic_aromatic_interactions(sequence: str, start: int, stop: int) -> List[Dict]:
    """
    Find hydrophobic-aromatic interactions with two conditions:
    1. At least 2 hydrophobic-aromatic adjacent pairs
    2. OR 1 hydrophobic-aromatic pair + at least 1 hydrophobic residue within 3 positions
    """
    region_seq = sequence[start-1:stop]
    hydrophobic_set = {'V', 'I', 'L', 'A', 'M'}
    aromatic_set = {'F', 'Y', 'W'}
    
    pairs = []
    aromatic_positions = []
    hydrophobic_positions = []
    
    # Collect all aromatic and hydrophobic positions
    for i, aa in enumerate(region_seq):
        pos = start + i
        if aa in aromatic_set:
            aromatic_positions.append(pos)
        if aa in hydrophobic_set:
            hydrophobic_positions.append(pos)
    
    # Find adjacent hydrophobic-aromatic pairs
    for i in range(len(region_seq) - 1):
        pos1 = start + i
        pos2 = start + i + 1
        
        aa1 = region_seq[i]
        aa2 = region_seq[i + 1]
        
        # Check for hydrophobic-aromatic or aromatic-hydrophobic pairs
        if (aa1 in hydrophobic_set and aa2 in aromatic_set) or \
           (aa1 in aromatic_set and aa2 in hydrophobic_set):
            pairs.append({
                'positions': [pos1, pos2],
                'residues': [aa1, aa2],
                'type': 'hydrophobic_aromatic_pair',
                'interaction': f"{aa1}{pos1}-{aa2}{pos2}",
                'is_hydrophobic_aromatic': True
            })
    
    # If no pairs found, return empty
    if not pairs:
        return []
    
    # Now check the two conditions
    valid_clusters = []
    
    # Condition 1: At least 2 pairs (original condition)
    if len(pairs) >= 2:
        # Check if pairs are close to each other (within 3 positions)
        all_pair_positions = []
        for pair in pairs:
            all_pair_positions.extend(pair['positions'])
        all_pair_positions = sorted(set(all_pair_positions))
        
        # Find clusters of these positions
        position_clusters = find_clusters(all_pair_positions, max_gap=3)
        
        for cluster in position_clusters:
            # Find pairs within this cluster
            cluster_pairs = []
            for pair in pairs:
                if pair['positions'][0] in cluster and pair['positions'][1] in cluster:
                    cluster_pairs.append(pair)
            
            if len(cluster_pairs) >= 2:
                # Get all residues in cluster
                cluster_residues = [sequence[pos-1] for pos in cluster]
                valid_clusters.append({
                    'positions': cluster,
                    'residues': cluster_residues,
                    'pairs': cluster_pairs,
                    'size': len(cluster),
                    'pair_count': len(cluster_pairs),
                    'span': max(cluster) - min(cluster) + 1,
                    'condition': 'at_least_2_pairs',
                    'hydrophobic_count': sum(1 for res in cluster_residues if res in hydrophobic_set),
                    'aromatic_count': sum(1 for res in cluster_residues if res in aromatic_set)
                })
    
    # Condition 2: 1 pair + at least 1 hydrophobic residue within 3 positions (CORRECTED)
    # Check each pair separately
    for pair in pairs:
        pair_positions = set(pair['positions'])
        
        # Look for hydrophobic residues within 3 positions of the pair
        nearby_hydrophobics = []
        for hydrophobic_pos in hydrophobic_positions:
            if hydrophobic_pos in pair_positions:
                continue  # Skip if it's already part of the pair
                
            # Check distance to either position in the pair
            min_distance = min(abs(hydrophobic_pos - pos) for pos in pair['positions'])
            if min_distance <= 3:  # Within 3 positions
                nearby_hydrophobics.append({
                    'position': hydrophobic_pos,
                    'residue': sequence[hydrophobic_pos-1],
                    'distance': min_distance
                })
        
        # If we have at least 1 hydrophobic nearby
        if nearby_hydrophobics:
            # Create cluster with pair + nearby hydrophobics
            all_positions = sorted(list(pair_positions) + [h['position'] for h in nearby_hydrophobics])
            cluster_residues = [sequence[pos-1] for pos in all_positions]
            
            # Count totals in cluster
            hydrophobic_count = sum(1 for res in cluster_residues if res in hydrophobic_set)
            aromatic_count = sum(1 for res in cluster_residues if res in aromatic_set)
            
            valid_clusters.append({
                'positions': all_positions,
                'residues': cluster_residues,
                'pairs': [pair],
                'size': len(all_positions),
                'pair_count': 1,
                'nearby_hydrophobics': nearby_hydrophobics,
                'span': max(all_positions) - min(all_positions) + 1,
                'condition': '1_pair_plus_hydrophobic',
                'hydrophobic_count': hydrophobic_count,
                'aromatic_count': aromatic_count
            })
    
    # Remove duplicate clusters (clusters with identical positions)
    unique_clusters = []
    seen_positions = set()
    
    for cluster in valid_clusters:
        pos_tuple = tuple(cluster['positions'])
        if pos_tuple not in seen_positions:
            seen_positions.add(pos_tuple)
            unique_clusters.append(cluster)
    
    return unique_clusters

def analyze_region(sequence: str, start: int, stop: int) -> Dict:
    """Analyze a region for rule compliance with clustering requirement"""
    region_seq = sequence[start-1:stop]
    results = {
        'region': (start, stop),
        'sequence': region_seq,
        'length': len(region_seq),
        'rules': {},
        'multi_rule_clusters': [],
        'aggregation_hotspots': []
    }
    
    # Analyze standard rules first
    for rule_name, rule in RULES.items():
        if rule_name == 'hydrophobic_and_aromatic':
            continue  # Handle separately
        
        # Find positions of matching residues in this region
        matching_positions = []
        matching_residues = []
        
        for i, aa in enumerate(region_seq):
            if aa in rule['residues']:
                matching_positions.append(i + start)
                matching_residues.append(aa)
        
        # Find clusters of matching residues
        clusters = find_clusters(matching_positions, rule['max_gap'])
        
        # Check if any cluster meets the size requirement
        condition_met = False
        qualifying_clusters = []
        
        for cluster in clusters:
            if len(cluster) >= rule['min_cluster_size']:
                condition_met = True
                cluster_residues = [sequence[pos-1] for pos in cluster]
                qualifying_clusters.append({
                    'positions': cluster,
                    'residues': cluster_residues,
                    'size': len(cluster),
                    'span': max(cluster) - min(cluster) + 1,
                    'rule_name': rule_name,
                    'aggregation_score': rule['aggregation_score']
                })
        
        results['rules'][rule_name] = {
            'description': rule['description'],
            'matching_residues': matching_residues,
            'matching_positions': matching_positions,
            'total_count': len(matching_positions),
            'clusters': clusters,
            'qualifying_clusters': qualifying_clusters,
            'condition_met': condition_met,
            'min_cluster_size': rule['min_cluster_size'],
            'max_gap': rule['max_gap'],
            'priority': rule['priority'],
            'aggregation_score': rule['aggregation_score']
        }
    
    # Special handling for hydrophobic_and_aromatic rule
    hydrophobic_aromatic_clusters = find_hydrophobic_aromatic_interactions(sequence, start, stop)
    if hydrophobic_aromatic_clusters:
        results['rules']['hydrophobic_and_aromatic'] = {
            'description': RULES['hydrophobic_and_aromatic']['description'],
            'matching_residues': [res for cluster in hydrophobic_aromatic_clusters for res in cluster['residues']],
            'matching_positions': [pos for cluster in hydrophobic_aromatic_clusters for pos in cluster['positions']],
            'total_count': len([pos for cluster in hydrophobic_aromatic_clusters for pos in cluster['positions']]),
            'clusters': [cluster['positions'] for cluster in hydrophobic_aromatic_clusters],
            'qualifying_clusters': [],
            'condition_met': True,
            'min_cluster_size': RULES['hydrophobic_and_aromatic']['min_cluster_size'],
            'max_gap': RULES['hydrophobic_and_aromatic']['max_gap'],
            'priority': RULES['hydrophobic_and_aromatic']['priority'],
            'aggregation_score': RULES['hydrophobic_and_aromatic']['aggregation_score'],
            'special_clusters': hydrophobic_aromatic_clusters
        }
        
        # Add qualifying clusters for this rule
        for cluster in hydrophobic_aromatic_clusters:
            # Create a proper cluster dictionary with all required fields
            qual_cluster = {
                'positions': cluster['positions'],
                'residues': cluster['residues'],
                'size': cluster['size'],
                'span': cluster['span'],
                'rule_name': 'hydrophobic_and_aromatic',
                'aggregation_score': RULES['hydrophobic_and_aromatic']['aggregation_score'],
                'pair_count': cluster.get('pair_count', 0),
                'condition': cluster.get('condition', ''),
                'hydrophobic_count': cluster.get('hydrophobic_count', 0),
                'aromatic_count': cluster.get('aromatic_count', 0),
                'pairs': cluster.get('pairs', []),
                'nearby_hydrophobics': cluster.get('nearby_hydrophobics', [])
            }
            results['rules']['hydrophobic_and_aromatic']['qualifying_clusters'].append(qual_cluster)
    else:
        results['rules']['hydrophobic_and_aromatic'] = {
            'description': RULES['hydrophobic_and_aromatic']['description'],
            'matching_residues': [],
            'matching_positions': [],
            'total_count': 0,
            'clusters': [],
            'qualifying_clusters': [],
            'condition_met': False,
            'min_cluster_size': RULES['hydrophobic_and_aromatic']['min_cluster_size'],
            'max_gap': RULES['hydrophobic_and_aromatic']['max_gap'],
            'priority': RULES['hydrophobic_and_aromatic']['priority'],
            'aggregation_score': RULES['hydrophobic_and_aromatic']['aggregation_score'],
            'special_clusters': []
        }
    
    # Now identify clusters that match multiple rules (optimized O(n log n))
    all_clusters = []
    for rule_name, rule_data in results['rules'].items():
        for cluster in rule_data['qualifying_clusters']:
            # Use the rule_name from the cluster
            all_clusters.append((cluster['rule_name'], cluster))
    
    # Use optimized overlap detection
    if all_clusters:
        merged_multi_clusters = find_overlapping_clusters_optimized(all_clusters, sequence)
        for merged_cluster in merged_multi_clusters:
            if merged_cluster.get('is_multi_rule', False):
                results['multi_rule_clusters'].append(merged_cluster)
                results['aggregation_hotspots'].extend(merged_cluster['positions'])
    
    # Remove duplicates from aggregation hotspots
    results['aggregation_hotspots'] = sorted(set(results['aggregation_hotspots']))
    
    return results

def find_overlapping_clusters_optimized(all_clusters: List[Tuple[str, Dict]], sequence: str) -> List[Dict]:
    """
    Optimized O(n log n) overlap detection using spatial sorting.
    Returns list of multi-rule clusters.
    """
    if not all_clusters:
        return []
    
    # Create cluster intervals with metadata (start, end, rule_name, cluster_data)
    intervals = []
    for rule_name, cluster_data in all_clusters:
        positions = cluster_data.get('positions', [])
        if positions:
            # Clean and validate rule_name
            if not rule_name or not isinstance(rule_name, str):
                rule_name = 'unknown'
            
            intervals.append({
                'start': min(positions),
                'end': max(positions),
                'rule_name': rule_name,
                'cluster_data': cluster_data,
                'positions_set': set(positions)
            })
    
    if not intervals:
        return []
    
    # Sort by start position (O(n log n))
    intervals.sort(key=lambda x: x['start'])
    
    # Single pass merge (O(n))
    merged_clusters = []
    current_group = [intervals[0]]
    
    for interval in intervals[1:]:
        # Check if current interval overlaps or is close to the last in current group
        last_in_group = current_group[-1]
        gap = interval['start'] - last_in_group['end'] - 1
        
        if interval['positions_set'].intersection(last_in_group['positions_set']) or gap <= MAX_GAP_FOR_MERGING:
            # Merge into current group
            current_group.append(interval)
        else:
            # Finalize current group and start new one
            if len(current_group) > 1:
                merged_clusters.append(_merge_cluster_group(current_group, sequence))
            current_group = [interval]
    
    # Don't forget the last group
    if len(current_group) > 1:
        merged_clusters.append(_merge_cluster_group(current_group, sequence))
    
    return merged_clusters

def _merge_cluster_group(group: List[Dict], sequence: str) -> Dict:
    """Merge a group of overlapping clusters into a single cluster"""
    all_positions = set()
    all_rule_names = set()  # Use set to collect all individual rule names
    max_score = 0
    total_score = 0
    
    for item in group:
        all_positions.update(item['positions_set'])
        # Handle merged rule names - split by '+' to get individual rules
        rule_name = item.get('rule_name', '')
        if not rule_name:
            rule_name = 'unknown'
        
        # Filter out any empty strings AND only keep valid rule names
        individual_rules = []
        for r in rule_name.split('+'):
            r = r.strip()
            if r and r in RULES:
                individual_rules.append(r)
        
        all_rule_names.update(individual_rules)
        score = item['cluster_data'].get('aggregation_score', 0)
        max_score = max(max_score, score)
        total_score += score
    
    sorted_positions = sorted(all_positions)
    combined_residues = [sequence[pos-1] for pos in sorted_positions]
    
    # Get unique rules (already unique since we used a set)
    unique_rules = sorted(list(all_rule_names))
    
    # If no valid rules found, check if we have any non-standard rule names
    if not unique_rules:
        # Try to extract any rule-like names that might have been created
        for item in group:
            rule_name = item.get('rule_name', '')
            # Look for patterns like rule names in the string
            for possible_rule in RULES.keys():
                if possible_rule in rule_name:
                    unique_rules.append(possible_rule)
        
        # Remove duplicates again
        unique_rules = sorted(list(set(unique_rules)))
    
    # If still no rules, use a placeholder but ensure it's valid for downstream processing
    if not unique_rules:
        # Use the first valid rule as a fallback
        unique_rules = [list(RULES.keys())[0]]
    
    if len(unique_rules) > 1:
        combined_score = total_score + len(unique_rules) - 1
    else:
        combined_score = max_score
    
    return {
        'positions': sorted_positions,
        'residues': combined_residues,
        'size': len(sorted_positions),
        'span': max(sorted_positions) - min(sorted_positions) + 1 if sorted_positions else 0,
        'rules': unique_rules,
        'combined_aggregation_score': combined_score,
        'is_multi_rule': len(unique_rules) > 1,
        'rule_name': '+'.join(unique_rules) if len(unique_rules) > 1 else unique_rules[0],
        'aggregation_score': combined_score
    }

def resolve_overlapping_clusters(all_clusters: List[Tuple[str, Dict]], sequence: str) -> List[Tuple[str, Dict]]:
    """
    Resolve overlapping clusters by keeping the union of motifs when overlaps occur.
    Additionally, merge clusters if the distance between two adjacent motifs is ≤ 2.
    Returns a filtered list of clusters with no overlapping positions.
    """
    if not all_clusters:
        return []
    
    # Extract clusters and their positions
    clusters = []
    for rule_name, cluster_data in all_clusters:
        # Clean up rule_name before storing
        # Filter out empty strings and non-valid rule names
        valid_rules = [r for r in rule_name.split('+') if r and r in RULES]
        if not valid_rules:
            # If no valid rules, check if rule_name contains any valid rule
            for possible_rule in RULES.keys():
                if possible_rule in rule_name:
                    valid_rules.append(possible_rule)
        
        # Create cleaned rule name
        cleaned_rule_name = '+'.join(sorted(set(valid_rules))) if valid_rules else list(RULES.keys())[0]
        
        clusters.append({
            'rule_name': cleaned_rule_name,
            'positions': set(cluster_data['positions']),
            'data': cluster_data
        })
    
    # First pass: Merge clusters that overlap or are close to each other
    merged = True
    while merged:
        merged = False
        new_clusters = []
        
        while clusters:
            current = clusters.pop(0)
            merged_current = False
            
            for i, other in enumerate(clusters):
                # Check for overlap or closeness
                min_current = min(current['positions'])
                max_current = max(current['positions'])
                min_other = min(other['positions'])
                max_other = max(other['positions'])
                
                # Calculate the gap between clusters
                if min_current < min_other:
                    gap = min_other - max_current - 1
                else:
                    gap = min_current - max_other - 1
                
                # Merge if overlapping OR gap ≤ 2
                if current['positions'].intersection(other['positions']) or gap <= 2:
                    # Merge the two clusters (union of positions)
                    merged_positions = current['positions'].union(other['positions'])
                    
                    # Combine data from both clusters
                    if len(merged_positions) > len(current['positions']):
                        # Use the cluster with more residues as base
                        if len(current['positions']) >= len(other['positions']):
                            base_data = current['data'].copy()
                        else:
                            base_data = other['data'].copy()
                        
                        # Update with merged positions
                        base_data['positions'] = sorted(merged_positions)
                        base_data['residues'] = [sequence[pos-1] for pos in sorted(merged_positions)]
                        base_data['size'] = len(merged_positions)
                        base_data['span'] = max(merged_positions) - min(merged_positions) + 1
                        
                        # Get rule sets for comparison (filtering out empty strings and invalid rules)
                        current_rules_set = set([r for r in current['rule_name'].split('+') if r and r in RULES])
                        other_rules_set = set([r for r in other['rule_name'].split('+') if r and r in RULES])
                        
                        # If merging different rules, create a multi-rule cluster
                        if current_rules_set != other_rules_set:
                            # Merge the rules
                            merged_rules = sorted(current_rules_set.union(other_rules_set))
                            # Ensure we have at least one valid rule
                            if not merged_rules:
                                # Try to extract rules from original rule names
                                for rule_str in [current['rule_name'], other['rule_name']]:
                                    for possible_rule in RULES.keys():
                                        if possible_rule in rule_str and possible_rule not in merged_rules:
                                            merged_rules.append(possible_rule)
                                merged_rules = sorted(list(set(merged_rules)))
                            
                            # If still no rules, use a default
                            if not merged_rules:
                                merged_rules = [list(RULES.keys())[0]]
                                
                            base_data['rule_name'] = '+'.join(merged_rules)
                            # Sum aggregation scores
                            score1 = current['data'].get('aggregation_score', 0)
                            score2 = other['data'].get('aggregation_score', 0)
                            base_data['aggregation_score'] = score1 + score2 + 1
                        else:
                            base_data['aggregation_score'] = max(
                                current['data'].get('aggregation_score', 0),
                                other['data'].get('aggregation_score', 0)
                            )
                        
                        # Create merged cluster
                        merged_cluster = {
                            'rule_name': base_data['rule_name'],
                            'positions': merged_positions,
                            'data': base_data
                        }
                        
                        # Remove the other cluster from the list
                        clusters.pop(i)
                        
                        # Add the merged cluster to be processed again
                        clusters.append(merged_cluster)
                        merged_current = True
                        merged = True
                        break
                    
                    # If we're merging identical or very similar clusters, just remove one
                    else:
                        clusters.pop(i)
                        merged_current = True
                        merged = True
                        break
            
            if not merged_current:
                new_clusters.append(current)
        
        clusters = new_clusters
    
    # Convert back to original format
    result = []
    for cluster_info in clusters:
        result.append((cluster_info['rule_name'], cluster_info['data']))
    
    return result

def create_mutated_sequence(sequence: str, position: int, new_aa: str) -> str:
    """
    Efficiently create a mutated sequence by replacing one amino acid.
    Optimized to avoid unnecessary list creation for large sequences.
    """
    if position < 1 or position > len(sequence):
        raise ValueError(f"Position {position} out of range (1-{len(sequence)})")
    
    # For small sequences, list conversion is fine
    # For large sequences, use slicing (more memory efficient)
    if len(sequence) < 1000:
        seq_list = list(sequence)
        seq_list[position - 1] = new_aa
        return ''.join(seq_list)
    else:
        # For large sequences, use string slicing (avoids full list creation)
        return sequence[:position-1] + new_aa + sequence[position:]

def get_mutations_for_position(pos: int, cluster_positions: List[int], 
                               region_start: int, region_end: int,
                               mutations: List[str], gatekeeping_aas: List[str]) -> List[str]:
    """Get list of mutations to apply based on position type (edge vs internal)"""
    if is_edge_position(pos, cluster_positions, region_start, region_end):
        return list(set(mutations + gatekeeping_aas))
    return mutations

def is_edge_position(pos: int, cluster_positions: List[int], region_start: int, region_end: int) -> bool:
    """
    Check if a position is at the edge of a cluster or directly adjacent to it.
    
    Args:
        pos: Position to check
        cluster_positions: All positions in the cluster
        region_start: Start of the region (1-indexed)
        region_end: End of the region (1-indexed)
    
    Returns:
        True if position is at edge of cluster or adjacent to region boundary
    """
    if not cluster_positions:
        return False
    
    # Check if position is at the edge of the cluster
    min_pos = min(cluster_positions)
    max_pos = max(cluster_positions)
    
    # Position is at the edge if it's the minimum or maximum position in the cluster
    if pos == min_pos or pos == max_pos:
        return True
    
    # Check if position is adjacent to region boundaries
    # (within 1 position of region start or end)
    if pos == region_start or pos == region_end:
        return True
    
    # Check if position is directly adjacent to region boundaries
    if pos == region_start + 1 or pos == region_end - 1:
        return True
    
    return False

def apply_rule_mutations(sequence: str, region_analysis: Dict, mutations: List[str], 
                        selected_rules: List[str] = None, gatekeeping_aas: List[str] = None) -> List[Tuple[str, str, int]]:
    """Apply mutations based on rules in analyzed region"""
    results = []
    
    if gatekeeping_aas is None:
        gatekeeping_aas = GATEKEEPING_AAS
    
    # Filter rules if specific ones are selected
    if selected_rules:
        rules_to_apply = {k: v for k, v in region_analysis['rules'].items() if k in selected_rules}
    else:
        rules_to_apply = region_analysis['rules']
    
    # Get region boundaries
    region_start, region_end = region_analysis['region']
    
    # Collect all clusters from all rules
    all_clusters = []
    for rule_name, rule_data in rules_to_apply.items():
        if not rule_data['condition_met']:
            continue
        
        if rule_name == 'hydrophobic_and_aromatic' and 'special_clusters' in rule_data:
            for cluster in rule_data['special_clusters']:
                # Ensure cluster has rule_name field
                cluster_with_rule = cluster.copy()
                cluster_with_rule['rule_name'] = rule_name
                if 'aggregation_score' not in cluster_with_rule:
                    cluster_with_rule['aggregation_score'] = rule_data['aggregation_score']
                all_clusters.append((rule_name, cluster_with_rule))
        else:
            for cluster in rule_data['qualifying_clusters']:
                all_clusters.append((cluster['rule_name'], cluster))
    
    # Resolve overlapping clusters (keep union of motifs)
    non_overlapping_clusters = resolve_overlapping_clusters(all_clusters, sequence)
    
    # Track mutated positions to avoid duplicates (performance optimization)
    mutated_positions = set()
    
    # Apply mutations from non-overlapping clusters
    for rule_name, cluster in non_overlapping_clusters:
        positions_to_mutate = cluster['positions']
        
        # Check if this cluster is part of a multi-rule cluster
        is_part_of_multi = any(
            set(positions_to_mutate).issubset(set(mc['positions']))
            for mc in region_analysis['multi_rule_clusters']
        )
        
        # Get aggregation score for this cluster
        agg_score = cluster.get('aggregation_score') or cluster.get('combined_aggregation_score', 0)
        
        # Apply each mutation to each position in the cluster
        for pos in positions_to_mutate:
            original_aa = sequence[pos-1]
            
            # Determine which mutations to apply based on edge status
            all_mutations = get_mutations_for_position(
                pos, positions_to_mutate, region_start, region_end, 
                mutations, gatekeeping_aas
            )
            
            for new_aa in all_mutations:
                if new_aa == original_aa:
                    continue
                
                # Track this position as mutated
                mutated_positions.add(pos)
                
                mutated_seq = create_mutated_sequence(sequence, pos, new_aa)
                
                # Check if this is a gatekeeping mutation
                is_gatekeeping = new_aa in gatekeeping_aas and new_aa not in mutations
                
                # Build description based on rule type
                if '+' in rule_name:  # Merged multi-rule cluster
                    description = f"{original_aa}{pos}{new_aa} | MERGED RULES {rule_name} (agg_score={agg_score})"
                elif is_part_of_multi:
                    # Find the specific multi-rule cluster
                    for mc in region_analysis['multi_rule_clusters']:
                        if set(positions_to_mutate).issubset(set(mc['positions'])):
                            description = f"{original_aa}{pos}{new_aa} | MERGED RULES {'+'.join(mc['rules'])} (agg_score={mc['combined_aggregation_score']})"
                            agg_score = mc['combined_aggregation_score']  # Use the multi-rule score
                            break
                    else:
                        # Fallback if not found in multi-rule clusters
                        description = f"{original_aa}{pos}{new_aa} | Rule '{rule_name}' (agg_score={agg_score})"
                else:
                    # Standard rule-based mutation
                    description = f"{original_aa}{pos}{new_aa} | Rule '{rule_name}' (agg_score={agg_score})"
                
                # Add gatekeeping indicator if applicable
                if is_gatekeeping:
                    description += f" | GATEKEEPING ({new_aa})"
                
                results.append((description, mutated_seq, agg_score))
    
    # Apply mutations from multi-rule clusters (ensure no overlap with single-rule clusters)
    for multi_cluster in region_analysis['multi_rule_clusters']:
        positions_to_mutate = multi_cluster['positions']
        
        # Check if any position in this multi-rule cluster has already been mutated
        # Use set intersection for efficient checking (O(1) average case)
        if any(pos in mutated_positions for pos in positions_to_mutate):
            continue  # Skip if already mutated
        
        # Track all positions that will be mutated
        mutated_positions.update(positions_to_mutate)
        
        # Get aggregation score for this cluster
        agg_score = multi_cluster.get('combined_aggregation_score', 0)
        
        for pos in positions_to_mutate:
            original_aa = sequence[pos-1]
            
            # Determine which mutations to apply based on edge status
            all_mutations = get_mutations_for_position(
                pos, positions_to_mutate, region_start, region_end,
                mutations, gatekeeping_aas
            )
            
            for new_aa in all_mutations:
                if new_aa == original_aa:
                    continue
                
                mutated_seq = create_mutated_sequence(sequence, pos, new_aa)
                
                # Check if this is a gatekeeping mutation
                is_gatekeeping = new_aa in gatekeeping_aas and new_aa not in mutations
                
                # Create description for multi-rule mutation
                description = f"{original_aa}{pos}{new_aa} | MERGED RULES {'+'.join(multi_cluster['rules'])} (agg_score={agg_score})"
                
                # Add gatekeeping indicator if applicable
                if is_gatekeeping:
                    description += f" | GATEKEEPING ({new_aa})"
                
                results.append((description, mutated_seq, agg_score))
    
    return results

def mutate_sequence(sequence: str, positions: List[int], mutations: List[str], 
                   regions: List[str] = None, selected_rules: List[str] = None,
                   insertion_positions: List[int] = None, insertion_aas: List[str] = None,
                   gatekeeping_aas: List[str] = None, verbose: bool = False) -> Tuple[List[Tuple[str, str]], List[Dict]]:
    """Generate mutated sequences with point mutations, rule-based mutations, and insertions"""
    all_results = []  # Will store (description, sequence, agg_score)
    region_analyses = []
    seq_len = len(sequence)
    
    # Use default gatekeeping amino acids if not provided
    if gatekeeping_aas is None:
        gatekeeping_aas = GATEKEEPING_AAS
    
    # Validate positions for direct mutations
    for pos in positions:
        if pos < 1 or pos > seq_len:
            raise ValueError(f"Position {pos} is out of range (sequence length: {seq_len})")
    
    # Analyze regions for rule-based mutations
    if regions:
        for region_str in regions:
            start, stop = parse_region(region_str, seq_len)
            analysis = analyze_region(sequence, start, stop)
            region_analyses.append(analysis)
            
            if verbose:
                print(f"\n{'='*60}")
                print(f"Analyzing region {start}:{stop}:")
                print(f"{'='*60}")
                print(f"Sequence: {analysis['sequence']}")
                print(f"Length: {analysis['length']} residues")
                
                for rule_name, rule_data in analysis['rules'].items():
                    if rule_name == 'hydrophobic_and_aromatic' and rule_data['condition_met']:
                        print(f"\n{rule_name.upper()}:")
                        print(f"  Description: {rule_data['description']}")
                        print(f"  ✓ Rule triggered!")
                        
                        for i, cluster in enumerate(rule_data['special_clusters'], 1):
                            print(f"\n  Cluster {i} ({cluster['condition']}):")
                            print(f"    Positions: {cluster['positions']}")
                            print(f"    Residues: {''.join(cluster['residues'])}")
                            print(f"    Hydrophobic count: {cluster['hydrophobic_count']}")
                            print(f"    Aromatic count: {cluster['aromatic_count']}")
                            
                            if cluster['condition'] == 'at_least_2_pairs':
                                print(f"    Pairs found: {cluster['pair_count']}")
                                for j, pair in enumerate(cluster.get('pairs', []), 1):
                                    print(f"      Pair {j}: {pair['interaction']}")
                            else:  # 1_pair_plus_hydrophobic
                                if cluster.get('pairs'):
                                    pair = cluster['pairs'][0]
                                    print(f"    Pair: {pair['interaction']}")
                                if cluster.get('nearby_hydrophobics'):
                                    print(f"    Nearby hydrophobic residues:")
                                    for h in cluster['nearby_hydrophobics']:
                                        print(f"      {h['residue']}{h['position']} (distance: {h['distance']})")
                            
                            print(f"    Size: {cluster['size']}, Span: {cluster['span']}aa")
                    
                    elif rule_data['matching_positions']:
                        print(f"\n{rule_name.upper()}:")
                        print(f"  Description: {rule_data['description']}")
                        residues_str = ''.join(rule_data['matching_residues'])
                        print(f"  Matching residues: {residues_str}")
                        print(f"  Positions: {rule_data['matching_positions']}")
                        
                        if rule_data['clusters']:
                            print(f"  Found {len(rule_data['clusters'])} cluster(s):")
                            for i, cluster in enumerate(rule_data['clusters'], 1):
                                cluster_residues = [sequence[pos-1] for pos in cluster]
                                span = max(cluster) - min(cluster) + 1
                                print(f"    Cluster {i}: positions {cluster}, "
                                      f"residues {''.join(cluster_residues)}, "
                                      f"size={len(cluster)}, span={span}aa")
                        
                        if rule_data['qualifying_clusters']:
                            print(f"  ✓ QUALIFYING CLUSTERS:")
                            for i, cluster in enumerate(rule_data['qualifying_clusters'], 1):
                                print(f"    Cluster {i}: positions {cluster['positions']}, "
                                      f"residues {''.join(cluster['residues'])}, "
                                      f"size={cluster['size']}, span={cluster['span']}aa, "
                                      f"agg_score={cluster['aggregation_score']}")
                        else:
                            print(f"  ✗ No qualifying clusters")
                
                if analysis['multi_rule_clusters']:
                    print(f"\n{'*'*60}")
                    print(f"MULTI-RULE CLUSTERS (HIGH AGGREGATION RISK):")
                    print(f"{'*'*60}")
                    for i, cluster in enumerate(analysis['multi_rule_clusters'], 1):
                        print(f"\n  Multi-Rule Cluster {i}:")
                        print(f"    Rules: {', '.join(cluster['rules'])}")
                        print(f"    Positions: {cluster['positions']}")
                        print(f"    Residues: {''.join(cluster['residues'])}")
                        print(f"    Combined Aggregation Score: {cluster['combined_aggregation_score']}/8")
    
    # Apply rule-based mutations
    for analysis in region_analyses:
        rule_mutations = apply_rule_mutations(sequence, analysis, mutations, selected_rules, gatekeeping_aas)
        all_results.extend(rule_mutations)
    
    # Generate direct point mutations
    for pos in positions:
        original_aa = sequence[pos-1]
        
        # For direct mutations, apply only regular mutations (not gatekeeping)
        for new_aa in mutations:
            if new_aa == original_aa:
                continue
            
            mutated_seq = create_mutated_sequence(sequence, pos, new_aa)
            
            # Standardized description for direct mutations with agg_score=0
            description = f"{original_aa}{pos}{new_aa} | Direct mutation (agg_score=0)"
            all_results.append((description, mutated_seq, 0))
    
    # Generate insertions if specified
    if insertion_positions and insertion_aas:
        if len(insertion_positions) != len(insertion_aas):
            raise ValueError("insertion_positions and insertion_aas must have the same length")
        
        for ins_pos, ins_aa in zip(insertion_positions, insertion_aas):
            if ins_pos < 1 or ins_pos > seq_len + 1:
                raise ValueError(f"Insertion position {ins_pos} is out of range (1-{seq_len+1})")
            
            mutated_seq = sequence[:ins_pos-1] + ins_aa + sequence[ins_pos-1:]
            # Standardized description for insertions with agg_score=0
            description = f"Insertion: {ins_aa} inserted before position {ins_pos} (agg_score=0)"
            all_results.append((description, mutated_seq, 0))
    
    # Sort all results by aggregation score (descending)
    all_results.sort(key=lambda x: x[2], reverse=True)
    
    # Convert to format expected by downstream functions (description, sequence)
    sorted_results = [(desc, seq) for desc, seq, _ in all_results]
    
    return sorted_results, region_analyses

def write_fasta(output_file: str, original_header: str, original_seq: str, 
                mutations: List[Tuple[str, str]], include_original: bool = True):
    """Write results to FASTA file"""
    with open(output_file, 'w') as f:
        if include_original:
            f.write(f"{original_header}\n")
            for i in range(0, len(original_seq), FASTA_LINE_LENGTH):
                f.write(f"{original_seq[i:i+FASTA_LINE_LENGTH]}\n")
        
        for i, (description, mutated_seq) in enumerate(mutations, 1):
            # Extract protein name from original header (remove '>')
            protein_name = original_header[1:].strip()
            # Keep spaces in description for readability
            f.write(f">{protein_name}_{description}\n")
            for j in range(0, len(mutated_seq), FASTA_LINE_LENGTH):
                f.write(f"{mutated_seq[j:j+FASTA_LINE_LENGTH]}\n")

def setup_argument_parser() -> argparse.ArgumentParser:
    """Setup and return argument parser"""
    parser = argparse.ArgumentParser(
        description='Perform rule-based in silico mutagenesis on protein sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        usage='python mutagenesis.py <input_file> [options]'
    )
    
    # Required arguments
    required = parser.add_argument_group('REQUIRED ARGUMENTS')
    required.add_argument('input_file', nargs='?', help='Input FASTA file containing protein sequence')
    
    # Region and rule arguments
    region_args = parser.add_argument_group('REGION-BASED MUTAGENESIS')
    region_args.add_argument('-r', '--regions', type=str, nargs='+',
                            help='Regions to analyze for rule-based mutations (format: start:stop)')
    region_args.add_argument('--rules', type=str, nargs='+', choices=list(RULES.keys()),
                            help=f'Specific rules to apply (choices: {", ".join(RULES.keys())})')
    
    # Direct mutation arguments
    mutation_args = parser.add_argument_group('DIRECT MUTATIONS')
    mutation_args.add_argument('-p', '--positions', type=int, nargs='+', default=[],
                              help='Specific positions to mutate (1-indexed)')
    mutation_args.add_argument('-m', '--mutations', type=str, nargs='+', default=DEFAULT_MUTATIONS,
                              help=f'Amino acids to mutate to (default: {DEFAULT_MUTATIONS})')
    
    # Gatekeeping amino acids argument
    gatekeeping_args = parser.add_argument_group('GATEKEEPING AMINO ACIDS')
    gatekeeping_args.add_argument('-g', '--gatekeeping', type=str, nargs='+', default=GATEKEEPING_AAS,
                                 help=f'Amino acids to use as gatekeepers for edge positions (default: {GATEKEEPING_AAS})')
    
    # Insertion arguments
    insertion_args = parser.add_argument_group('INSERTIONS')
    insertion_args.add_argument('--insert-positions', type=int, nargs='+',
                               help='Positions for insertions (insertion happens BEFORE this position)')
    insertion_args.add_argument('--insert-aas', type=str, nargs='+',
                               help='Amino acids to insert at insertion positions')
    
    # Aggregation analysis arguments
    agg_args = parser.add_argument_group('AGGREGATION ANALYSIS')
    agg_args.add_argument('--agg-only', action='store_true',
                         help='Only identify aggregation hotspots without generating mutations')
    agg_args.add_argument('--min-agg-score', type=int, default=4,
                         help='Minimum aggregation score to flag as hotspot (default: 4)')
    
    # Output arguments
    output_args = parser.add_argument_group('OUTPUT')
    output_args.add_argument('-o', '--output', default='mutated_sequences.fasta',
                            help='Output FASTA file (default: mutated_sequences.fasta)')
    output_args.add_argument('--no-original', action='store_true',
                            help='Do not include original sequence in output')
    
    # Verbosity and help
    other_args = parser.add_argument_group('OTHER OPTIONS')
    other_args.add_argument('-v', '--verbose', action='store_true',
                           help='Show detailed analysis of regions including clusters')
    other_args.add_argument('-h', '--help', action='store_true',
                           help='Show this help message and exit')
    
    return parser

def print_help_info(parser: argparse.ArgumentParser):
    """Print detailed help information"""
    print("=" * 70)
    print("RULE-BASED IN SILICO MUTAGENESIS SCRIPT (WITH CLUSTERING)")
    print("=" * 70)
    print("\nDESCRIPTION:")
    print("Performs rule-based mutagenesis on protein sequences.")
    print("Rules apply ONLY when amino acids are clustered together.")
    print("Multiple rules can apply simultaneously to the same motif.")
    print("Motifs with multiple rule matches are flagged as aggregation hotspots.")
    print("\n" + "=" * 70)
    
    print_usage_example()
    parser.print_help()
    
    print("\n" + "=" * 70)
    print("RULE DETAILS AND AGGREGATION SCORES:")
    print("=" * 70)
    for rule_name, rule in RULES.items():
        print(f"\n{rule_name}:")
        print(f"  {rule['description']}")
        if rule_name == 'hydrophobic_and_aromatic':
            print(f"  Conditions:")
            print(f"    1. At least 2 hydrophobic-aromatic adjacent pairs")
            print(f"    2. OR 1 hydrophobic-aromatic pair + at least 1 hydrophobic within 3 positions")
        else:
            print(f"  Residues: {', '.join(sorted(rule['residues']))}")
        print(f"  Min cluster size: {rule['min_cluster_size']}")
        print(f"  Max gap: {rule['max_gap']} positions")
        print(f"  Aggregation score: {rule['aggregation_score']}")
    
    print("\n" + "=" * 70)
    print("OVERLAP HANDLING PROCEDURE:")
    print("=" * 70)
    print("Overlapping motifs are resolved by keeping the UNION of motif positions.")
    print("Additionally, motifs within 2 positions of each other are also merged.")
    print("When clusters overlap or are close (≤2 positions apart), the script will:")
    print("1. Merge all overlapping clusters into a single unified cluster")
    print("2. Merge clusters that are within 2 positions of each other")
    print("3. Apply mutations to the unified motif")
    print("This creates larger mutagenesis regions covering all adjacent risk areas.")
    print("=" * 70)
    
    print("\n" + "=" * 70)
    print("GATEKEEPING AMINO ACIDS:")
    print("=" * 70)
    print("Gatekeeping amino acids are only applied to positions at the edge of motifs")
    print("or directly adjacent to region boundaries (within 1 position).")
    print(f"Default gatekeeping amino acids: {GATEKEEPING_AAS}")
    print("Use --gatekeeping option to specify custom gatekeeping amino acids.")
    print("=" * 70)

def validate_arguments(args) -> None:
    """Validate command line arguments"""
    # Check if at least one type of mutation is requested
    if not args.agg_only and not args.positions and not args.regions and not args.insert_positions:
        print("\nERROR: You must specify at least one of:")
        print("  • --positions for direct mutations")
        print("  • --regions for rule-based mutations")
        print("  • --insert-positions for insertions")
        print("  • --agg-only for aggregation analysis only")
        print("\nUse --help for more information.")
        sys.exit(1)
    
    # Validate mutation list
    try:
        validate_amino_acids(args.mutations, "mutation amino acids", strict=True)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate gatekeeping amino acids list
    try:
        validate_amino_acids(args.gatekeeping, "gatekeeping amino acids", strict=True)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate that insert-positions and insert-aas are both provided or both omitted
    if bool(args.insert_positions) != bool(args.insert_aas):
        print("\nERROR: Both --insert-positions and --insert-aas must be provided together", file=sys.stderr)
        print("Example: --insert-positions 10 20 --insert-aas K M")
        sys.exit(1)

def print_mutation_summary(mutations: List[Tuple[str, str]]):
    """Print summary of generated mutations"""
    if not mutations:
        print("\nNo mutations generated. Check your criteria.")
        return
    
    print(f"\n{'='*70}")
    print("MUTATION SUMMARY (Sorted by aggregation score)")
    print(f"{'='*70}")
    
    # First, extract aggregation scores and group mutations
    mutation_data = []
    for desc, seq in mutations:
        # Extract agg_score from description
        agg_score = 0
        if "(agg_score=" in desc:
            try:
                agg_str = desc.split("(agg_score=")[1].split(")")[0]
                agg_score = int(agg_str)
            except (IndexError, ValueError):
                agg_score = 0
        mutation_data.append((desc, seq, agg_score))
    
    # Count by type using the new standardized headers
    direct = sum(1 for d, _, _ in mutation_data if "Direct mutation" in d)
    rule_based = sum(1 for d, _, _ in mutation_data if "Rule '" in d and "GATEKEEPING" not in d)
    merged_rules = sum(1 for d, _, _ in mutation_data if "MERGED RULES" in d and "GATEKEEPING" not in d)
    insertions = sum(1 for d, _, _ in mutation_data if "Insertion" in d)
    gatekeeping = sum(1 for d, _, _ in mutation_data if "GATEKEEPING" in d)
    
    # Count specific rule types for detailed breakdown
    hydrophobic = sum(1 for d, _, _ in mutation_data if "Rule 'hydrophobic_aliphatic'" in d and "GATEKEEPING" not in d)
    aromatic = sum(1 for d, _, _ in mutation_data if "Rule 'aromatic'" in d and "GATEKEEPING" not in d)
    amide = sum(1 for d, _, _ in mutation_data if "Rule 'amide'" in d and "GATEKEEPING" not in d)
    hydrophobic_arom = sum(1 for d, _, _ in mutation_data if "Rule 'hydrophobic_and_aromatic'" in d and "GATEKEEPING" not in d)
    
    print(f"Direct mutations: {direct} (agg_score: 0)")
    print(f"Rule-based mutations: {rule_based}")
    print(f"  • Hydrophobic mutations: {hydrophobic} (agg_score: 3)")
    print(f"  • Hydrophobic-aromatic mutations: {hydrophobic_arom} (agg_score: 2)")
    print(f"  • Aromatic mutations: {aromatic} (agg_score: 2)")
    print(f"  • Amide mutations: {amide} (agg_score: 1)")
    print(f"Merged rule mutations: {merged_rules} (agg_score: 4-8)")
    print(f"Insertions: {insertions} (agg_score: 0)")
    print(f"Gatekeeping mutations: {gatekeeping}")
    print(f"TOTAL: {len(mutations)}")
    
    # Show top 5 mutations by aggregation score
    if mutation_data:
        print(f"\nTOP 5 MUTATIONS BY AGGREGATION SCORE:")
        for i, (desc, _, agg_score) in enumerate(mutation_data[:5], 1):
            if len(desc) > 80:
                desc = desc[:77] + "..."
            print(f"{i}. {desc}")
    
    if hydrophobic_arom > 0:
        print(f"\nHydrophobic-aromatic mutations found:")
        for desc, _, _ in mutation_data:
            if "Rule 'hydrophobic_and_aromatic'" in desc and "GATEKEEPING" not in desc:
                if len(desc) > 80:
                    desc = desc[:77] + "..."
                print(f"  • {desc}")
                break
    
    if gatekeeping > 0:
        print(f"\nGatekeeping mutations (edge positions only):")
        gatekeeping_shown = 0
        for i, (desc, _, _) in enumerate(mutation_data):
            if "GATEKEEPING" in desc:
                if len(desc) > 80:
                    desc = desc[:77] + "..."
                print(f"  {gatekeeping_shown+1}. {desc}")
                gatekeeping_shown += 1
                if gatekeeping_shown >= 5:  # Show first 5 gatekeeping mutations
                    if gatekeeping > 5:
                        print(f"  ... and {gatekeeping - 5} more gatekeeping mutations")
                    break
    
    print(f"\nAll mutations (first 10, sorted by agg_score):")
    shown = 0
    for i, (desc, _, agg_score) in enumerate(mutation_data):
        print(f"{i+1}. {desc}")
        shown += 1
        if shown >= 10:
            if len(mutation_data) > 10:
                print(f"... and {len(mutation_data) - 10} more mutations")
            break

def main():
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # Show help if requested or no input file provided
    if args.help or not args.input_file:
        print_help_info(parser)
        if not args.input_file:
            print("\nERROR: Input FASTA file is required!")
            sys.exit(1)
        sys.exit(0)
    
    # Validate arguments
    validate_arguments(args)
    
    try:
        # Read and validate input sequence
        print(f"\n{'='*70}")
        print("RULE-BASED MUTAGENESIS WITH UNION OVERLAP RESOLUTION")
        print(f"{'='*70}")
        print(f"Input file: {args.input_file}")
        
        header, sequence = read_fasta(args.input_file)
        
        if not validate_sequence(sequence):
            print("Warning: Sequence contains non-standard amino acid codes", file=sys.stderr)
        
        print(f"Original sequence length: {len(sequence)} amino acids")
        
        if args.regions:
            print(f"Regions to analyze: {args.regions}")
        
        if args.positions:
            print(f"Direct mutation positions: {args.positions}")
        
        print(f"Regular mutations to apply: {args.mutations}")
        print(f"Gatekeeping amino acids: {args.gatekeeping}")
        
        if args.rules:
            print(f"Rules selected: {args.rules}")
        else:
            print(f"Rules selected: All rules")
        
        if args.agg_only:
            print(f"\nMODE: Aggregation hotspot analysis only")
            print(f"Minimum aggregation score: {args.min_agg_score}")
        
        print(f"\n{'='*70}")
        print("ANALYZING REGIONS FOR AMYLOIDOGENIC MOTIFS...")
        print(f"{'='*70}")
        
        # Generate mutations or analyze
        if not args.agg_only:
            mutations, region_analyses = mutate_sequence(
                sequence, 
                args.positions,
                [m.upper() for m in args.mutations],
                args.regions,
                args.rules,
                args.insert_positions,
                [aa.upper() for aa in args.insert_aas] if args.insert_aas else None,
                [aa.upper() for aa in args.gatekeeping],
                args.verbose
            )
            
            # NOTE: mutations are now sorted by aggregation score (descending)
            print(f"\n✓ Generated {len(mutations)} mutated sequences (sorted by aggregation score)")
        else:
            region_analyses = []
            mutations = []
            if args.regions:
                for region_str in args.regions:
                    start, stop = parse_region(region_str, len(sequence))
                    analysis = analyze_region(sequence, start, stop)
                    region_analyses.append(analysis)
        
        # MODIFIED SECTION: Print region summary in the requested format
        if region_analyses and not args.verbose:
            print("\nREGION ANALYSIS SUMMARY:")
            print("-" * 50)
            
            # Track overlapping clusters across regions
            all_region_clusters = []
            
            for analysis in region_analyses:
                start, stop = analysis['region']
                print(f"\nRegion {start}:{stop} (length: {analysis['length']}):")
                
                triggered_rules = False
                for rule_name, rule_data in analysis['rules'].items():
                    if rule_name == 'hydrophobic_and_aromatic':
                        # Special handling for hydrophobic_and_aromatic rule
                        if rule_data['condition_met'] and rule_data.get('special_clusters'):
                            triggered_rules = True
                            print(f"  ✓ {rule_name}: {len(rule_data['special_clusters'])} qualifying cluster(s)")
                            for cluster in rule_data['special_clusters']:
                                print(f"      • Positions: {cluster['positions']}")
                                print(f"        Residues: {''.join(cluster['residues'])}")
                                print(f"        Size: {cluster['size']}, Span: {cluster['span']}aa")
                                print(f"        Agg score: {cluster.get('aggregation_score', rule_data['aggregation_score'])}")
                                if cluster['condition'] == 'at_least_2_pairs':
                                    print(f"        Condition: {cluster['pair_count']} hydrophobic-aromatic pairs")
                                else:
                                    print(f"        Condition: 1 pair + {len(cluster.get('nearby_hydrophobics', []))} hydrophobic(s)")
                                
                                # Store for overlap analysis
                                all_region_clusters.append({
                                    'region': (start, stop),
                                    'rule': rule_name,
                                    'positions': cluster['positions'],
                                    'size': cluster['size']
                                })
                    elif rule_data['qualifying_clusters']:
                        triggered_rules = True
                        print(f"  ✓ {rule_name}: {len(rule_data['qualifying_clusters'])} qualifying cluster(s)")
                        for cluster in rule_data['qualifying_clusters']:
                            print(f"      • Positions: {cluster['positions']}")
                            print(f"        Residues: {''.join(cluster['residues'])}")
                            print(f"        Size: {cluster['size']}, Span: {cluster['span']}aa")
                            print(f"        Agg score: {cluster['aggregation_score']}")
                            
                            # Store for overlap analysis
                            all_region_clusters.append({
                                'region': (start, stop),
                                'rule': rule_name,
                                'positions': cluster['positions'],
                                'size': cluster['size']
                            })
                
                if not triggered_rules:
                    print("  ✗ No rules triggered (no qualifying clusters found)")
            
            # Check for overlaps across regions
            if len(all_region_clusters) > 1:
                # Sort clusters by size (largest first)
                sorted_clusters = sorted(all_region_clusters, key=lambda x: x['size'], reverse=True)
                
                # Find overlapping or close clusters
                overlapping_clusters = []
                
                for i, cluster1 in enumerate(sorted_clusters):
                    for j, cluster2 in enumerate(sorted_clusters):
                        if i >= j:
                            continue
                        
                        set1 = set(cluster1['positions'])
                        set2 = set(cluster2['positions'])
                        
                        # Calculate the gap between clusters
                        min1, max1 = min(set1), max(set1)
                        min2, max2 = min(set2), max(set2)
                        
                        if min1 < min2:
                            gap = min2 - max1 - 1
                        else:
                            gap = min1 - max2 - 1
                        
                        # Check for overlap OR gap ≤ 2
                        if set1.intersection(set2) or gap <= 2:
                            overlapping_clusters.append((cluster1, cluster2, gap))
                
                if overlapping_clusters:
                    print(f"\n{'!'*70}")
                    print(f"OVERLAPPING/CLOSE MOTIFS DETECTED:")
                    print(f"{'!'*70}")
                    print(f"Found {len(overlapping_clusters)} overlapping/close cluster pairs")
                    print("Motifs will be merged into unified mutagenesis regions.")
                    
                    for cluster1, cluster2, gap in overlapping_clusters[:3]:  # Show first 3 overlaps
                        print(f"\n{'Overlap' if gap < 0 else 'Proximity'} detected (gap: {abs(gap)}):")
                        print(f"  • Region {cluster1['region'][0]}:{cluster1['region'][1]} - {cluster1['rule']}")
                        print(f"    Positions: {cluster1['positions']} (size: {cluster1['size']})")
                        print(f"  • Region {cluster2['region'][0]}:{cluster2['region'][1]} - {cluster2['rule']}")
                        print(f"    Positions: {cluster2['positions']} (size: {cluster2['size']})")
                        print(f"  → Will create unified cluster covering {len(set(cluster1['positions']).union(set(cluster2['positions'])))} positions")
                    
                    if len(overlapping_clusters) > 3:
                        print(f"  ... and {len(overlapping_clusters) - 3} more overlapping/close pairs")
        
        if not args.agg_only:
            write_fasta(args.output, header, sequence, mutations, not args.no_original)
            print(f"✓ Results written to {args.output}")
            
            print_mutation_summary(mutations)
        else:
            print(f"\n✓ Aggregation analysis completed")
        
        print(f"\n{'='*70}")
        print("PROCESS COMPLETED SUCCESSFULLY")
        print(f"{'='*70}")
            
    except Exception as e:
        print(f"\n{'='*70}")
        print("ERROR OCCURRED")
        print(f"{'='*70}")
        print(f"Error: {e}", file=sys.stderr)
        print(f"\nFor help, use: python {sys.argv[0]} --help")
        sys.exit(1)

if __name__ == '__main__':
    main()
