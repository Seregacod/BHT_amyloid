#!/usr/bin/env python3
"""
Extract Amyloid Region Sequences (Enhanced)

This script reads the unified amyloid database and generates FASTA files
with actual region sequences fetched from UniProt and NCBI databases.

Features:
- Fetches sequences from UniProt API (by accession)
- Fetches sequences from NCBI Protein database (by protein name search)
- Extracts specific amyloidogenic regions based on start/end coordinates
- Generates separate FASTA files for regions and full proteins
- Comprehensive caching to avoid redundant API calls

Output files:
- {prefix}_regions.fasta     - Amyloidogenic region sequences
- {prefix}_full_proteins.fasta - Full protein sequences
- {prefix}_cache.json        - Sequence cache for reuse

Author: Claude (Anthropic)
Date: 2025
"""

import csv
import re
import time
import json
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
from urllib.parse import quote, urlencode
import xml.etree.ElementTree as ET


# =============================================================================
# Configuration
# =============================================================================

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

REQUEST_DELAY = 0.4  # seconds between API requests (NCBI requires 0.33s min)
TIMEOUT = 30  # request timeout in seconds


# =============================================================================
# Sequence Fetching Functions
# =============================================================================

def fetch_uniprot_sequence(uniprot_id: str, verbose: bool = False) -> Optional[str]:
    """
    Fetch protein sequence from UniProt REST API.
    
    Args:
        uniprot_id: UniProt accession (e.g., P05067)
        verbose: Print progress messages
    
    Returns:
        Protein sequence or None if not found
    """
    if not uniprot_id:
        return None
    
    # Clean UniProt ID
    uniprot_id = uniprot_id.strip().upper()
    if not re.match(r'^[A-Z][A-Z0-9]{5,9}$', uniprot_id):
        return None
    
    try:
        url = f"{UNIPROT_API}/{uniprot_id}.fasta"
        req = Request(url, headers={'User-Agent': 'Python/AmyloidDB-Fetcher'})
        
        with urlopen(req, timeout=TIMEOUT) as response:
            content = response.read().decode('utf-8')
            
            # Parse FASTA
            lines = content.strip().split('\n')
            if lines and lines[0].startswith('>'):
                sequence = ''.join(lines[1:]).replace(' ', '')
                if verbose:
                    print(f"    UniProt {uniprot_id}: {len(sequence)} aa")
                return sequence
    
    except HTTPError as e:
        if e.code == 404:
            if verbose:
                print(f"    UniProt {uniprot_id}: not found")
        else:
            if verbose:
                print(f"    UniProt {uniprot_id}: HTTP {e.code}")
    except URLError as e:
        if verbose:
            print(f"    UniProt {uniprot_id}: network error")
    except Exception as e:
        if verbose:
            print(f"    UniProt {uniprot_id}: error - {e}")
    
    return None


def search_ncbi_protein(query: str, organism: str = None, verbose: bool = False) -> Optional[str]:
    """
    Search NCBI Protein database and return first matching GI/accession.
    
    Args:
        query: Protein name to search
        organism: Optional organism filter
        verbose: Print progress messages
    
    Returns:
        NCBI Protein ID or None
    """
    if not query:
        return None
    
    # Build search query
    search_term = f'"{query}"[Protein Name]'
    if organism:
        search_term += f' AND "{organism}"[Organism]'
    
    params = {
        'db': 'protein',
        'term': search_term,
        'retmax': 1,
        'retmode': 'json',
        'sort': 'relevance'
    }
    
    try:
        url = f"{NCBI_ESEARCH}?{urlencode(params)}"
        req = Request(url, headers={'User-Agent': 'Python/AmyloidDB-Fetcher'})
        
        with urlopen(req, timeout=TIMEOUT) as response:
            data = json.loads(response.read().decode('utf-8'))
            
            id_list = data.get('esearchresult', {}).get('idlist', [])
            if id_list:
                return id_list[0]
    
    except Exception as e:
        if verbose:
            print(f"    NCBI search '{query[:30]}...': error - {e}")
    
    return None


def fetch_ncbi_sequence(protein_id: str, verbose: bool = False) -> Optional[Tuple[str, str]]:
    """
    Fetch protein sequence from NCBI by ID.
    
    Args:
        protein_id: NCBI Protein ID (GI or accession)
        verbose: Print progress messages
    
    Returns:
        Tuple of (accession, sequence) or None
    """
    if not protein_id:
        return None
    
    params = {
        'db': 'protein',
        'id': protein_id,
        'rettype': 'fasta',
        'retmode': 'text'
    }
    
    try:
        url = f"{NCBI_EFETCH}?{urlencode(params)}"
        req = Request(url, headers={'User-Agent': 'Python/AmyloidDB-Fetcher'})
        
        with urlopen(req, timeout=TIMEOUT) as response:
            content = response.read().decode('utf-8')
            
            lines = content.strip().split('\n')
            if lines and lines[0].startswith('>'):
                # Extract accession from header
                header = lines[0]
                acc_match = re.search(r'>(\S+)', header)
                accession = acc_match.group(1) if acc_match else protein_id
                
                sequence = ''.join(lines[1:]).replace(' ', '')
                if verbose:
                    print(f"    NCBI {protein_id}: {len(sequence)} aa")
                return (accession, sequence)
    
    except Exception as e:
        if verbose:
            print(f"    NCBI {protein_id}: error - {e}")
    
    return None


def fetch_sequence_by_name(protein_name: str, organism: str = None, 
                           verbose: bool = False) -> Optional[Tuple[str, str]]:
    """
    Search and fetch sequence from NCBI by protein name.
    
    Args:
        protein_name: Protein name to search
        organism: Optional organism filter
        verbose: Print progress messages
    
    Returns:
        Tuple of (accession, sequence) or None
    """
    # Clean protein name
    name = protein_name.strip()
    if not name or name.upper() in ['N.A.', 'NA', 'UNKNOWN', '']:
        return None
    
    # Remove common suffixes that might interfere with search
    name = re.sub(r'_HUMAN$|_MOUSE$|_RAT$|_YEAST$', '', name)
    
    # Search NCBI
    ncbi_id = search_ncbi_protein(name, organism, verbose)
    if ncbi_id:
        time.sleep(REQUEST_DELAY)  # Rate limiting
        return fetch_ncbi_sequence(ncbi_id, verbose)
    
    return None


# =============================================================================
# Cache Management
# =============================================================================

class SequenceCache:
    """Manages sequence caching to/from JSON file."""
    
    def __init__(self, filepath: str = None):
        self.filepath = filepath
        self.uniprot: Dict[str, Optional[str]] = {}  # uniprot_id -> sequence
        self.ncbi: Dict[str, Optional[Tuple[str, str]]] = {}  # protein_name -> (acc, seq)
        self.failed: Set[str] = set()  # IDs that failed to fetch
        
        if filepath and Path(filepath).exists():
            self.load()
    
    def load(self):
        """Load cache from JSON file."""
        try:
            with open(self.filepath, 'r') as f:
                data = json.load(f)
                self.uniprot = data.get('uniprot', {})
                self.ncbi = {k: tuple(v) if v else None 
                            for k, v in data.get('ncbi', {}).items()}
                self.failed = set(data.get('failed', []))
                print(f"  Loaded cache: {len(self.uniprot)} UniProt, {len(self.ncbi)} NCBI")
        except Exception as e:
            print(f"  Warning: Could not load cache - {e}")
    
    def save(self):
        """Save cache to JSON file."""
        if not self.filepath:
            return
        
        try:
            data = {
                'uniprot': {k: v for k, v in self.uniprot.items() if v is not None},
                'ncbi': {k: list(v) if v else None for k, v in self.ncbi.items()},
                'failed': list(self.failed)
            }
            with open(self.filepath, 'w') as f:
                json.dump(data, f, indent=2)
            print(f"  Saved cache: {len(data['uniprot'])} UniProt, {len(data['ncbi'])} NCBI")
        except Exception as e:
            print(f"  Warning: Could not save cache - {e}")
    
    def get_uniprot(self, uniprot_id: str) -> Optional[str]:
        """Get sequence from UniProt cache."""
        return self.uniprot.get(uniprot_id)
    
    def set_uniprot(self, uniprot_id: str, sequence: Optional[str]):
        """Set sequence in UniProt cache."""
        self.uniprot[uniprot_id] = sequence
        if sequence is None:
            self.failed.add(f"uniprot:{uniprot_id}")
    
    def has_uniprot(self, uniprot_id: str) -> bool:
        """Check if UniProt ID is in cache."""
        return uniprot_id in self.uniprot
    
    def get_ncbi(self, protein_name: str) -> Optional[Tuple[str, str]]:
        """Get sequence from NCBI cache."""
        return self.ncbi.get(protein_name.lower())
    
    def set_ncbi(self, protein_name: str, result: Optional[Tuple[str, str]]):
        """Set sequence in NCBI cache."""
        self.ncbi[protein_name.lower()] = result
        if result is None:
            self.failed.add(f"ncbi:{protein_name}")
    
    def has_ncbi(self, protein_name: str) -> bool:
        """Check if protein name is in NCBI cache."""
        return protein_name.lower() in self.ncbi
    
    def is_failed(self, key: str) -> bool:
        """Check if this ID previously failed to fetch."""
        return key in self.failed


# =============================================================================
# Sequence Extraction
# =============================================================================

def extract_region(full_sequence: str, start: int, end: int) -> Optional[str]:
    """
    Extract region from full sequence using 1-based coordinates.
    
    Args:
        full_sequence: Full protein sequence
        start: 1-based start position
        end: 1-based end position (inclusive)
    
    Returns:
        Extracted subsequence or None if invalid
    """
    if not full_sequence or not start or not end:
        return None
    
    # Convert to 0-based indexing
    start_idx = start - 1
    end_idx = end  # Python slicing is exclusive, so end stays as-is
    
    # Validate coordinates
    if start_idx < 0:
        start_idx = 0
    if end_idx > len(full_sequence):
        end_idx = len(full_sequence)
    if start_idx >= end_idx:
        return None
    
    return full_sequence[start_idx:end_idx]


def clean_sequence(seq: str) -> str:
    """Clean and validate amino acid sequence."""
    if not seq:
        return ""
    seq = seq.upper().strip()
    # Remove non-amino acid characters
    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq)
    return seq


# =============================================================================
# Main Processing
# =============================================================================

def load_unified_database(filepath: str) -> List[Dict]:
    """Load unified amyloid database."""
    entries = []
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            entries.append(dict(row))
    return entries


def fetch_all_sequences(entries: List[Dict], cache: SequenceCache,
                        fetch_uniprot: bool = True, fetch_ncbi: bool = True,
                        verbose: bool = False) -> Dict[str, str]:
    """
    Fetch all sequences needed for the database entries.
    
    Returns:
        Dictionary mapping entry key to full protein sequence
    """
    sequences = {}  # (source, id) -> sequence
    
    # Collect unique IDs to fetch
    uniprot_ids = set()
    ncbi_queries = []  # (protein_name, organism)
    
    for entry in entries:
        # Skip entries without coordinates
        if not entry.get('region_start') or not entry.get('region_end'):
            continue
        
        uniprot_id = entry.get('uniprot_id', '').strip()
        protein_name = entry.get('protein_name', '').strip()
        organism = entry.get('organism', '').strip()
        
        if uniprot_id and uniprot_id not in ['N.A.', 'NA', '']:
            uniprot_ids.add(uniprot_id)
        elif protein_name and protein_name not in ['N.A.', 'NA', '']:
            ncbi_queries.append((protein_name, organism))
    
    # Fetch from UniProt
    if fetch_uniprot and uniprot_ids:
        to_fetch = [uid for uid in uniprot_ids if not cache.has_uniprot(uid)]
        if to_fetch:
            print(f"\nFetching {len(to_fetch)} sequences from UniProt...")
            for i, uid in enumerate(to_fetch):
                if (i + 1) % 50 == 0:
                    print(f"  Progress: {i+1}/{len(to_fetch)}")
                
                seq = fetch_uniprot_sequence(uid, verbose)
                cache.set_uniprot(uid, seq)
                time.sleep(REQUEST_DELAY)
    
    # Fetch from NCBI (for entries without UniProt ID)
    if fetch_ncbi and ncbi_queries:
        # Deduplicate
        seen = set()
        unique_queries = []
        for name, org in ncbi_queries:
            key = name.lower()
            if key not in seen and not cache.has_ncbi(name):
                seen.add(key)
                unique_queries.append((name, org))
        
        if unique_queries:
            print(f"\nFetching {len(unique_queries)} sequences from NCBI...")
            for i, (name, org) in enumerate(unique_queries):
                if (i + 1) % 20 == 0:
                    print(f"  Progress: {i+1}/{len(unique_queries)}")
                
                result = fetch_sequence_by_name(name, org, verbose)
                cache.set_ncbi(name, result)
                time.sleep(REQUEST_DELAY)
    
    return cache


def process_entries(entries: List[Dict], cache: SequenceCache,
                    amyloid_only: bool = True, min_length: int = 4,
                    include_direct: bool = True) -> Tuple[Dict, Dict]:
    """
    Process database entries and extract sequences.
    
    Returns:
        Tuple of (regions_dict, full_proteins_dict)
        Each dict maps sequence to (header, metadata)
    """
    regions = {}  # seq -> (header, entry)
    full_proteins = {}  # seq -> (header, entry)
    
    stats = defaultdict(int)
    
    for entry in entries:
        # Filter by amyloid status
        if amyloid_only and entry.get('is_amyloid', '').lower() != 'yes':
            continue
        
        stats['total_amyloid'] += 1
        
        region_seq = None
        full_seq = None
        source = None
        
        # Priority 1: Direct peptide sequence (for short peptides)
        direct_seq = clean_sequence(entry.get('sequence', ''))
        if include_direct and direct_seq and len(direct_seq) >= min_length:
            region_seq = direct_seq
            source = 'direct'
            stats['from_direct'] += 1
        
        # Priority 2: Extract from UniProt
        uniprot_id = entry.get('uniprot_id', '').strip()
        start = entry.get('region_start')
        end = entry.get('region_end')
        
        if uniprot_id and start and end:
            try:
                start = int(start)
                end = int(end)
                
                full_seq = cache.get_uniprot(uniprot_id)
                if full_seq:
                    extracted = extract_region(full_seq, start, end)
                    if extracted and len(extracted) >= min_length:
                        region_seq = extracted
                        source = 'uniprot'
                        stats['from_uniprot'] += 1
            except (ValueError, TypeError):
                stats['invalid_coords'] += 1
        
        # Priority 3: Extract from NCBI (if no UniProt)
        if not region_seq and not uniprot_id:
            protein_name = entry.get('protein_name', '').strip()
            if protein_name and start and end:
                try:
                    start = int(start)
                    end = int(end)
                    
                    ncbi_result = cache.get_ncbi(protein_name)
                    if ncbi_result:
                        acc, full_seq = ncbi_result
                        extracted = extract_region(full_seq, start, end)
                        if extracted and len(extracted) >= min_length:
                            region_seq = extracted
                            source = 'ncbi'
                            stats['from_ncbi'] += 1
                except (ValueError, TypeError):
                    stats['invalid_coords'] += 1
        
        if not region_seq:
            stats['no_sequence'] += 1
            continue
        
        # Build headers
        parts = [entry.get('source_db', 'Unknown')]
        if entry.get('record_id'):
            parts.append(entry['record_id'])
        if entry.get('protein_name'):
            parts.append(entry['protein_name'].replace(' ', '_').replace('|', '-'))
        if entry.get('uniprot_id'):
            parts.append(f"UniProt:{entry['uniprot_id']}")
        if start and end:
            parts.append(f"region:{start}-{end}")
        if entry.get('category'):
            parts.append(f"cat:{entry['category']}")
        parts.append(f"src:{source}")
        
        region_header = "|".join(parts)
        
        # Store region (deduplicate by sequence)
        if region_seq not in regions:
            regions[region_seq] = (region_header, entry)
        
        # Store full protein
        if full_seq and full_seq not in full_proteins:
            protein_header_parts = []
            if entry.get('uniprot_id'):
                protein_header_parts.append(f"UniProt:{entry['uniprot_id']}")
            if entry.get('protein_name'):
                protein_header_parts.append(entry['protein_name'].replace(' ', '_'))
            if entry.get('organism'):
                protein_header_parts.append(f"[{entry['organism']}]")
            
            protein_header = "|".join(protein_header_parts) if protein_header_parts else f"protein_{len(full_proteins)}"
            full_proteins[full_seq] = (protein_header, entry)
    
    # Print statistics
    print("\n" + "="*60)
    print("SEQUENCE EXTRACTION STATISTICS")
    print("="*60)
    for key, count in sorted(stats.items()):
        print(f"  {key}: {count}")
    print(f"\n  Unique region sequences: {len(regions)}")
    print(f"  Unique full proteins: {len(full_proteins)}")
    print("="*60)
    
    return regions, full_proteins


def write_fasta(sequences: Dict[str, Tuple[str, Dict]], output_path: str,
                wrap_width: int = 60):
    """Write sequences to FASTA file."""
    with open(output_path, 'w') as f:
        for seq, (header, _) in sorted(sequences.items(), key=lambda x: x[1][0]):
            f.write(f">{header}\n")
            for i in range(0, len(seq), wrap_width):
                f.write(f"{seq[i:i+wrap_width]}\n")
    
    print(f"Wrote {len(sequences)} sequences to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Extract amyloid region sequences from unified database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output files:
  {prefix}_regions.fasta       - Amyloidogenic region sequences
  {prefix}_full_proteins.fasta - Full parent protein sequences  
  {prefix}_cache.json          - Sequence cache for reuse

Example usage:
  # Fetch from both UniProt and NCBI
  python extract_amyloid_sequences.py unified_amyloid_db.tsv \\
      -o amyloid_sequences --fetch-all -v

  # Use existing cache (no network requests)
  python extract_amyloid_sequences.py unified_amyloid_db.tsv \\
      -o amyloid_sequences --cache existing_cache.json

  # Fetch only from UniProt
  python extract_amyloid_sequences.py unified_amyloid_db.tsv \\
      -o amyloid_sequences --fetch-uniprot
        """
    )
    
    parser.add_argument('-i', '--input', help='Unified amyloid database TSV file')
    parser.add_argument('-o', '--output', default='amyloid',
                        help='Output file prefix (default: amyloid)')
    parser.add_argument('--cache', help='Sequence cache JSON file (default: {prefix}_cache.json)')
    
    # Fetching options
    fetch_group = parser.add_argument_group('Fetching options')
    fetch_group.add_argument('--fetch-all', action='store_true',
                             help='Fetch from both UniProt and NCBI')
    fetch_group.add_argument('--fetch-uniprot', action='store_true',
                             help='Fetch sequences from UniProt API')
    fetch_group.add_argument('--fetch-ncbi', action='store_true',
                             help='Fetch sequences from NCBI API')
    
    # Filtering options
    filter_group = parser.add_argument_group('Filtering options')
    filter_group.add_argument('--all', action='store_true',
                              help='Include all entries (not just amyloid-positive)')
    filter_group.add_argument('--no-direct', action='store_true',
                              help='Exclude direct peptide sequences')
    filter_group.add_argument('--min-length', type=int, default=4,
                              help='Minimum sequence length (default: 4)')
    
    # Other options
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    parser.add_argument('--stats-only', action='store_true',
                        help='Only print statistics, no output files')
    
    args = parser.parse_args()
    
    # Determine cache file path
    cache_path = args.cache or f"{args.output}_cache.json"
    
    # Load data
    print(f"Loading {args.input}...")
    entries = load_unified_database(args.input)
    print(f"  Loaded {len(entries)} entries")
    
    # Initialize cache
    cache = SequenceCache(cache_path)
    
    # Determine what to fetch
    fetch_uniprot = args.fetch_uniprot or args.fetch_all
    fetch_ncbi = args.fetch_ncbi or args.fetch_all
    
    # Fetch sequences if requested
    if fetch_uniprot or fetch_ncbi:
        fetch_all_sequences(
            entries, cache,
            fetch_uniprot=fetch_uniprot,
            fetch_ncbi=fetch_ncbi,
            verbose=args.verbose
        )
        cache.save()
    
    if args.stats_only:
        # Just print cache statistics
        print(f"\nCache contains:")
        print(f"  UniProt sequences: {len([v for v in cache.uniprot.values() if v])}")
        print(f"  NCBI sequences: {len([v for v in cache.ncbi.values() if v])}")
        print(f"  Failed lookups: {len(cache.failed)}")
        return
    
    # Process entries
    regions, full_proteins = process_entries(
        entries, cache,
        amyloid_only=not args.all,
        min_length=args.min_length,
        include_direct=not args.no_direct
    )
    
    # Write output files
    regions_file = f"{args.output}_regions.fasta"
    proteins_file = f"{args.output}_full_proteins.fasta"
    
    write_fasta(regions, regions_file)
    write_fasta(full_proteins, proteins_file)
    
    print(f"\nDone! Output files:")
    print(f"  {regions_file}")
    print(f"  {proteins_file}")
    print(f"  {cache_path}")


if __name__ == '__main__':
    main()