#!/usr/bin/env python3
"""
Unified Amyloid Database Merger

Merges data from multiple amyloid databases into a single unified table
and generates multi-FASTA file with amyloidogenic sequences.

Supported databases:
- AmyloidExplorer (cryo-EM/NMR structures)
- AmyloidAtlas (structural atlas)
- AmyLoad (experimental peptide validation)
- Amylobase (kinetic data)
- AmyPro (functional amyloids/prions)
- AmyloGraph (cross-seeding interactions)
- WALTZ-DB (hexapeptide predictions)

Author: Claude (Anthropic)
Date: 2025
"""

import csv
import re
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict


@dataclass
class AmyloidRegion:
    """Represents a single amyloidogenic region."""
    # Core identifiers
    record_id: str = ""
    source_db: str = ""
    
    # Protein info
    protein_name: str = ""
    uniprot_id: str = ""
    organism: str = ""
    
    # Region coordinates
    region_start: Optional[int] = None
    region_end: Optional[int] = None
    
    # Sequence data
    sequence: str = ""
    full_protein_sequence: str = ""
    
    # Classification
    is_amyloid: bool = True
    category: str = ""  # pathogenic, functional, prion, etc.
    
    # Structural data
    pdb_id: str = ""
    emdb_id: str = ""
    method: str = ""
    resolution: str = ""
    
    # Experimental evidence
    disease: str = ""
    tissue: str = ""
    mutation: str = ""
    
    # References
    doi: str = ""
    pmid: str = ""
    reference: str = ""
    
    # Additional metadata
    notes: str = ""


def clean_value(val) -> str:
    """Clean and normalize a value."""
    if val is None:
        return ""
    val = str(val).strip()
    if val.lower() in ['na', 'n/a', 'none', 'null', '-', '']:
        return ""
    return val


def parse_range(range_str: str) -> List[Tuple[Optional[int], Optional[int]]]:
    """
    Parse range string into list of (start, end) tuples.
    Handles formats: "1-100", "1-100,200-300", "1-10 20-30", etc.
    """
    if not range_str:
        return [(None, None)]
    
    ranges = []
    # Split by comma or space
    parts = re.split(r'[,\s]+', range_str.strip())
    
    for part in parts:
        part = part.strip()
        if not part:
            continue
            
        # Match patterns like "1-100", "1-100 (99)", etc.
        match = re.match(r'(\d+)\s*[-–]\s*(\d+)', part)
        if match:
            start, end = int(match.group(1)), int(match.group(2))
            if start <= end:  # Valid range
                ranges.append((start, end))
        else:
            # Try single number
            match = re.match(r'^(\d+)$', part)
            if match:
                pos = int(match.group(1))
                ranges.append((pos, pos))
    
    return ranges if ranges else [(None, None)]


def extract_uniprot(uniprot_str: str) -> str:
    """Extract clean UniProt ID from various formats."""
    if not uniprot_str:
        return ""
    
    # Handle formats like "P12345 20-100" or "P12345"
    match = re.match(r'^([A-Z][A-Z0-9]{5,9})', uniprot_str.strip())
    if match:
        return match.group(1)
    return clean_value(uniprot_str)


def parse_amyloid_explorer(filepath: str) -> List[AmyloidRegion]:
    """Parse AmyloidExplorer TSV file."""
    regions = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            try:
                start = int(row.get('resStart', 0)) if row.get('resStart') else None
                end = int(row.get('resEnd', 0)) if row.get('resEnd') else None
            except (ValueError, TypeError):
                start, end = None, None
            
            region = AmyloidRegion(
                record_id=clean_value(row.get('Name', '')),
                source_db="AmyloidExplorer",
                protein_name=clean_value(row.get('Protein', '')),
                uniprot_id=extract_uniprot(row.get('UniprotId', '')),
                region_start=start,
                region_end=end,
                pdb_id=clean_value(row.get('Name', '')),  # PDB code in Name field
                emdb_id=clean_value(row.get('EMDB', '')),
                method=clean_value(row.get('Method', '')),
                resolution=clean_value(row.get('Resolution', '')),
                disease=clean_value(row.get('Disease', '')),
                tissue=clean_value(row.get('Tissue', '')),
                mutation=clean_value(row.get('Mutant', '')),
                doi=clean_value(row.get('DOI', '')),
                pmid=clean_value(row.get('PMID', '')),
                category="pathogenic" if row.get('Patient') == 'True' else "in_vitro",
                is_amyloid=True
            )
            regions.append(region)
    
    return regions


def parse_amyloid_atlas(filepath: str) -> List[AmyloidRegion]:
    """Parse AmyloidAtlas TSV file."""
    regions = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            residues = clean_value(row.get('Residues_Ordered', ''))
            range_list = parse_range(residues)
            
            for start, end in range_list:
                region = AmyloidRegion(
                    record_id=clean_value(row.get('PDB_ID', '')),
                    source_db="AmyloidAtlas",
                    protein_name=clean_value(row.get('Protein', '')),
                    region_start=start,
                    region_end=end,
                    pdb_id=clean_value(row.get('PDB_ID', '')),
                    method=clean_value(row.get('Method', '')),
                    resolution=clean_value(row.get('Resolution_A', '')),
                    doi=clean_value(row.get('DOI_URL', '')).replace('https://doi.org/', ''),
                    reference=clean_value(row.get('Reference', '')),
                    notes=clean_value(row.get('Fibril_Origins', '')),
                    is_amyloid=True
                )
                regions.append(region)
    
    return regions


def parse_amyload(filepath: str) -> List[AmyloidRegion]:
    """Parse AmyLoad CSV file (no header)."""
    regions = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        
        for row in reader:
            if len(row) < 5:
                continue
            
            amy_id, protein, fragment, sequence, is_amyloid = row[:5]
            
            region = AmyloidRegion(
                record_id=clean_value(amy_id),
                source_db="AmyLoad",
                protein_name=clean_value(protein),
                sequence=clean_value(sequence),
                is_amyloid=clean_value(is_amyloid).lower() == 'yes',
                notes=f"Fragment: {clean_value(fragment)}",
                category="experimental"
            )
            regions.append(region)
    
    return regions


def parse_amylobase(filepath: str) -> List[AmyloidRegion]:
    """Parse Amylobase tab-delimited file."""
    regions = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            try:
                start = int(row.get('Start', 0)) if row.get('Start') else None
                end = int(row.get('Stop', 0)) if row.get('Stop') else None
            except (ValueError, TypeError):
                start, end = None, None
            
            mutation = ""
            if row.get('Mutation') and row.get('Mutation') != 'na':
                mutation = clean_value(row.get('Mutation', ''))
            
            region = AmyloidRegion(
                record_id=f"AB_{extract_uniprot(row.get('Protein ID', ''))}",
                source_db="Amylobase",
                protein_name=clean_value(row.get('Protein name', '')),
                uniprot_id=extract_uniprot(row.get('Protein ID', '')),
                organism=clean_value(row.get('Organism', '')),
                region_start=start,
                region_end=end,
                mutation=mutation,
                category="kinetic_data",
                is_amyloid=True,
                notes=f"pH={row.get('pH', '')}, T={row.get('Temperature (Kelvin Degrees)', '')}K"
            )
            regions.append(region)
    
    return regions


def parse_amypro(filepath: str) -> List[AmyloidRegion]:
    """Parse AmyPro tab-delimited file."""
    regions = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            regions_str = clean_value(row.get('regions', ''))
            range_list = parse_range(regions_str)
            
            # Determine category
            category = clean_value(row.get('category', ''))
            is_prion = row.get('prion', '').lower() == 'true'
            if is_prion:
                category = f"prion/{category}" if category else "prion"
            
            for start, end in range_list:
                region = AmyloidRegion(
                    record_id=clean_value(row.get('ID', '')).replace('#', ''),
                    source_db="AmyPro",
                    protein_name=clean_value(row.get('protein', '')),
                    uniprot_id=extract_uniprot(row.get('uniprot', '')),
                    organism=clean_value(row.get('species', '')),
                    region_start=start,
                    region_end=end,
                    sequence=clean_value(row.get('sequence', '')),
                    full_protein_sequence=clean_value(row.get('sequence', '')),
                    pdb_id=clean_value(row.get('pdb', '')),
                    pmid=clean_value(row.get('pubmed', '')),
                    category=category,
                    mutation=clean_value(row.get('mutations', '')),
                    is_amyloid=True
                )
                regions.append(region)
    
    return regions


def parse_amylograph(filepath: str) -> List[AmyloidRegion]:
    """Parse AmyloGraph CSV file (cross-seeding data)."""
    regions = []
    seen_sequences = set()
    
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            # Extract unique sequences from interactor and interactee
            for prefix in ['interactor', 'interactee']:
                name = clean_value(row.get(f'{prefix}_name', ''))
                seq = clean_value(row.get(f'{prefix}_sequence', ''))
                
                if seq and seq not in seen_sequences:
                    seen_sequences.add(seq)
                    
                    region = AmyloidRegion(
                        record_id=f"AG_{name}_{len(seen_sequences)}",
                        source_db="AmyloGraph",
                        protein_name=name,
                        sequence=seq,
                        doi=clean_value(row.get('doi', '')),
                        category="cross_seeding",
                        is_amyloid=True,
                        notes=f"AGID: {row.get('AGID', '')}"
                    )
                    regions.append(region)
    
    return regions


def parse_waltzdb(filepath: str) -> List[AmyloidRegion]:
    """Parse WALTZ-DB TSV file."""
    regions = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            # Handle both 'Sequence' and 'sequence' column names
            seq = clean_value(row.get('Sequence', '') or row.get('sequence', ''))
            # Handle both 'Classification' and 'classification'
            classification = clean_value(row.get('Classification', '') or row.get('classification', ''))
            
            if not seq:  # Skip empty sequences
                continue
            
            # Extract scores
            waltz = clean_value(row.get('WALTZ', '') or row.get('waltz', ''))
            tango = clean_value(row.get('TANGO', '') or row.get('tango', ''))
            pasta_p = clean_value(row.get('PASTA parallel', '') or row.get('pasta_parallel', ''))
            pasta_ap = clean_value(row.get('PASTA antiparallel', '') or row.get('pasta_antiparallel', ''))
            
            # Extract UniProt info
            # Note: WALTZ-DB has non-standard naming - 'UniProt AC' contains entry name, 'UniProt ID' contains accession
            uniprot_accession = clean_value(row.get('UniProt ID', '') or row.get('uniprot id', ''))
            uniprot_entry_name = clean_value(row.get('UniProt AC', '') or row.get('uniprot ac', ''))
            
            # Position in parent protein
            position = clean_value(row.get('Position', '') or row.get('position', ''))
            
            region = AmyloidRegion(
                record_id=f"WDB_{seq}",
                source_db="WALTZ-DB",
                sequence=seq,
                uniprot_id=extract_uniprot(uniprot_accession) if uniprot_accession and uniprot_accession != 'N.A.' else '',
                protein_name=uniprot_entry_name if uniprot_entry_name and uniprot_entry_name != 'N.A.' else '',
                region_start=int(position) if position and position.isdigit() else None,
                region_end=int(position) + len(seq) - 1 if position and position.isdigit() else None,
                is_amyloid=classification.lower() == 'amyloid',
                category="hexapeptide",
                pdb_id=clean_value(row.get('Structures', '')) if row.get('Structures', '') != 'N.A.' else '',
                notes=f"WALTZ={waltz}, TANGO={tango}, PASTA_P={pasta_p}, PASTA_AP={pasta_ap}"
            )
            regions.append(region)
    
    return regions


def merge_databases(input_files: Dict[str, str]) -> List[AmyloidRegion]:
    """
    Merge all database files into unified list.
    
    Args:
        input_files: Dict mapping database name to filepath
    
    Returns:
        List of AmyloidRegion objects
    """
    all_regions = []
    
    parsers = {
        'amyloid_explorer': parse_amyloid_explorer,
        'amyloid_atlas': parse_amyloid_atlas,
        'amyload': parse_amyload,
        'amylobase': parse_amylobase,
        'amypro': parse_amypro,
        'amylograph': parse_amylograph,
        'waltzdb': parse_waltzdb,
    }
    
    for db_name, filepath in input_files.items():
        if not filepath or not Path(filepath).exists():
            print(f"  Skipping {db_name}: file not found")
            continue
        
        parser = parsers.get(db_name.lower().replace('-', '_').replace(' ', '_'))
        if parser:
            try:
                regions = parser(filepath)
                all_regions.extend(regions)
                print(f"  {db_name}: {len(regions)} regions loaded")
            except Exception as e:
                print(f"  Error parsing {db_name}: {e}")
        else:
            print(f"  Unknown database: {db_name}")
    
    return all_regions


def write_unified_table(regions: List[AmyloidRegion], output_path: str, 
                        amyloid_only: bool = False):
    """Write unified table to TSV file."""
    
    if amyloid_only:
        regions = [r for r in regions if r.is_amyloid]
    
    headers = [
        'record_id', 'source_db', 'protein_name', 'uniprot_id', 'organism',
        'region_start', 'region_end', 'sequence', 'is_amyloid', 'category',
        'pdb_id', 'emdb_id', 'method', 'resolution',
        'disease', 'tissue', 'mutation',
        'doi', 'pmid', 'reference', 'notes'
    ]
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        
        for region in regions:
            row = [
                region.record_id,
                region.source_db,
                region.protein_name,
                region.uniprot_id,
                region.organism,
                region.region_start if region.region_start else '',
                region.region_end if region.region_end else '',
                region.sequence,
                'yes' if region.is_amyloid else 'no',
                region.category,
                region.pdb_id,
                region.emdb_id,
                region.method,
                region.resolution,
                region.disease,
                region.tissue,
                region.mutation,
                region.doi,
                region.pmid,
                region.reference,
                region.notes
            ]
            writer.writerow(row)
    
    print(f"Wrote {len(regions)} records to {output_path}")


def generate_fasta(regions: List[AmyloidRegion], output_path: str,
                   amyloid_only: bool = True, deduplicate: bool = True,
                   min_length: int = 4):
    """
    Generate multi-FASTA file from regions.
    
    Args:
        regions: List of AmyloidRegion objects
        output_path: Output FASTA file path
        amyloid_only: Only include amyloid-positive sequences
        deduplicate: Remove duplicate sequences
        min_length: Minimum sequence length
    """
    if amyloid_only:
        regions = [r for r in regions if r.is_amyloid]
    
    # Collect sequences with metadata
    sequences = {}  # seq -> (header, region)
    
    for region in regions:
        seq = region.sequence.upper().strip()
        if not seq or len(seq) < min_length:
            continue
        
        # Clean sequence (remove non-amino acid characters)
        seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq)
        if len(seq) < min_length:
            continue
        
        # Build header
        parts = [region.source_db]
        if region.record_id:
            parts.append(region.record_id)
        if region.protein_name:
            parts.append(region.protein_name.replace(' ', '_'))
        if region.uniprot_id:
            parts.append(f"UniProt:{region.uniprot_id}")
        if region.region_start and region.region_end:
            parts.append(f"region:{region.region_start}-{region.region_end}")
        if region.category:
            parts.append(f"category:{region.category}")
        
        header = "|".join(parts)
        
        if deduplicate:
            if seq not in sequences:
                sequences[seq] = (header, region)
        else:
            # Add suffix for duplicates
            base_seq = seq
            counter = 1
            while seq in sequences:
                seq = f"{base_seq}_{counter}"
                counter += 1
            sequences[seq] = (header, region)
    
    # Write FASTA
    with open(output_path, 'w', encoding='utf-8') as f:
        for seq, (header, region) in sequences.items():
            f.write(f">{header}\n")
            # Wrap sequence at 60 characters
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")
    
    print(f"Wrote {len(sequences)} sequences to {output_path}")


def print_statistics(regions: List[AmyloidRegion]):
    """Print summary statistics."""
    print("\n" + "="*60)
    print("UNIFIED AMYLOID DATABASE STATISTICS")
    print("="*60)
    
    # By source
    by_source = defaultdict(int)
    for r in regions:
        by_source[r.source_db] += 1
    
    print("\nRecords by source database:")
    for db, count in sorted(by_source.items()):
        print(f"  {db}: {count}")
    
    # Amyloid vs non-amyloid
    amyloid = sum(1 for r in regions if r.is_amyloid)
    non_amyloid = len(regions) - amyloid
    print(f"\nClassification:")
    print(f"  Amyloid-positive: {amyloid}")
    print(f"  Amyloid-negative: {non_amyloid}")
    
    # With coordinates
    with_coords = sum(1 for r in regions if r.region_start and r.region_end)
    print(f"\nWith region coordinates: {with_coords}")
    
    # With sequences
    with_seq = sum(1 for r in regions if r.sequence)
    print(f"With peptide sequences: {with_seq}")
    
    # With UniProt
    with_uniprot = sum(1 for r in regions if r.uniprot_id)
    print(f"With UniProt ID: {with_uniprot}")
    
    # With PDB
    with_pdb = sum(1 for r in regions if r.pdb_id)
    print(f"With PDB ID: {with_pdb}")
    
    # Unique proteins
    unique_proteins = len(set(r.protein_name.lower() for r in regions if r.protein_name))
    print(f"Unique protein names: {unique_proteins}")
    
    # Unique UniProt IDs
    unique_uniprot = len(set(r.uniprot_id for r in regions if r.uniprot_id))
    print(f"Unique UniProt IDs: {unique_uniprot}")
    
    # Categories
    by_category = defaultdict(int)
    for r in regions:
        if r.category:
            by_category[r.category] += 1
    
    if by_category:
        print("\nBy category:")
        for cat, count in sorted(by_category.items(), key=lambda x: -x[1]):
            print(f"  {cat}: {count}")
    
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Merge multiple amyloid databases into unified format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python merge_amyloid_databases.py \\
      --amyloid-explorer AmyloidExplorer.tsv \\
      --amyloid-atlas AmyloidAtlas.tsv \\
      --amyload AmyLoad.csv \\
      --amylobase amylobase.txt \\
      --amypro amypro.txt \\
      --amylograph amylograph.txt \\
      --waltzdb waltzdb_peptides.tsv \\
      -o unified_amyloid_db.tsv \\
      --fasta amyloid_sequences.fasta
        """
    )
    
    # Input files
    parser.add_argument('--amyloid-explorer', help='AmyloidExplorer TSV file')
    parser.add_argument('--amyloid-atlas', help='AmyloidAtlas TSV file')
    parser.add_argument('--amyload', help='AmyLoad CSV file')
    parser.add_argument('--amylobase', help='Amylobase TXT file')
    parser.add_argument('--amypro', help='AmyPro TXT file')
    parser.add_argument('--amylograph', help='AmyloGraph TXT file')
    parser.add_argument('--waltzdb', help='WALTZ-DB TSV file')
    
    # Output options
    parser.add_argument('-o', '--output', default='unified_amyloid_db.tsv',
                        help='Output TSV file (default: unified_amyloid_db.tsv)')
    parser.add_argument('--fasta', help='Output FASTA file for sequences')
    parser.add_argument('--amyloid-only', action='store_true',
                        help='Only include amyloid-positive entries')
    parser.add_argument('--no-deduplicate', action='store_true',
                        help='Do not deduplicate sequences in FASTA')
    parser.add_argument('--min-length', type=int, default=4,
                        help='Minimum sequence length for FASTA (default: 4)')
    parser.add_argument('--stats', action='store_true',
                        help='Print statistics')
    
    args = parser.parse_args()
    
    # Collect input files
    input_files = {
        'amyloid_explorer': args.amyloid_explorer,
        'amyloid_atlas': args.amyloid_atlas,
        'amyload': args.amyload,
        'amylobase': args.amylobase,
        'amypro': args.amypro,
        'amylograph': args.amylograph,
        'waltzdb': args.waltzdb,
    }
    
    # Filter to only provided files
    input_files = {k: v for k, v in input_files.items() if v}
    
    if not input_files:
        parser.error("At least one input database file is required")
    
    print("Loading databases...")
    regions = merge_databases(input_files)
    
    if not regions:
        print("No data loaded!")
        return
    
    print(f"\nTotal regions loaded: {len(regions)}")
    
    # Write unified table
    write_unified_table(regions, args.output, amyloid_only=args.amyloid_only)
    
    # Generate FASTA if requested
    if args.fasta:
        generate_fasta(
            regions, args.fasta,
            amyloid_only=args.amyloid_only,
            deduplicate=not args.no_deduplicate,
            min_length=args.min_length
        )
    
    # Print statistics
    if args.stats:
        print_statistics(regions)


if __name__ == '__main__':
    main()