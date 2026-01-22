#!/usr/bin/env python3
"""
Amyloid Atlas HTML Parser - Fixed for malformed HTML
=====================================================
The original HTML has missing <tr> tags - uses <!--Row PDB_ID--> comments instead.

This parser splits by Row comments and extracts <td> cells for each entry.

Usage:
    python parse_atlas_fixed.py -i AmyloidAtlas.html -o amyloid_atlas.tsv
"""

import argparse
import csv
import re
import sys
from pathlib import Path

try:
    from bs4 import BeautifulSoup
    HAS_BS4 = True
except ImportError:
    HAS_BS4 = False


def clean_text(text):
    """Normalize whitespace, decode HTML entities, clean text."""
    if not text:
        return ''
    # Decode common HTML entities
    text = text.replace('&alpha;', 'α')
    text = text.replace('&beta;', 'β')
    text = text.replace('&gamma;', 'γ')
    text = text.replace('&kappa;', 'κ')
    text = text.replace('&lambda;', 'λ')
    text = text.replace('&tau;', 'τ')
    text = text.replace('&#x212B;', 'Å')
    # Normalize whitespace
    text = re.sub(r'\s+', ' ', text)
    return text.strip()


def extract_link_from_html(html_fragment):
    """Extract text and URL from HTML fragment containing an anchor."""
    if not HAS_BS4:
        # Regex fallback
        link_match = re.search(r'<a\s+href=["\']([^"\']+)["\'][^>]*>([^<]+)</a>', html_fragment, re.I)
        if link_match:
            return clean_text(link_match.group(2)), link_match.group(1)
        # No link, extract text
        text = re.sub(r'<[^>]+>', '', html_fragment)
        return clean_text(text), ''
    
    soup = BeautifulSoup(html_fragment, 'html.parser')
    link = soup.find('a')
    if link:
        return clean_text(link.get_text()), link.get('href', '')
    return clean_text(soup.get_text()), ''


def parse_html_by_comments(html_content):
    """
    Parse HTML by splitting on <!--Row PDB_ID--> comments.
    
    Each entry is marked by a comment like <!--Row 6cu8--> followed by 9 <td> cells.
    """
    entries = []
    seen_pdb_ids = set()
    
    # Split by Row comments
    # Pattern: <!--Row PDBID--> where PDBID is alphanumeric
    row_pattern = re.compile(r'<!--\s*Row\s+(\w+)\s*-->', re.IGNORECASE)
    
    # Find all Row comments and their positions
    matches = list(row_pattern.finditer(html_content))
    
    if not matches:
        print("Warning: No <!--Row--> comments found in HTML", file=sys.stderr)
        return entries
    
    print(f"Found {len(matches)} Row comments")
    
    for i, match in enumerate(matches):
        pdb_id_from_comment = match.group(1).lower()
        
        # Get content from this match to next match (or end of file)
        start_pos = match.end()
        end_pos = matches[i + 1].start() if i + 1 < len(matches) else len(html_content)
        
        row_content = html_content[start_pos:end_pos]
        
        # Extract all <td>...</td> cells from this segment
        # Use non-greedy matching
        td_pattern = re.compile(r'<td[^>]*>(.*?)</td>', re.DOTALL | re.IGNORECASE)
        cells = td_pattern.findall(row_content)
        
        if len(cells) < 9:
            # Try alternative: some cells may not be properly closed
            # Look for <td...> and capture until next <td or </tr
            td_pattern2 = re.compile(r'<td[^>]*>(.*?)(?=<td|</tr|$)', re.DOTALL | re.IGNORECASE)
            cells = td_pattern2.findall(row_content)
        
        if len(cells) < 9:
            # Try to recover partial data for malformed rows
            print(f"Warning: Row {pdb_id_from_comment} has only {len(cells)} cells, attempting recovery", file=sys.stderr)
            
            # Try to extract at least basic info
            if len(cells) >= 5:
                protein = clean_text(re.sub(r'<[^>]+>', '', cells[0]))
                fibril_origins = clean_text(re.sub(r'<[^>]+>', '', cells[1]))
                residues = clean_text(re.sub(r'<[^>]+>', '', cells[2]))
                method = clean_text(re.sub(r'<[^>]+>', '', cells[3]))
                resolution = clean_text(re.sub(r'<[^>]+>', '', cells[4]))
                
                # Use PDB ID from comment
                pdb_id = pdb_id_from_comment
                pdb_url = f"https://www.rcsb.org/structure/{pdb_id}"
                
                # Try to find reference in remaining content
                ref_match = re.search(r'<a\s+href=["\']([^"\']+doi[^"\']+)["\'][^>]*>([^<]+)</a>', row_content, re.I)
                if ref_match:
                    doi_url = ref_match.group(1)
                    reference = clean_text(ref_match.group(2))
                else:
                    reference = ""
                    doi_url = ""
                
                pdb_key = pdb_id.lower()
                if pdb_key not in seen_pdb_ids:
                    seen_pdb_ids.add(pdb_key)
                    entries.append({
                        'Protein': protein,
                        'Fibril_Origins': fibril_origins,
                        'Residues_Ordered': residues,
                        'Method': method,
                        'Resolution_A': resolution,
                        'PDB_ID': pdb_id,
                        'PDB_URL': pdb_url,
                        'Reference': reference,
                        'DOI_URL': doi_url
                    })
                    print(f"  Recovered: {protein[:30]} | {pdb_id}", file=sys.stderr)
                continue
            else:
                print(f"  Cannot recover - too few cells", file=sys.stderr)
                continue
        
        # Extract data from cells
        # 0: Protein, 1: Fibril Origins, 2: Residues, 3: Method, 4: Resolution
        # 5: Polarity Map (skip), 6: Energy Map (skip)
        # 7: PDB ID, 8: Reference
        
        protein = clean_text(re.sub(r'<[^>]+>', '', cells[0]))
        fibril_origins = clean_text(re.sub(r'<[^>]+>', '', cells[1]))
        residues = clean_text(re.sub(r'<[^>]+>', '', cells[2]))
        method = clean_text(re.sub(r'<[^>]+>', '', cells[3]))
        resolution = clean_text(re.sub(r'<[^>]+>', '', cells[4]))
        
        # Skip cells 5 and 6 (Polarity and Energy maps)
        
        pdb_id, pdb_url = extract_link_from_html(cells[7])
        reference, doi_url = extract_link_from_html(cells[8])
        
        # Use PDB ID from cell if available, otherwise from comment
        if not pdb_id:
            pdb_id = pdb_id_from_comment
        
        # Skip duplicates
        pdb_key = pdb_id.lower()
        if pdb_key in seen_pdb_ids:
            continue
        seen_pdb_ids.add(pdb_key)
        
        entries.append({
            'Protein': protein,
            'Fibril_Origins': fibril_origins,
            'Residues_Ordered': residues,
            'Method': method,
            'Resolution_A': resolution,
            'PDB_ID': pdb_id,
            'PDB_URL': pdb_url,
            'Reference': reference,
            'DOI_URL': doi_url
        })
    
    return entries


def write_tsv(entries, output_path):
    """Write entries to TSV file."""
    fieldnames = [
        'Protein', 'Fibril_Origins', 'Residues_Ordered', 'Method',
        'Resolution_A', 'PDB_ID', 'PDB_URL', 'Reference', 'DOI_URL'
    ]
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(entries)
    
    return len(entries)


def print_summary(entries):
    """Print summary statistics."""
    if not entries:
        return
    
    methods = {}
    for e in entries:
        m = e['Method']
        methods[m] = methods.get(m, 0) + 1
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total entries: {len(entries)}")
    
    print(f"\nEntries by method:")
    for method, count in sorted(methods.items(), key=lambda x: -x[1]):
        print(f"  {method}: {count}")
    
    print(f"\nFirst 5 entries:")
    for e in entries[:5]:
        print(f"  {e['Protein'][:35]:35} | {e['PDB_ID']:6} | {e['Method']:8}")
    
    print(f"\nLast 5 entries:")
    for e in entries[-5:]:
        print(f"  {e['Protein'][:35]:35} | {e['PDB_ID']:6} | {e['Method']:8}")


def main():
    parser = argparse.ArgumentParser(
        description='Parse Amyloid Atlas HTML to TSV (handles malformed HTML)'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input HTML file')
    parser.add_argument('-o', '--output', default='amyloid_atlas.tsv',
                        help='Output TSV file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print summary')
    
    args = parser.parse_args()
    
    # Read input
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: File not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    content = input_path.read_text(encoding='utf-8', errors='replace')
    print(f"Read {len(content):,} characters from {args.input}")
    
    # Parse
    entries = parse_html_by_comments(content)
    
    if not entries:
        print("Error: No entries extracted!", file=sys.stderr)
        sys.exit(1)
    
    # Write output
    count = write_tsv(entries, args.output)
    print(f"\nExtracted {count} entries to {args.output}")
    
    if args.verbose:
        print_summary(entries)


if __name__ == '__main__':
    main()