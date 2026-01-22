#!/usr/bin/env python3
"""
WALTZ-DB 2.0 Scraper
====================
Downloads peptide data from WALTZ-DB 2.0 (http://waltzdb.switchlab.org/)

WALTZ-DB is a Drupal-based database containing experimentally characterized 
amyloidogenic hexapeptides with:
- ~1400+ peptides (515 amyloid-forming, 901 non-amyloidogenic)
- Experimental data: TEM images, FTIR spectra, ThT binding
- Predictions: WALTZ, TANGO, PASTA scores
- Structural models of fibril cores

The site has built-in CSV/Excel export buttons, but this script can also
parse the HTML table directly.

Requirements:
    pip install requests beautifulsoup4 pandas lxml

Usage:
    python waltzdb_scraper.py
    python waltzdb_scraper.py --output waltzdb_peptides.tsv --format tsv
    python waltzdb_scraper.py --amyloid-only  # Only amyloid-forming peptides

Author: Generated for Xenia's amyloid research  
Date: 2026-01-19
"""

import argparse
import json
import csv
import re
import sys
import time
from pathlib import Path
from typing import List, Dict, Any, Optional
from urllib.parse import urljoin, urlencode

try:
    import requests
    from bs4 import BeautifulSoup
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False
    print("Missing dependencies. Install with:")
    print("  pip install requests beautifulsoup4 lxml")

BASE_URL = "http://waltzdb.switchlab.org"
SEQUENCES_URL = f"{BASE_URL}/sequences"


def get_session() -> requests.Session:
    """Create a session with appropriate headers."""
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Connection': 'keep-alive',
    })
    return session


def try_direct_download(session: requests.Session, verbose: bool = False) -> Optional[str]:
    """
    Try to find and use the built-in CSV/Excel download functionality.
    WALTZ-DB has export buttons at the bottom of the table.
    """
    # Common Drupal/Views export paths
    export_paths = [
        '/sequences/export/csv',
        '/sequences?_format=csv',
        '/export/sequences/csv',
        '/sequences/csv',
        '/api/sequences',
        '/sequences?export=csv',
    ]
    
    for path in export_paths:
        url = BASE_URL + path
        if verbose:
            print(f"  Trying: {url}")
        
        try:
            resp = session.get(url, timeout=15)
            content_type = resp.headers.get('Content-Type', '')
            
            # Check if we got CSV/Excel data
            if resp.status_code == 200:
                if 'csv' in content_type or 'excel' in content_type or 'spreadsheet' in content_type:
                    print(f"Found direct download at: {url}")
                    return resp.text
                
                # Check content looks like CSV
                if resp.text.startswith('Sequence,') or resp.text.startswith('"Sequence"'):
                    print(f"Found CSV data at: {url}")
                    return resp.text
                    
        except Exception as e:
            if verbose:
                print(f"    Error: {e}")
    
    return None


def parse_sequences_page(session: requests.Session, page: int = 0, 
                         verbose: bool = False) -> tuple[List[Dict], bool]:
    """
    Parse a single page of the sequences table.
    Returns (entries, has_more_pages)
    """
    # WALTZ-DB uses Drupal Views pagination
    params = {}
    if page > 0:
        params['page'] = page
    
    url = SEQUENCES_URL
    if params:
        url = f"{SEQUENCES_URL}?{urlencode(params)}"
    
    if verbose:
        print(f"  Fetching page {page}: {url}")
    
    try:
        resp = session.get(url, timeout=30)
        resp.raise_for_status()
    except Exception as e:
        print(f"Error fetching page {page}: {e}")
        return [], False
    
    soup = BeautifulSoup(resp.text, 'lxml')
    entries = []
    
    # Find the main data table
    # Drupal Views tables often have class 'views-table' or similar
    table = soup.find('table', class_=re.compile(r'views|data|sequences'))
    if not table:
        table = soup.find('table')
    
    if not table:
        if verbose:
            print("    No table found on page")
        return [], False
    
    # Extract headers
    headers = []
    header_row = table.find('thead')
    if header_row:
        headers = [th.get_text(strip=True) for th in header_row.find_all(['th', 'td'])]
    else:
        # Try first row
        first_row = table.find('tr')
        if first_row:
            headers = [cell.get_text(strip=True) for cell in first_row.find_all(['th', 'td'])]
    
    # Normalize headers
    headers = [normalize_header(h) for h in headers]
    
    if verbose and page == 0:
        print(f"    Headers: {headers}")
    
    # Extract data rows
    rows = table.find_all('tr')
    for row in rows:
        cells = row.find_all('td')
        if not cells:
            continue
        
        entry = {}
        for idx, cell in enumerate(cells):
            key = headers[idx] if idx < len(headers) else f'col_{idx}'
            
            # Get text content
            text = cell.get_text(strip=True)
            entry[key] = text
            
            # Extract links (for peptide detail pages, UniProt, PDB)
            links = cell.find_all('a')
            for link in links:
                href = link.get('href', '')
                if 'uniprot' in href.lower():
                    entry['uniprot_url'] = href
                    # Extract UniProt ID
                    match = re.search(r'([A-Z][A-Z0-9]{5})', href)
                    if match:
                        entry['uniprot_id'] = match.group(1)
                elif 'pdb' in href.lower() or 'rcsb' in href.lower():
                    entry['pdb_url'] = href
                elif href.startswith('/sequence/'):
                    entry['detail_url'] = urljoin(BASE_URL, href)
            
            # Check for images (TEM available indicator)
            if cell.find('img'):
                entry[f'{key}_has_image'] = True
        
        if entry and entry.get('sequence', ''):
            entries.append(entry)
    
    # Check for pagination (next page)
    has_more = False
    pager = soup.find('ul', class_=re.compile(r'pager|pagination'))
    if pager:
        next_link = pager.find('a', text=re.compile(r'next|›|»', re.I))
        if not next_link:
            next_link = pager.find('li', class_=re.compile(r'next'))
        has_more = next_link is not None
    
    # Also check for "page=N" links
    if not has_more:
        next_page_link = soup.find('a', href=re.compile(f'page={page+1}'))
        has_more = next_page_link is not None
    
    if verbose:
        print(f"    Found {len(entries)} entries, has_more={has_more}")
    
    return entries, has_more


def normalize_header(header: str) -> str:
    """Normalize column header to a clean key name."""
    # Common transformations
    header = header.lower().strip()
    header = re.sub(r'[^\w\s\(\)]', '', header)  # Keep parentheses for PASTA
    header = re.sub(r'\s+', '_', header)
    
    # Known mappings based on actual WALTZ-DB HTML structure
    mappings = {
        'sequence': 'sequence',
        'classification': 'classification',
        'source': 'source',
        'subset': 'subset',
        'waltz': 'waltz',
        'tango': 'tango',
        'pasta_(p)': 'pasta_parallel',
        'pasta_(ap)': 'pasta_antiparallel',
        'hp': 'hydrophobicity',
        'cf_h': 'chou_fasman_helix',
        'cf_s': 'chou_fasman_strand',
        # Legacy mappings for compatibility
        'peptide_sequence': 'sequence',
        'seq': 'sequence',
        'amyloid': 'classification',
        'amyloid_forming': 'classification',
        'morphology': 'classification',
        'morphology_decision': 'classification',
        'source_protein': 'source',
        'protein': 'source',
        'waltz_score': 'waltz',
        'tango_score': 'tango',
        'pasta_score': 'pasta_parallel',
        'hydrophobicity': 'hydrophobicity',
    }
    
    return mappings.get(header, header)


def fetch_peptide_details(session: requests.Session, url: str, 
                          verbose: bool = False) -> Dict:
    """
    Fetch additional details from a peptide's detail page.
    This includes experimental data, predictions, structural info.
    """
    try:
        resp = session.get(url, timeout=15)
        resp.raise_for_status()
    except Exception as e:
        if verbose:
            print(f"    Error fetching details: {e}")
        return {}
    
    soup = BeautifulSoup(resp.text, 'lxml')
    details = {}
    
    # Look for data tables on the detail page
    tables = soup.find_all('table')
    
    for table in tables:
        rows = table.find_all('tr')
        for row in rows:
            cells = row.find_all(['th', 'td'])
            if len(cells) >= 2:
                key = cells[0].get_text(strip=True)
                value = cells[1].get_text(strip=True)
                if key and value:
                    details[normalize_header(key)] = value
    
    # Look for specific fields in divs (common in Drupal)
    for field_div in soup.find_all('div', class_=re.compile(r'field')):
        label = field_div.find(class_=re.compile(r'label|title'))
        value = field_div.find(class_=re.compile(r'value|content|item'))
        if label and value:
            key = label.get_text(strip=True)
            val = value.get_text(strip=True)
            if key and val:
                details[normalize_header(key)] = val
    
    return details


def scrape_all_peptides(session: requests.Session, verbose: bool = False,
                        fetch_details: bool = False, max_pages: int = 100) -> List[Dict]:
    """
    Scrape all peptides from WALTZ-DB.
    """
    all_entries = []
    page = 0
    
    print("Scraping WALTZ-DB peptide sequences...")
    
    while page < max_pages:
        entries, has_more = parse_sequences_page(session, page, verbose)
        
        if not entries:
            if page == 0:
                print("No entries found on first page. The site structure may have changed.")
            break
        
        all_entries.extend(entries)
        print(f"  Page {page}: {len(entries)} entries (total: {len(all_entries)})")
        
        if not has_more:
            break
        
        page += 1
        time.sleep(0.5)  # Be polite
    
    # Optionally fetch details for each peptide
    if fetch_details and all_entries:
        print(f"\nFetching details for {len(all_entries)} peptides...")
        for i, entry in enumerate(all_entries):
            if 'detail_url' in entry:
                if verbose or (i % 50 == 0):
                    print(f"  {i+1}/{len(all_entries)}: {entry.get('sequence', '?')}")
                
                details = fetch_peptide_details(session, entry['detail_url'], verbose)
                entry.update(details)
                time.sleep(0.3)
    
    return all_entries


def filter_entries(entries: List[Dict], amyloid_only: bool = False,
                   non_amyloid_only: bool = False,
                   source_filter: str = None) -> List[Dict]:
    """Filter entries based on criteria."""
    filtered = entries
    
    if amyloid_only:
        filtered = [e for e in filtered if is_amyloid_positive(e)]
        print(f"Filtered to {len(filtered)} amyloid-positive entries")
    
    if non_amyloid_only:
        filtered = [e for e in filtered if not is_amyloid_positive(e)]
        print(f"Filtered to {len(filtered)} non-amyloidogenic entries")
    
    if source_filter:
        filtered = [e for e in filtered 
                   if source_filter.lower() in e.get('source', '').lower()]
        print(f"Filtered to {len(filtered)} entries from '{source_filter}'")
    
    return filtered


def is_amyloid_positive(entry: Dict) -> bool:
    """Check if entry is classified as amyloid-forming."""
    # Primary check: 'classification' field (actual WALTZ-DB column)
    classification = str(entry.get('classification', '')).lower().strip()
    if classification == 'amyloid':
        return True
    if classification == 'non-amyloid':
        return False
    
    # Fallback checks for other possible field names
    for key in ['is_amyloid', 'amyloid', 'morphology_decision', 'morphology', 'classification', 'Classification']:
        val = str(entry.get(key, '')).lower().strip()
        if val in ['yes', 'true', '1', 'positive', 'amyloid', 'fibril', 'fiber']:
            return True
        if val in ['no', 'false', '0', 'negative', 'non-amyloid', 'non_amyloid']:
            return False
        if 'amyloid' in val and 'non' not in val:
            return True
    
    return False


def save_data(entries: List[Dict], output_path: str, format: str = 'tsv'):
    """Save data to file."""
    path = Path(output_path)
    
    if not entries:
        print("No data to save")
        return
    
    # Ensure classification column exists (compute if missing)
    for entry in entries:
        if 'classification' not in entry or not entry['classification']:
            entry['classification'] = 'amyloid' if is_amyloid_positive(entry) else 'non-amyloid'
    
    # Collect all keys
    all_keys = set()
    for entry in entries:
        all_keys.update(entry.keys())
    
    # Remove internal keys and sort
    all_keys = sorted(k for k in all_keys if not k.startswith('_'))
    
    # Prioritize important columns based on actual WALTZ-DB structure
    priority = ['sequence', 'classification', 'source', 'subset',
                'waltz', 'tango', 'pasta_parallel', 'pasta_antiparallel',
                'hydrophobicity', 'chou_fasman_helix', 'chou_fasman_strand',
                'uniprot_id', 'pdb_id', 'detail_url']
    headers = [k for k in priority if k in all_keys]
    headers += [k for k in all_keys if k not in headers]
    
    if format == 'json':
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(entries, f, indent=2, ensure_ascii=False)
    
    elif format in ['tsv', 'csv']:
        delimiter = '\t' if format == 'tsv' else ','
        with open(path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter=delimiter,
                                   extrasaction='ignore')
            writer.writeheader()
            writer.writerows(entries)
    
    print(f"Saved {len(entries)} entries to {output_path}")


def print_summary(entries: List[Dict]):
    """Print summary statistics."""
    if not entries:
        return
    
    print("\n" + "="*50)
    print("WALTZ-DB 2.0 Summary")
    print("="*50)
    print(f"Total peptides: {len(entries)}")
    
    # Count by classification
    amyloid_count = sum(1 for e in entries if is_amyloid_positive(e))
    non_amyloid_count = len(entries) - amyloid_count
    
    print(f"\nClassification:")
    print(f"  Amyloid-forming:    {amyloid_count:4d} ({100*amyloid_count/len(entries):.1f}%)")
    print(f"  Non-amyloidogenic:  {non_amyloid_count:4d} ({100*non_amyloid_count/len(entries):.1f}%)")
    
    # Sources distribution
    sources = {}
    for e in entries:
        src = e.get('source', 'Unknown')
        if not src:
            src = 'Unknown'
        sources[src] = sources.get(src, 0) + 1
    
    print(f"\nBy source (top 10):")
    for src, count in sorted(sources.items(), key=lambda x: -x[1])[:10]:
        print(f"  {src:20s}: {count:4d}")
    
    # Morphology decisions if available
    morphologies = {}
    for e in entries:
        morph = e.get('morphology_decision', e.get('morphology', ''))
        if morph:
            morphologies[morph] = morphologies.get(morph, 0) + 1
    
    if morphologies:
        print(f"\nBy morphology decision:")
        for morph, count in sorted(morphologies.items(), key=lambda x: -x[1]):
            print(f"  {morph:20s}: {count:4d}")
    
    # Sample sequences from each class
    print("\nSample amyloid-forming sequences:")
    amyloid_samples = [e for e in entries if is_amyloid_positive(e)][:5]
    for e in amyloid_samples:
        seq = e.get('sequence', '?')
        src = e.get('source', '?')[:15]
        print(f"  {seq:8s} (source: {src})")
    
    print("\nSample non-amyloidogenic sequences:")
    non_amyloid_samples = [e for e in entries if not is_amyloid_positive(e)][:5]
    for e in non_amyloid_samples:
        seq = e.get('sequence', '?')
        src = e.get('source', '?')[:15]
        print(f"  {seq:8s} (source: {src})")


def main():
    parser = argparse.ArgumentParser(
        description='Download peptide data from WALTZ-DB 2.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
    # Download ALL peptides (both amyloid and non-amyloid)
    python %(prog)s -o waltzdb_all.tsv
    
    # Only amyloid-forming peptides (positive class)
    python %(prog)s --amyloid-only -o waltzdb_amyloid.tsv
    
    # Only non-amyloidogenic peptides (negative class)  
    python %(prog)s --non-amyloid-only -o waltzdb_non_amyloid.tsv
    
    # Filter by source protein
    python %(prog)s --source "Sup35" -o waltzdb_sup35.tsv
    python %(prog)s --source "TDP-43" -o waltzdb_tdp43.tsv
    
    # With detailed information (slower)
    python %(prog)s --fetch-details -v

Data fields:
    sequence          - Hexapeptide sequence (6 amino acids)
    is_amyloid        - Amyloid classification (Yes/No)
    morphology        - Experimental morphology decision
    source            - Source protein/dataset
    waltz/tango/pasta - Prediction scores
    uniprot_id        - UniProt accession
    detail_url        - Link to full entry

Dataset composition (WALTZ-DB 2.0):
    - ~512 amyloid-forming hexapeptides
    - ~892 non-amyloidogenic hexapeptides
    - Sources: FUS, TDP-43, SOD-1, Sup35, functional amyloids, literature

Reference:
    Louros et al. (2020) NAR. doi:10.1093/nar/gkz758
'''
    )
    
    parser.add_argument('-o', '--output', default='waltzdb_peptides.tsv',
                        help='Output file path')
    parser.add_argument('-f', '--format', choices=['json', 'tsv', 'csv'],
                        default='tsv', help='Output format')
    parser.add_argument('--amyloid-only', action='store_true',
                        help='Only download amyloid-forming peptides')
    parser.add_argument('--non-amyloid-only', action='store_true',
                        help='Only download non-amyloidogenic peptides')
    parser.add_argument('--source', default=None,
                        help='Filter by source protein (e.g., "Sup35", "TDP-43")')
    parser.add_argument('--fetch-details', action='store_true',
                        help='Fetch detailed info for each peptide (slow)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    parser.add_argument('--max-pages', type=int, default=100,
                        help='Maximum pages to fetch')
    
    args = parser.parse_args()
    
    if not HAS_DEPS:
        return 1
    
    session = get_session()
    
    # Try direct CSV download first
    print("Checking for direct download option...")
    csv_data = try_direct_download(session, args.verbose)
    
    if csv_data:
        # Parse CSV data
        import io
        reader = csv.DictReader(io.StringIO(csv_data))
        entries = list(reader)
        print(f"Downloaded {len(entries)} entries directly")
    else:
        # Scrape HTML tables
        print("Direct download not available, scraping HTML...")
        entries = scrape_all_peptides(
            session, 
            verbose=args.verbose,
            fetch_details=args.fetch_details,
            max_pages=args.max_pages
        )
    
    if not entries:
        print("\nFailed to extract data. Suggestions:")
        print("1. Visit http://waltzdb.switchlab.org/sequences in browser")
        print("2. Use the CSV/Excel download buttons at bottom of table")
        print("3. Check if the site is accessible")
        return 1
    
    # Apply filters
    entries = filter_entries(entries, args.amyloid_only, args.non_amyloid_only, args.source)
    
    # Save
    save_data(entries, args.output, args.format)
    
    # Summary
    print_summary(entries)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())