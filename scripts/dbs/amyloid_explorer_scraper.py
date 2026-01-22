#!/usr/bin/env python3
"""
Amyloid Explorer Database Scraper v2
=====================================
Downloads data from the Amyloid Explorer database.
https://amyloid-explorer.switchlab.org/ (or https://stamp.switchlab.org/)

Key improvements over v1:
- Focus on Playwright with network interception (React SPAs need this)
- Better handling of dynamically loaded data
- JSON-based data extraction from XHR/Fetch requests
- Fallback to DOM scraping with improved selectors

Requirements:
    pip install playwright
    playwright install chromium

Usage:
    python amyloid_explorer_scraper_v2.py
    python amyloid_explorer_scraper_v2.py --output data.json --verbose
    python amyloid_explorer_scraper_v2.py --pdb-list  # Extract just PDB IDs

Author: Generated for Xenia's amyloid research
Date: 2026-01-19
"""

import argparse
import json
import sys
import time
import re
from pathlib import Path
from typing import Optional, Dict, List, Any
from dataclasses import dataclass, asdict
from datetime import datetime

BASE_URL = "https://amyloid-explorer.switchlab.org"
ALT_URL = "https://stamp.switchlab.org"


@dataclass
class FibrilEntry:
    """Structure for a single fibril entry."""
    pdb_id: str
    protein_name: str = ""
    method: str = ""
    resolution: str = ""
    origin: str = ""  # recombinant, ex vivo, seeded
    organism: str = ""
    tissue: str = ""
    mutations: str = ""
    ptms: str = ""
    emdb_id: str = ""
    bmrb_id: str = ""
    # Additional metadata
    raw_data: dict = None
    
    def to_dict(self):
        d = asdict(self)
        if d['raw_data'] is None:
            del d['raw_data']
        return d


def scrape_with_playwright(base_url: str, verbose: bool = False, 
                           timeout: int = 60000) -> Dict[str, Any]:
    """
    Use Playwright to scrape the SPA with network request interception.
    
    This is the primary method for React SPAs - we intercept the actual
    API calls the frontend makes to get the data.
    """
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("ERROR: Playwright not installed.")
        print("Install with: pip install playwright && playwright install chromium")
        return None
    
    # Storage for intercepted data
    api_responses = []
    captured_requests = []
    
    def handle_response(response):
        """Intercept and capture API responses."""
        url = response.url
        content_type = response.headers.get('content-type', '')
        
        # Capture potentially interesting responses
        is_api = any(x in url.lower() for x in [
            'api', 'data', '.json', 'graphql', 'query', 
            'fibril', 'structure', 'protein', 'pdb',
            'supabase', 'firebase', 'hasura'  # Common backend services
        ])
        is_json = 'json' in content_type
        
        if is_api or is_json:
            try:
                data = response.json()
                size = len(json.dumps(data)) if data else 0
                
                api_responses.append({
                    'url': url,
                    'content_type': content_type,
                    'data': data,
                    'size': size,
                    'status': response.status
                })
                
                if verbose:
                    print(f"  [API] {response.status} {url[:80]}...")
                    print(f"        Size: {size} bytes, Type: {content_type[:30]}")
                    
            except Exception as e:
                if verbose:
                    print(f"  [API] Could not parse JSON from {url}: {e}")
    
    def handle_request(request):
        """Log outgoing requests for debugging."""
        if verbose and any(x in request.url.lower() for x in ['api', 'data', 'graphql']):
            captured_requests.append({
                'url': request.url,
                'method': request.method,
                'headers': dict(request.headers)
            })
    
    print(f"Starting browser automation for {base_url}...")
    
    results = {
        'api_data': None,
        'table_data': [],
        'card_data': [],
        'pdb_ids': [],
        'raw_responses': []
    }
    
    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=['--no-sandbox', '--disable-dev-shm-usage']
        )
        
        context = browser.new_context(
            viewport={'width': 1920, 'height': 1080},
            user_agent='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        )
        
        page = context.new_page()
        
        # Set up interceptors
        page.on('response', handle_response)
        if verbose:
            page.on('request', handle_request)
        
        # Navigate to the database page
        database_url = f"{base_url}/database"
        print(f"Navigating to: {database_url}")
        
        try:
            page.goto(database_url, timeout=timeout, wait_until='networkidle')
            print("Page loaded, waiting for data...")
        except Exception as e:
            print(f"Navigation note (SPAs often timeout): {e}")
        
        # Wait for React to render
        time.sleep(5)
        
        # Try clicking various UI elements to trigger data loading
        print("Attempting to trigger data loading...")
        
        load_triggers = [
            'button:has-text("Load")',
            'button:has-text("All")',
            'button:has-text("Show")',
            '[data-testid*="load"]',
            '[class*="load-more"]',
            '[class*="pagination"] button',
            'select option:has-text("All")',  # Page size selector
            '[aria-label*="page"]',
        ]
        
        for selector in load_triggers:
            try:
                elements = page.query_selector_all(selector)
                for el in elements[:3]:  # Limit clicks
                    try:
                        el.click()
                        time.sleep(1)
                    except:
                        pass
            except:
                pass
        
        # Wait for any additional data
        time.sleep(3)
        
        # Extract data from DOM
        print("Extracting data from page DOM...")
        
        # Method 1: Look for tables
        table_data = page.evaluate('''() => {
            const results = [];
            const tables = document.querySelectorAll('table');
            
            tables.forEach((table, tableIdx) => {
                const headers = [];
                const rows = [];
                
                // Get headers
                table.querySelectorAll('thead th, tr:first-child th').forEach(th => {
                    headers.push(th.textContent.trim());
                });
                
                // Get data rows
                const dataRows = table.querySelectorAll('tbody tr');
                dataRows.forEach(row => {
                    const cells = row.querySelectorAll('td');
                    if (cells.length > 0) {
                        const rowData = { _source: 'table_' + tableIdx };
                        cells.forEach((cell, idx) => {
                            const key = headers[idx] || 'col_' + idx;
                            rowData[key] = cell.textContent.trim();
                            
                            // Extract PDB links
                            const pdbLink = cell.querySelector('a[href*="rcsb.org"], a[href*="pdb"]');
                            if (pdbLink) {
                                rowData[key + '_pdb_url'] = pdbLink.href;
                                const pdbMatch = pdbLink.href.match(/([A-Za-z0-9]{4})$/);
                                if (pdbMatch) rowData['pdb_id'] = pdbMatch[1].toUpperCase();
                            }
                            
                            // Extract any PDB ID from text (4 alphanumeric chars)
                            const pdbText = cell.textContent.match(/\\b([0-9][A-Za-z0-9]{3})\\b/);
                            if (pdbText && !rowData['pdb_id']) {
                                rowData['pdb_id'] = pdbText[1].toUpperCase();
                            }
                        });
                        rows.push(rowData);
                    }
                });
                
                if (rows.length > 0) {
                    results.push({ headers, rows, table_index: tableIdx });
                }
            });
            
            return results;
        }''')
        
        if table_data:
            for table in table_data:
                results['table_data'].extend(table.get('rows', []))
            print(f"  Found {len(results['table_data'])} table rows")
        
        # Method 2: Look for cards/list items
        card_data = page.evaluate('''() => {
            const results = [];
            
            // Common card/item selectors for React apps
            const selectors = [
                '[class*="card"]',
                '[class*="item"]',
                '[class*="entry"]',
                '[class*="fibril"]',
                '[class*="structure"]',
                '[class*="row"]:not(thead *):not(tbody *)',
                '[data-testid*="entry"]',
                '[data-testid*="row"]',
                'li[class*="list"]'
            ];
            
            const seen = new Set();
            
            selectors.forEach(selector => {
                document.querySelectorAll(selector).forEach(item => {
                    // Avoid duplicates
                    const key = item.textContent.trim().substring(0, 100);
                    if (seen.has(key) || key.length < 10) return;
                    seen.add(key);
                    
                    const data = {
                        text: item.textContent.trim().substring(0, 1000),
                        className: item.className,
                        _source: 'card'
                    };
                    
                    // Get data attributes
                    Array.from(item.attributes).forEach(attr => {
                        if (attr.name.startsWith('data-')) {
                            data[attr.name] = attr.value;
                        }
                    });
                    
                    // Extract PDB IDs from text
                    const pdbMatches = item.textContent.match(/\\b([0-9][A-Za-z0-9]{3})\\b/g);
                    if (pdbMatches) {
                        data['pdb_ids'] = [...new Set(pdbMatches.map(p => p.toUpperCase()))];
                    }
                    
                    // Extract links
                    const links = item.querySelectorAll('a[href]');
                    if (links.length > 0) {
                        data['links'] = Array.from(links).slice(0, 10).map(l => ({
                            text: l.textContent.trim(),
                            href: l.href
                        }));
                    }
                    
                    results.push(data);
                });
            });
            
            return results;
        }''')
        
        if card_data:
            results['card_data'] = card_data
            print(f"  Found {len(card_data)} card/item elements")
        
        # Method 3: Extract all PDB IDs from the page
        all_pdb_ids = page.evaluate('''() => {
            const text = document.body.textContent;
            // PDB IDs are 4 characters: 1 digit followed by 3 alphanumeric
            const matches = text.match(/\\b([0-9][A-Za-z0-9]{3})\\b/g) || [];
            // Filter to likely PDB IDs and deduplicate
            return [...new Set(matches.map(m => m.toUpperCase()))];
        }''')
        
        results['pdb_ids'] = all_pdb_ids
        print(f"  Found {len(all_pdb_ids)} potential PDB IDs")
        
        # Store captured API responses
        results['raw_responses'] = api_responses
        
        # Find the best API response (likely the main data)
        if api_responses:
            # Sort by size and find ones that look like lists
            for resp in sorted(api_responses, key=lambda x: x['size'], reverse=True):
                data = resp['data']
                if isinstance(data, list) and len(data) > 5:
                    results['api_data'] = data
                    print(f"  Found main API data: {len(data)} entries from {resp['url'][:60]}...")
                    break
                elif isinstance(data, dict):
                    # Check for common patterns
                    for key in ['data', 'results', 'entries', 'fibrils', 'structures', 'items']:
                        if key in data and isinstance(data[key], list):
                            results['api_data'] = data[key]
                            print(f"  Found API data in '{key}': {len(data[key])} entries")
                            break
                    if results['api_data']:
                        break
        
        browser.close()
    
    return results


def process_results(results: Dict[str, Any]) -> List[Dict]:
    """
    Process and merge data from different extraction methods.
    """
    entries = []
    seen_pdb = set()
    
    # Priority 1: API data (most reliable)
    if results.get('api_data'):
        for item in results['api_data']:
            if isinstance(item, dict):
                entries.append(item)
                # Track PDB IDs
                for key in ['pdb_id', 'pdbId', 'pdb', 'PDB', 'id']:
                    if key in item:
                        seen_pdb.add(str(item[key]).upper())
    
    # Priority 2: Table data
    if results.get('table_data') and not entries:
        entries = results['table_data']
        for item in entries:
            if 'pdb_id' in item:
                seen_pdb.add(item['pdb_id'])
    
    # Priority 3: Card data (fallback)
    if not entries and results.get('card_data'):
        entries = results['card_data']
    
    return entries, list(seen_pdb)


def save_data(data: Any, output_path: str, format: str = 'json') -> None:
    """Save data to file."""
    output = Path(output_path)
    
    if format == 'json':
        with open(output, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False, default=str)
    
    elif format == 'tsv':
        import csv
        if isinstance(data, list) and len(data) > 0:
            # Get all possible keys
            all_keys = set()
            for item in data:
                if isinstance(item, dict):
                    all_keys.update(item.keys())
            
            # Remove internal keys
            all_keys = {k for k in all_keys if not k.startswith('_')}
            headers = sorted(all_keys)
            
            with open(output, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t',
                                       extrasaction='ignore')
                writer.writeheader()
                for item in data:
                    if isinstance(item, dict):
                        # Convert complex values to strings
                        row = {}
                        for k, v in item.items():
                            if k.startswith('_'):
                                continue
                            if isinstance(v, (list, dict)):
                                row[k] = json.dumps(v)
                            else:
                                row[k] = v
                        writer.writerow(row)
    
    elif format == 'pdb_list':
        # Just PDB IDs, one per line
        if isinstance(data, list):
            with open(output, 'w') as f:
                for pdb in sorted(set(data)):
                    f.write(f"{pdb}\n")
    
    print(f"Saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Download data from Amyloid Explorer database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
    # Full scrape with default settings
    python %(prog)s
    
    # Extract just PDB IDs
    python %(prog)s --pdb-list -o pdb_ids.txt
    
    # Verbose mode
    python %(prog)s -v -o amyloid_data.json
    
    # Use alternative URL
    python %(prog)s --url https://stamp.switchlab.org

Requirements:
    pip install playwright
    playwright install chromium
'''
    )
    
    parser.add_argument('-o', '--output', default='amyloid_explorer_data.json',
                        help='Output file path')
    parser.add_argument('-f', '--format', choices=['json', 'tsv', 'pdb_list'],
                        default='json', help='Output format')
    parser.add_argument('--pdb-list', action='store_true',
                        help='Extract only PDB IDs (shortcut for --format pdb_list)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    parser.add_argument('--url', default=BASE_URL,
                        help=f'Base URL (default: {BASE_URL})')
    parser.add_argument('--timeout', type=int, default=60000,
                        help='Page load timeout in ms (default: 60000)')
    parser.add_argument('--save-raw', action='store_true',
                        help='Save raw API responses for debugging')
    
    args = parser.parse_args()
    
    if args.pdb_list:
        args.format = 'pdb_list'
        if args.output == 'amyloid_explorer_data.json':
            args.output = 'amyloid_explorer_pdb_ids.txt'
    
    # Try primary URL
    results = scrape_with_playwright(args.url, args.verbose, args.timeout)
    
    # Fall back to alternative URL
    if not results or (not results.get('api_data') and not results.get('table_data')):
        alt_url = ALT_URL if args.url == BASE_URL else BASE_URL
        print(f"\nTrying alternative URL: {alt_url}")
        results = scrape_with_playwright(alt_url, args.verbose, args.timeout)
    
    if not results:
        print("\nFailed to extract data. Suggestions:")
        print("1. Check if the website is accessible in your browser")
        print("2. Open browser DevTools > Network tab and manually identify API endpoints")
        print("3. Contact authors: nikolaos.louros@utsouthwestern.edu")
        return 1
    
    # Process results
    entries, pdb_ids = process_results(results)
    
    # Save raw responses for debugging
    if args.save_raw and results.get('raw_responses'):
        raw_path = Path(args.output).stem + '_raw.json'
        save_data(results['raw_responses'], raw_path)
        print(f"Raw API responses saved to {raw_path}")
    
    # Determine what to save
    if args.format == 'pdb_list':
        # Combine PDB IDs from all sources
        all_pdbs = set(pdb_ids)
        for entry in entries:
            if isinstance(entry, dict):
                for key in ['pdb_id', 'pdbId', 'pdb', 'PDB']:
                    if key in entry:
                        all_pdbs.add(str(entry[key]).upper())
        save_data(list(all_pdbs), args.output, 'pdb_list')
        print(f"\nExtracted {len(all_pdbs)} unique PDB IDs")
    else:
        if entries:
            save_data(entries, args.output, args.format)
            print(f"\nExtracted {len(entries)} entries")
            if entries and isinstance(entries[0], dict):
                print(f"Sample fields: {list(entries[0].keys())[:8]}")
        else:
            print("\nNo structured data extracted.")
            print("Raw PDB IDs found:", pdb_ids[:20], "..." if len(pdb_ids) > 20 else "")
    
    # Summary
    print("\n=== Summary ===")
    print(f"API data entries: {len(results.get('api_data') or [])}")
    print(f"Table rows: {len(results.get('table_data', []))}")
    print(f"Card elements: {len(results.get('card_data', []))}")
    print(f"PDB IDs found: {len(pdb_ids)}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
