# pdac_pathway_pubmed_enrichment.ipynb
#
# Jupyter Notebook style workflow for identifying pathways from a GSEA result list
# that are significantly associated with pancreatic ductal adenocarcinoma (PDAC)
# based on PubMed literature evidence.

import time
import math
import requests
import sys
from collections import namedtuple

try:
    from scipy.stats import fisher_exact
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

Result = namedtuple('Result', ['pathway', 'pathway_count', 'pdac_count', 'overlap', 'union_count', 'pvalue', 'adj_pvalue'])

ESEARCH_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
PDAC_QUERY = '(("pancreatic ductal adenocarcinoma"[Title/Abstract]) OR (PDAC[Title/Abstract]) OR ("pancreatic cancer"[Title/Abstract]))'
MIN_PAUSE = 0.35

# Utility: Query PubMed with E-utilities (esearch)
def esearch_count(query, email=None, api_key=None):
    params = {
        'db': 'pubmed',
        'term': query,
        'retmode': 'json',
        'rettype': 'count'
    }
    if email:
        params['email'] = email
    if api_key:
        params['api_key'] = api_key
    try:
        r = requests.get(ESEARCH_URL, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()
        return int(data.get('esearchresult', {}).get('count', 0))
    except Exception as e:
        print(f"Error querying PubMed: {e}", file=sys.stderr)
        return 0

# Helpers for formatting queries and computing p-values
def safe_term(ts):
    ts = ts.strip().replace('"', '')
    return f'"{ts}"[Title/Abstract]'


def compute_pvalue(k, n, K, M):
    a, b, c, d = k, n-k, K-k, M-n-K+k
    if any(val < 0 for val in (a, b, c, d)):
        return 1.0
    if SCIPY_AVAILABLE:
        try:
            _, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            return float(p)
        except Exception:
            return 1.0
    # Fallback: conservative hypergeometric tail
    def comb(N, k): return math.comb(N, k)
    max_k, min_k = min(n, K), max(0, n+K-M)
    if not (min_k <= k <= max_k):
        return 1.0
    denom = comb(M, n)
    tail = sum(comb(K, i) * comb(M-K, n-i) / denom for i in range(k, max_k+1))
    return min(1.0, 2*tail)

# Benjamini-Hochberg FDR correction
def benjamini_hochberg(pvalues):
    m = len(pvalues)
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adj = [0.0] * m
    prev_adj = 0.0
    for rank, (idx, p) in enumerate(indexed, start=1):
        adj_p = min(p * m / rank, 1.0)
        if adj_p < prev_adj:
            adj_p = prev_adj
        else:
            prev_adj = adj_p
        adj[idx] = adj_p
    return adj

# Step 1: Load pathways
input_file = "pathways.txt"   # replace with your file
with open(input_file, 'r', encoding='utf-8') as f:
    pathways = [line.strip() for line in f if line.strip()]
print(f"Loaded {len(pathways)} pathways")

# Step 2: Get PDAC count
pdac_count = esearch_count(PDAC_QUERY)
time.sleep(MIN_PAUSE)
print("PDAC PubMed count:", pdac_count)

# Step 3: Iterate through pathways
results, pvalues = [], []
for idx, pw in enumerate(pathways, start=1):
    print(f"[{idx}/{len(pathways)}] {pw}")
    pw_term = safe_term(pw)
    count_pw = esearch_count(pw_term)
    time.sleep(MIN_PAUSE)
    count_and = esearch_count(f'({pw_term}) AND ({PDAC_QUERY})')
    time.sleep(MIN_PAUSE)
    count_or = esearch_count(f'({pw_term}) OR ({PDAC_QUERY})')
    time.sleep(MIN_PAUSE)

    p = 1.0 if count_or <= 0 else compute_pvalue(count_and, count_pw, pdac_count, count_or)
    results.append(Result(pw, count_pw, pdac_count, count_and, count_or, p, None))
    pvalues.append(p)

# Step 4: Adjust p-values
adj = benjamini_hochberg(pvalues)
results = [r._replace(adj_pvalue=a) for r, a in zip(results, adj)]

# Filter significant results
threshold = 0.05
significant = [r for r in results if r.adj_pvalue <= threshold]

print(f"Significant pathways: {len(significant)}")
for r in significant:
    print(r.pathway, r.adj_pvalue)

# Step 5: Save output
output_file = "significant_pathways.txt"
with open(output_file, 'w', encoding='utf-8') as out:
    out.write(f"# Pathways significantly associated with PDAC (adj p <= {threshold})\n")
    out.write("# pathway\tpathway_count\tpdac_count\toverlap\tunion_count\tpvalue\tadj_pvalue\n")
    for r in sorted(significant, key=lambda x: x.adj_pvalue):
        out.write(f"{r.pathway}\t{r.pathway_count}\t{r.pdac_count}\t{r.overlap}\t{r.union_count}\t{r.pvalue:.3e}\t{r.adj_pvalue:.3e}\n")

print("Results written to:", output_file)
