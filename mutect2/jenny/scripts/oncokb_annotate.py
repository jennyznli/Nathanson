"""
Add OncoKB somatic oncogenicity / mutation effect / hotspot + gene type to the merged
VEP csv (output of annotate.smk merge_samples). Unique variants are queried once via
POST /annotate/mutations/byGenomicChange (GRCh38, no tumor type) and cached in a TSV
keyed on Chr_Start_REF_ALT, so reruns and other callers/labels only query new variants.

    python scripts/oncokb_annotate.py -i all_samples.mutect2.no_pon.vep.csv \
        -o all_samples.mutect2.no_pon.vep.oncokb.csv --cache oncokb/cache.GRCh38.tsv

Token: $ONCOKB_TOKEN, else --token-file (default ~/.oncokb_token). Never put it in config/git.
Stdlib only (urllib) so the vep_parser env needs nothing extra.
"""
import os
import csv
import json
import time
import argparse
import urllib.request
import urllib.error
from collections import Counter

API='https://www.oncokb.org/api/v1'
REFERENCE_GENOME='GRCh38'
BATCH_SIZE=100
MAX_RETRIES=5
ONCOKB_FIELDS=['OncoKB.Oncogenic','OncoKB.MutationEffect','OncoKB.Hotspot','OncoKB.VariantExist','OncoKB.Alteration']
OUT_FIELDS=ONCOKB_FIELDS+['OncoKB.GeneType','OncoKB.DataVersion']
CACHE_FIELDS=['key']+ONCOKB_FIELDS+['OncoKB.DataVersion']


def api_call(path,token,body=None):
    url=f'{API}{path}'
    data=json.dumps(body).encode() if body is not None else None
    headers={'Accept':'application/json','Authorization':f'Bearer {token}'}
    if data is not None:
        headers['Content-Type']='application/json'
    for attempt in range(MAX_RETRIES):
        req=urllib.request.Request(url,data=data,headers=headers)
        try:
            with urllib.request.urlopen(req,timeout=300) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            # 401/403 = bad token, 400 = bad query -- retrying won't help
            if e.code in (429,500,502,503,504) and attempt<MAX_RETRIES-1:
                wait=2**attempt*5
                print(f'  HTTP {e.code} on {path}, retrying in {wait}s')
                time.sleep(wait)
                continue
            raise SystemExit(f'OncoKB {path} failed: HTTP {e.code} {e.read()[:500]!r}')
        except urllib.error.URLError as e:
            if attempt<MAX_RETRIES-1:
                wait=2**attempt*5
                print(f'  {e.reason} on {path}, retrying in {wait}s')
                time.sleep(wait)
                continue
            raise SystemExit(f'OncoKB {path} unreachable: {e.reason}')


def read_token(token_file):
    token=os.environ.get('ONCOKB_TOKEN','').strip()
    if token:
        return token
    path=os.path.expanduser(token_file)
    if not os.path.isfile(path):
        raise SystemExit(f'No OncoKB token: set $ONCOKB_TOKEN or create {path} (chmod 600)')
    with open(path) as fh:
        return fh.read().strip()


def variant_key(row):
    return f"{row['Chr']}_{row['Start']}_{row['REF']}_{row['ALT']}"


def to_maf_location(chrom,pos,ref,alt):
    """
    VCF (anchor-base) -> TCGA MAF genomicLocation 'chr,start,end,ref,alt'.
    Insertions: ref '-', start/end flank the insertion. Deletions: alt '-', 1-based span of deleted bases.
    Returns None for symbolic/spanning alleles.
    """
    if alt in ['','.','*'] or alt.startswith('<') or ref in ['','.']:
        return None
    chrom=chrom[3:] if chrom.lower().startswith('chr') else chrom
    pos=int(pos)
    i=0
    while i<min(len(ref),len(alt)) and ref[i]==alt[i]:
        i+=1
    ref,alt,start=ref[i:],alt[i:],pos+i
    if not ref:
        return f'{chrom},{start-1},{start},-,{alt}'
    end=start+len(ref)-1
    return f'{chrom},{start},{end},{ref},{alt or "-"}'


def is_queryable(row):
    # protein-changing or splice -- everything else (intronic/UTR/intergenic) is left '.'
    return row['HGVSp'] not in ['','.'] or 'splice' in row['Variant.Consequence']


def load_cache(path,data_version):
    cache={}
    stale=0
    if not os.path.isfile(path):
        return cache,stale
    with open(path) as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            # re-query anything annotated under an older OncoKB data release
            if r.get('OncoKB.DataVersion')!=data_version:
                stale+=1
                continue
            cache[r['key']]=r
    return cache,stale


def save_cache(path,cache):
    os.makedirs(os.path.dirname(path) or '.',exist_ok=True)
    tmp=f'{path}.tmp.{os.getpid()}'
    with open(tmp,'w') as fh:
        w=csv.DictWriter(fh,fieldnames=CACHE_FIELDS,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        for k in sorted(cache):
            w.writerow(cache[k])
    os.replace(tmp,path)


def parse_response(r,data_version):
    q=r.get('query') or {}
    me=r.get('mutationEffect') or {}
    return {
        'key':q.get('id'),
        'OncoKB.Oncogenic':r.get('oncogenic') or '.',
        'OncoKB.MutationEffect':me.get('knownEffect') or '.',
        'OncoKB.Hotspot':str(r.get('hotspot')) if r.get('hotspot') is not None else '.',
        'OncoKB.VariantExist':str(r.get('variantExist')) if r.get('variantExist') is not None else '.',
        'OncoKB.Alteration':q.get('alteration') or '.',
        'OncoKB.DataVersion':data_version
    }


def get_args(argv):
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input_csv',required=True,help='Input: merged vep csv from annotate.smk merge_samples.')
    p.add_argument('-o','--output_csv',required=True,help='Output: input csv + OncoKB columns.')
    p.add_argument('-c','--cache',required=True,help='TSV cache of OncoKB results, shared across callers/labels.')
    p.add_argument('-t','--token-file',default='~/.oncokb_token',help='File holding the OncoKB API token (used if $ONCOKB_TOKEN unset).')
    return p.parse_args(argv)


def main(argv=None):
    args=get_args(argv)
    token=read_token(args.token_file)

    info=api_call('/info',token)
    data_version=(info.get('dataVersion') or {}).get('version','.')
    print(f'OncoKB data version: {data_version}')

    gene_type={g['hugoSymbol']:g.get('geneType') or '.' for g in api_call('/utils/cancerGeneList',token)}

    with open(args.input_csv) as fh:
        reader=csv.DictReader(fh)
        in_fields=reader.fieldnames
        rows=list(reader)

    locations={}
    for row in rows:
        if is_queryable(row):
            loc=to_maf_location(row['Chr'],row['Start'],row['REF'],row['ALT'])
            if loc:
                locations[variant_key(row)]=loc

    cache,stale=load_cache(args.cache,data_version)
    todo=[k for k in locations if k not in cache]
    print(f'Rows: {len(rows)}  unique queryable variants: {len(locations)}  '
          f'cache hits: {len(locations)-len(todo)}  stale cache entries: {stale}  to query: {len(todo)}')

    n_errors=0
    for b in range(0,len(todo),BATCH_SIZE):
        batch=todo[b:b+BATCH_SIZE]
        body=[{'id':k,'genomicLocation':locations[k],'referenceGenome':REFERENCE_GENOME} for k in batch]
        for r in api_call('/annotate/mutations/byGenomicChange',token,body):
            if r.get('errors'):
                n_errors+=1
            res=parse_response(r,data_version)
            if res['key'] in locations:
                cache[res['key']]=res
        print(f'  queried {min(b+BATCH_SIZE,len(todo))}/{len(todo)}')
        # save as we go so an interrupted run keeps what it already paid for
        save_cache(args.cache,cache)

    oncogenic_counts=Counter()
    out_fields=in_fields+[f for f in OUT_FIELDS if f not in in_fields]
    with open(args.output_csv,'w') as outfile:
        writer=csv.DictWriter(outfile,fieldnames=out_fields,delimiter=',',restval='.',extrasaction='ignore',quoting=csv.QUOTE_NONNUMERIC,dialect='excel')
        writer.writeheader()
        for row in rows:
            k=variant_key(row)
            hit=cache.get(k) if k in locations else None
            for f in ONCOKB_FIELDS:
                row[f]=hit[f] if hit else '.'
            row['OncoKB.DataVersion']=data_version if hit else '.'
            row['OncoKB.GeneType']=gene_type.get(row['Gene'],'.')
            oncogenic_counts[row['OncoKB.Oncogenic']]+=1
            writer.writerow(row)

    print(f'Responses with OncoKB errors: {n_errors}')
    print('OncoKB.Oncogenic (rows):')
    for k,v in oncogenic_counts.most_common():
        print(f'  {k}: {v}')
    print(f'{args.output_csv} written.')


if __name__=='__main__':
    main()
