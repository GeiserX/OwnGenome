#!/usr/bin/env bash
# The Locityper output folders of the locityper-report cases: c<N>/loci/CYP2D6
# with res.json.gz and success as Locityper writes them, and the loci BED.
set -euo pipefail
python3 - <<'PY'
import gzip, json, os
def res(genotype='HG00733.1,NA19240.2', quality=27.8, total=900, unexpl=12, wd=0.2, warnings=None):
    r = {'total_reads': total, 'quality': quality, 'unexpl_reads': unexpl, 'genotype': genotype,
         'options': 120, 'dist_type': 'edit'}
    if wd is not None:
        r['weight_dist'] = wd
    if warnings is not None:
        r['warnings'] = warnings
    return r
def put(name, raw=None, success=True, cut=False, plain=False):
    d = f'{name}/loci/CYP2D6'
    os.makedirs(d)
    if raw is not None:
        body = raw if plain else gzip.compress(raw, mtime=0)
        with open(f'{d}/res.json.gz', 'wb') as f:
            f.write(body[:len(body) // 2] if cut else body)
    if success:
        open(f'{d}/success', 'w').close()
j = lambda r: json.dumps(r).encode()
put('clean', j(res()))
put('far', j(res(wd=31.5)))
put('unexplained', j(res(total=3000, unexpl=1200)))
put('share', j(res(total=10000, unexpl=1500)))
put('warned', j(res(warnings=['UnexplainedReads'])))
put('nodist', j(res(wd=None)))
put('nocounts', j({k: v for k, v in res().items() if k not in ('total_reads', 'unexpl_reads')}))
put('nulltotal', j(res(total=None)))
put('nosuccess', j(res()), success=False)
put('nores')
put('haploid', j(res(genotype='HG00733.1')))
put('notjson', b'{not json')
put('cut', j(res()), cut=True)
put('plain', j(res()), plain=True)
put('list', b'[1, 2]')
put('intgenotype', j(res(genotype=5)))
put('boolcounts', j(res(total=True, unexpl=False)))
put('strquality', j(res(quality='high', wd='x')))
put('nanquality', b'{"genotype": "A.1,B.2", "quality": NaN, "weight_dist": Infinity, "total_reads": 10, "unexpl_reads": 0}')
with open('loci.bed', 'w') as f:
    f.write('# a comment line\nchr22\t42122284\t42152960\tCYP2D6\n')
with open('three.bed', 'w') as f:
    f.write('# a comment line\n\nchr6\t32652603\t32669077\tHLA-DQB1\nchr22\t42122284\t42152960\tCYP2D6\tx\n'
            'chr1\t161506761\t161681047\tFCGR\n')
with open('bad.bed', 'w') as f:
    f.write('chr22\tx\t42152960\tCYP2D6\n')
with open('short.bed', 'w') as f:
    f.write('chr22\t42122284\t42152960\n')
with open('unicode.bed', 'w') as f:
    f.write('chr22\t\u00b2\t42152960\tCYP2D6\n')
PY
