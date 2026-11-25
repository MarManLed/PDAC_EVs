#!/usr/bin/env python3
"""Retrieve the public GEO processed 10x matrices and arrange PDAC1-6/ADJ1-6."""
from pathlib import Path
import tarfile, urllib.request, re, shutil, argparse
p=argparse.ArgumentParser()
p.add_argument('--archive',default='data/scRNAseq/GSE212966_RAW.tar')
p.add_argument('--download',action='store_true')
a=p.parse_args()
archive=Path(a.archive); archive.parent.mkdir(parents=True,exist_ok=True)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212966/suppl/GSE212966_RAW.tar'
if a.download and not archive.exists():
    print('Downloading ~440 MB from',url)
    urllib.request.urlretrieve(url,archive)
if not archive.exists(): raise SystemExit('Missing archive. Run with --download, or save GEO RAW.tar to '+str(archive))
base=archive.parent/'GSE212966_RAW'; base.mkdir(exist_ok=True)
with tarfile.open(archive) as tf:
    for member in tf:
        if not member.isfile(): continue
        name=Path(member.name).name
        match=re.search(r'(PDAC|ADJ)([1-6])',name,re.I)
        kind=re.search(r'(matrix\.mtx|genes\.tsv|features\.tsv|barcodes\.tsv)(\.gz)?$',name,re.I)
        if not match or not kind: continue
        tissue,n=match.group(1).upper(),match.group(2)
        folder=base/tissue/(tissue+n+'_filtered_feature_bc_matrix')
        folder.mkdir(parents=True,exist_ok=True)
        dest=folder/kind.group(0).lower()
        with tf.extractfile(member) as fin, dest.open('wb') as fout: shutil.copyfileobj(fin,fout)
for tissue in ('PDAC','ADJ'):
    dirs=sorted((base/tissue).glob('*_filtered_feature_bc_matrix'))
    print(tissue,len(dirs),'samples')
    for d in dirs:
        present={f.name for f in d.iterdir()}
        assert 'matrix.mtx.gz' in present and 'barcodes.tsv.gz' in present and ('genes.tsv.gz' in present or 'features.tsv.gz' in present),f'Missing components in {d}: {present}'
    assert len(dirs)==6, f'Expected six {tissue} samples; check raw archive names and metadata'
print('Done:',base)
