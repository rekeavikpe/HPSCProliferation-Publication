import pathlib,tarfile
root=pathlib.Path('results/david_comments_20261020')
if not (root/'COMPUTATIONS_COMPLETE').exists(): raise SystemExit('Workflow has not finished. Check progress.log.')
archive=pathlib.Path('david_comments_20261020_for_review.tar.gz')
with tarfile.open(archive,'w:gz') as tar:
    for p in sorted(root.rglob('*')):
        if not p.is_file() or any(x.startswith('run') and x[3:].isdigit() for x in p.parts) or p.suffix=='.bin': continue
        tar.add(p,arcname=str(p))
print(archive.resolve())
