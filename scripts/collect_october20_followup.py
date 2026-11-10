import pathlib,tarfile
root=pathlib.Path('results/david_october20_followup')
if not (root/'COMPUTATIONS_COMPLETE').is_file():raise SystemExit('Follow-up has not completed')
out=pathlib.Path('david_october20_followup_for_review.tar.gz')
with tarfile.open(out,'w:gz') as tar:
 for p in sorted(root.rglob('*')):
  if not p.is_file() or p.suffix=='.bin' or any(x.startswith('run') and x[3:].isdigit() for x in p.parts):continue
  tar.add(p,arcname=str(p))
print(out.resolve())
