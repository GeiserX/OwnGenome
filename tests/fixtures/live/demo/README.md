# tests/fixtures/live/demo

Invented inputs of the live-site demo (`tests/demo/build_live_demo.py`), added to the invented
sample DEMO-001 that `tests/demo/make_demo_sample.py` writes. No real sample, person or run.

- `DEMO-001/recheck/state.json`: a database re-check as `bin/owngenome_recheck.py` writes it: a newer
  ClinVar release and an Ensembl cache with nothing newer. The build sets `checked_utc` to one hour
  after the replayed run ends and each item's `steps` with `owngenome_recheck.affected_steps`, so the
  badges name the steps the real check would name.
- `DEMO-FAMILY/members.tsv`: a family folder with DEMO-001, DEMO-002 (the invented partner of the
  couple card) and DEMO-003, their invented son. The build writes DEMO-002 and DEMO-003 with
  `make_demo_sample.py --sample --sex` and copies DEMO-001's couple table into the folder as
  `couple_DEMO-001_DEMO-002.tsv`, as step 39 names it.
- `DEMO-FAMILY/family_plan.tsv`, `relationships.tsv`, `ibd_pairs.tsv`, `DEMO-003_segments.tsv`,
  `de_novo.tsv`, `compound_het.tsv`: step 38's files for the trio, invented, in the columns step 38
  writes, so the family page draws its grid, the son's painted stretches and his tables.
  `compound_het.tsv` also holds a pair in KCNH2 whose first variant is the son's planted result to act on
  (`make_demo_sample.py` PLANTED): the family page keeps the son's pairs closed until his results are opened there.

The run itself is `../demo-001.weblog.jsonl`, written by `tests/demo/make_demo_events.py`.
