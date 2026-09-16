### Run provenance

- **Reference bands:** PLAbDab (all pairings) + PLAbDab/TheraSAbDab (clinical-stage) for vh_germline_identity, vl_germline_identity (fit n=29147, validation n=29030, cuts p5/p95/p99/p99.9)
- **Numbering:** 20/20 candidates numbered (100.0%)
- **Metrics with no reference band (unmeasured, NOT typical):** none

| Candidate | Level | What to do | Driver |
|---|---|---|---|
| AFO05311_AFO05305 | **L1 Ready** | put it on the wet-lab list as it stands | oxidation-prone M at VH FR2 53 — one framework substitution |
| Micvotabart | **L1 Ready** | put it on the wet-lab list as it stands | isomerization motif 'DS' at VH FR3 69 — one framework substitution |
| Tamgiblimab | **L1 Ready** | put it on the wet-lab list as it stands | oxidation-prone M at VH FR2 39 — one framework substitution |
| Conatumumab | **L2 Point repair** | one or two substitutions first; re-measure affinity if any is in a CDR | deamidation motif 'NS' at VH CDR2 58 — one CDR substitution + re-measure affinity |
| Izastobart | **L2 Point repair** | one or two substitutions first; re-measure affinity if any is in a CDR | isomerization motif 'DS' at VL CDR3 105 — one CDR substitution + re-measure affinity |
| Petosemtamab_2 | **L2 Point repair** | one or two substitutions first; re-measure affinity if any is in a CDR | deamidation motif 'NT' at VH CDR2 57 — one CDR substitution + re-measure affinity |
| QHZ08643_QHZ08651 | **L2 Point repair** | one or two substitutions first; re-measure affinity if any is in a CDR | deamidation motif 'NG' at VH CDR1 35 — one CDR substitution + re-measure affinity |
| UOS82608_UOS82606 | **L2 Point repair** | one or two substitutions first; re-measure affinity if any is in a CDR | deamidation motif 'NG' at VH CDR2 62 — one CDR substitution + re-measure affinity |
| AFQ83141_AFQ83148 | **L3 Redesign** | a loop or a surface has to be redesigned — budget a design cycle, not a mutation | unusual: vh_aggregation=0.2576 — above the p5-p95 range [0.1559, 0.2436] of PLAbDab (all pairings) (n=29147) [AGGRESCAN Na4vSS — a regional property, not one residue] |
| Lutikizumab | **L3 Redesign** | a loop or a surface has to be redesigned — budget a design cycle, not a mutation | unusual: cdr_net_charge=2.268 — above the p5-p95 range [-5.069, 2.159] of PLAbDab (all pairings) (n=29147) |
| XEF76981_XEF76997 | **L3 Redesign** | a loop or a surface has to be redesigned — budget a design cycle, not a mutation | VH has an odd cysteine count (3) — at least one cannot be canonically paired — one CDR substitution + re-measure affinity |
| Anivovetmab | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | VH closest germline is mouse, not human (IGHV2-6*02) — outside this pipeline's declared scope; every band was fitted on human antibodies and does not describe this — redesign a loop or a surface — not a substitution |
| Donitabart | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | strongly unusual: cdrh3_len=6 — below the p5-p95 range [9, 22] of PLAbDab (all pairings) (n=29147) |
| Ibentatug | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | VH closest germline is mouse, not human (IGHV5-9-3*01) — outside this pipeline's declared scope; every band was fitted on human antibodies and does not describe this — redesign a loop or a surface — not a substitution |
| QKO95500_QKO95510 | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | strongly unusual: cdrl3_len=14 — above the p5-p95 range [8, 12] of PLAbDab (all pairings) (n=29147) |
| QRK71970_QRK71974 | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | VH closest germline is mouse, not human (IGHV1-69-1*01) — outside this pipeline's declared scope; every band was fitted on human antibodies and does not describe this — redesign a loop or a surface — not a substitution |
| UOS76540_UOS76539 | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | VH closest germline is mouse, not human (IGHV1-9*02) — outside this pipeline's declared scope; every band was fitted on human antibodies and does not describe this — redesign a loop or a surface — not a substitution |
| Umizortamig1_2 | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | VH closest germline is mouse, not human (IGHV3-2*02) — outside this pipeline's declared scope; every band was fitted on human antibodies and does not describe this — redesign a loop or a surface — not a substitution |
| WFI50829_WFI50834 | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | strongly unusual: fv_net_charge=-9.035 — below the p5-p95 range [-4.0487, 5.848] of PLAbDab (all pairings) (n=29147) |
| WYM58511_WYM58510 | **L4 Human decision** | nothing routes this automatically; a person decides whether it is worth the bench time | VH closest germline is mouse, not human (IGHV1-47-2*01) — outside this pipeline's declared scope; every band was fitted on human antibodies and does not describe this — redesign a loop or a surface — not a substitution |
