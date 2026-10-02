# Saved-path terminal diagnosis

The reviewed policy source is `/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a`. This readout advances **saved policies and distributions only**; it does not run a Bellman problem, stationary equilibrium, or dated market root. All twelve mappings had passed their dated market/fiscal roots. The six originally reported failing cases remain failed terminal-gate cases; the recovered numbers are diagnostic, not accepted policy effects.

`recover.py` loads each final mapping's SHA-checked last-date diagnostic packet, advances its `evaluation` with the native `advance_sequential_calendar_distribution`, and adds the entrant cohort using the mapping's recorded flow. It seeds both birth-entry queues from the SHA-checked stationary 80% reference packet and advances them using every saved adjusted and raw birth row. Each due flow matches its saved mapping row exactly; the forward mass-accounting gap is below `1e-10`. For permanent cases, the accepted fresh 100% terminal one-step packet supplies the stationary endpoint. It evaluates the unchanged native `terminal_convergence_diagnostics` at the production tolerance `1e-3`, then applies the unchanged raw-queue gate. The three forward-helper modules load from the frozen project tree mounted by the production launcher. Their SHA-256 hashes match `current_source_files` in every case's saved native preparation receipt, and none is replaced by the production base, floor, or reviewed mechanism overlays. Missing or differing receipt pins stop the replay.

All four 80% controls reproduce every metric in their official `completed.json` within `1e-10`, including the raw queue metric, using the stricter `1e-6` control tolerance. The remaining cases fail as follows (all numbers are relative gaps except normalized distribution L1):

| Arm | Policy path | Horizon | Normalized distribution L1 | Renter-price gap | Other failed gates |
| --- | --- | ---: | ---: | ---: | --- |
| Hard | Temporary | 12 | 0.002710 | 0.001704 | None |
| Hard | Temporary | 16 | 0.002909 | 0.002055 | None |
| Quarter-saving | Temporary | 12 | 0.002457 | 0.001358 | None |
| Quarter-saving | Temporary | 16 | 0.002765 | 0.001749 | None |
| Hard | Permanent | 12 | 0.064322 | 0.040347 | Population 0.017492, adjusted/raw queues, asset price |
| Hard | Permanent | 16 | 0.053655 | 0.029150 | Population 0.012765, adjusted/raw queues, asset price |
| Quarter-saving | Permanent | 12 | 0.057573 | 0.036760 | Population 0.015835, adjusted/raw queues, asset price |
| Quarter-saving | Permanent | 16 | 0.046129 | 0.026823 | Population 0.011626, adjusted/raw queues, asset price |

The complete metrics, source paths, and packet digests are in `hard_recovered_terminal.json` and `quarter_recovered_terminal.json`. The temporary paths therefore miss the terminal distribution and rental-price gates even though their population, adjusted/raw queue, asset-price, and preference gates pass. The permanent paths remain substantially farther from their new stationary endpoints. The terminal gate compares the final endogenous date's rent with the stationary rent, so the rental gap is not caused by an omitted forward-distribution reconstruction.

To reproduce on Torch without model solves, stage `recover.py` at `/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a/terminal_diagnosis_recover.py` and run it separately for `hard` and `quarter` using `/share/apps/anaconda3/2025.06/bin/python recover.py OUT ARM`, with one BLAS/Numba thread. The resulting JSONs were copied into this folder. No production source or policy result was modified.
