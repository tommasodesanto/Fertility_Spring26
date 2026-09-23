**APPROVE.** No replacement needed.

- **Wording:** The new clause ("Housing moments use the national household sample in the 2005--2006 ACS") changes only the geography, from 42 metros to the whole U.S. The years and the row definitions stay as they were, which is what the author decided. "Household" is the right unit for all four ACS rows; ownership is measured for household heads aged 30–55. "Use" matches "wealth moments use the 2003 and 2005 PSID" later in the same sentence. The clause is also consistent with the `recomputed.national` block of the receipt and with the September 23 14:31 entry in `CALIBRATION_STATUS.md`.
- **Minimality:** This is the smallest appropriate edit. That clause was the only mention of geography in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_mock/sections/04_quantification.tex:13`. Lines 204–206, which say which moments come from the ACS and which from the PSID event study, need no change. The fit-table footnote at lines 208–214 still correctly says the weighting and other row choices are under review.
- **Values:** I compared the table against `recomputed.national`:

| Row | Table | National value in the file |
|---|---|---|
| Mean rooms | 5.608 | 5.607886 ✓ |
| Ownership, ages 30–55 | 67.63% | 0.676260 ✓ |
| Family-room gap | 0.385 | 0.385100 ✓ |
| Recent-parent ownership gap | 12.76 pp | 0.127608 ✓ |

  None of the four is the old 42-metro value (5.561, 64.83, 0.347, 16.29). The model column is still "---" and the other nine rows are unchanged.

**Caution:** In `target_recomputed.json`, the top-level `target_values` and `gate` blocks still hold the 42-metro numbers. The national numbers appear only under `recomputed.national`, so anyone pulling the targets from that file again needs to read that key.

**Scope:** I read only this mock section, so I didn't check other mock sections for 42-metro wording. There was no `memory/daily/` file, so for startup I used `AGENT_MEMORY.md` and `CALIBRATION_STATUS.md`. I made no writes and ran nothing.
