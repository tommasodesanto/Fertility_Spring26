# Moving-item codes (verbatim from official PSID codebooks)

## Codebook access note
`https://psidonline.isr.umich.edu/documents/psid/codebook/*.pdf` returns a
Cloudflare "Security Challenge" bot-detection page to both `WebFetch` and a
browser-UA `curl` (confirmed for FAM2019ER, FAM1984, FAM1997ER). Per policy,
bot-detection is not bypassed. The same three PDFs were instead retrieved,
unmodified, from the Internet Archive Wayback Machine
(`web.archive.org/web/.../psidonline.isr.umich.edu/...`), which mirrors the
identical PDF without the challenge. Text below is copied verbatim from those
archived PDFs via `pdftotext -layout`; nothing is invented. The raw .dta's own
embedded Stata value-label metadata was also checked (`value_label_names()`):
only 41 of 87,146 raw variables carry any value label at all, and none of the
three target raw items are among them, so the .dta gives no shortcut here.

## MOVEDFREF_ raw item ("moved since ...")

**2019 (ER72156, "A49 WTR MOVED SINCE JAN 1 OF PRIOR YEAR"), from FAM2019ER_codebook.pdf:**
> A49. What is the street address and move-in date of (your/Reference Person's)
> current residence? (Have you/Has [he/she]) lived anywhere else since January
> 2017?
> 1 Yes | 5 No | 8 DK | 9 NA; refused

**1997 (ER10072, "A42 MOVED SINCE SPG?"), from FAM1997ER_codebook.pdf:**
> A42. Have you (HEAD) moved any time since the spring of 1996?
> 1 Yes | 5 No | 8 DK | 9 NA; refused

**1984 (V10447, "B27 MOVD SINCE SPR 1983?"), from FAM1984_codebook.pdf:**
> B27. Have you (HEAD) moved any time since the spring of 1983?
> 1 Yes | 5 No | 9 NA; DK

**Reference period, by era:**
- 1969-1997/2001 (annual, and the last two biennial waves before the wording
  changed): "since the spring of [prior survey year]" -- confirmed by the
  1984 and 1997 question text above, and by the header label pattern "MOVD/
  MOVED SINCE SPR..." or "MOVED SINCE SPG?" running continuously through 2001
  in the raw-file header scan.
- 2003-2019 (biennial): "since January 1 of [prior year]" -- confirmed by the
  2019 question text above; the raw-file label itself changes to "WTR MOVED
  SINCE JAN 1 OF PRIOR YEAR" starting exactly at the 2003 wave.
- The 1999 and 2001 waves were not independently codebook-checked (only a
  1997/1984/2019 codebook was fetched, per task scope); their raw label is
  still "A42 MOVED SINCE SPG?" (spring wording), so the switch to the
  January-1 reference period is inferred to fall at 2003, not earlier -- flag
  this as label-based, not codebook-quoted, for 1999/2001.
- 1968 has no "moved since" item: it is PSID's first wave, so there is no
  prior interview to reference.

**Code stability:** 1=Yes and 5=No are stable across all vintages checked.
The DK/NA split is NOT stable: 1984 pools DK and NA into a single code 9
("NA; DK"), while 1997 and 2019 split them into 8=DK and 9=NA;refused. Do not
assume 8/9 mean the same thing before 1997 without checking the specific
year's codebook.

## WHYMOVED1_ raw item ("why moved," first mention)

**2019 (ER72159, "A50 WHY MOVED 1ST MENTION"):**
> A50. Why did (you/he/she) move?--FIRST MENTION
> 1 Purposive productive reasons: to take another job; transfer; stopped going to school
> 2 To get nearer to work
> 3 Purposive consumptive reasons--expansion of housing: more space; more rent; better place
> 4 Purposive consumptive reasons--contraction of housing: less space; less rent
> 5 Purposive consumptive--other house-related: get own home/place; got married; physical conditions of the previous housing unit
> 6 Purposive consumptive--neighborhood-related: better neighborhood; go to school; to be closer to friends and/or relatives
> 7 Response to outside events (involuntary reasons): HU coming down; being evicted; armed services, etc.; health reasons; divorce; retiring because of health
> 8 Ambiguous, mixed, or other reasons, including reasons such as to save money, all my old neighbors moved away, retiring
> 9 Homeless
> 98 DK
> 99 NA; refused
> 0 Inap.: has not moved (ER72156=5); DK, NA, or RF whether moved (ER72156=8 or 9)

**1997 (ER10075, "A44 WHY MOVED 1ST"):**
> A44. Why did you (HEAD) move?--FIRST MENTION. The codes below are in priority order.
> 1 Purposive productive reasons: to take another job; transfer; stopped going to school
> 2 To get nearer to work
> 3 Purposive consumptive reasons--expansion of housing: more space; more rent; better place
> 4 Purposive consumptive reasons--contraction of housing: less space; less rent
> 5 Purposive consumptive--other house-related: want to own home; got married
> 6 Purposive consumptive--neighborhood-related: better neighborhood; go to school; to be closer to friends and/or relatives
> 7 Response to outside events (involuntary reasons): HU coming down; being evicted; armed services, etc.; health reasons; divorce; retiring because of health
> 8 Ambiguous or mixed reasons: to save money; all my old neighbors moved away; retiring (NA why)
> 99 DK; NA; refused
> 0 Inap.: has not moved

**1984 (V10449, "B29 WHY MOVED"):**
> B29. Why did you (HEAD) move? The codes below are in priority order.
> 1 Purposive productive reasons: to take another job; transfer; stopped going to school
> 2 To get nearer to work
> 3 Purposive consumptive reasons--expansion of housing: more space; better place
> 4 Purposive consumptive reasons--contraction of housing: less space; less rent
> 5 Purposive consumptive--other house-related: want to own home; got married
> 6 Purposive consumptive--neighborhood-related: better neighborhood; go to school
> 7 Response to outside events (involuntary reasons): HU coming down; being evicted; armed services, etc.; health reasons; divorce; retiring because of health
> 8 Ambiguous or mixed reasons: to save money; all my old neighbors moved away; retiring (NA why)
> 9 NA; DK
> 0 Inap.: has not moved (V10447=5 or 9)

**Vintage differences to respect:**
- Codes 1-8 (the substantive reason categories) are worded almost identically
  across 1984/1997/2019 -- stable enough to pool, with minor wording drift
  (e.g. category 5 adds "physical conditions of the previous housing unit" by
  2019; category 3 adds "more rent" by 1997/2019).
- Code 9 changes meaning across vintages: in 1984 and 1997, 9 (or 99) is the
  DK/NA catch-all; in 2019, 9 is reassigned to a new substantive category
  ("Homeless"), and DK/NA moves to 98/99. Any pooled 1968-2019 construction
  must NOT treat "9" as a constant DK/NA flag after this split -- check the
  exact vintage of each wave's why-moved variable before recoding.
- 0 = inapplicable (did not move) in every vintage checked.
- 1999-2019: up to 3 (2003-2007) or up to 4 (1997-2001, and 1999) "mentions"
  are coded as parallel variables (WHY MOVED 1ST/2ND/3RD/4TH...); the task's
  WHYMOVED1_ crosswalk uses only the first-mention variable, consistent with
  the shelf naming.

## HOMEOWN raw item (own/rent) -- not from a fetched PDF codebook

This one was not looked up in a PDF; it comes directly from the project's own
`/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/Construction_Files/Code/01 Collect housing variables.do`
recode (global `var_homeown`, executed for every wave):
```
replace homeown`y'=1 if inlist(`var', 1)   -> Owns home
replace homeown`y'=2 if inlist(`var', 5)   -> Pays rent
replace homeown`y'=3 if inlist(`var', 8)   -> Neither owns nor rents
replace homeown`y'=. if inlist(`var', 0, 9) | missing   -> NA/inap.
```
So the raw code list is 1=Owns, 5=Rents, 8=Neither owns nor rents (e.g. free
housing), with 0 and 9 (and missing) treated as NA. This matches the shelf
file's own embedded value-label set for `HOMEOWN` (`homeown_3cat`: 1="Owns
home", 2="Pays rent", 3="Neither owns home nor pays rent"), which is exactly
the do-file's 1/5/8 -> 1/2/3 recode -- an internal cross-check, not a second
independent source.
