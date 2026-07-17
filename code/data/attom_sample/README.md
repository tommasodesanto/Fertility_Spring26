# ATTOM assessor-history sample

Downloaded 2026-05-27 from ATTOM. 64 snappy parquet shards merged into a
single Stata `.dta` for quick exploration.

## Files

- `assessor_history.parquet` (~38 MB) — merged shards, snappy-compressed.
  Original schema: 279 columns, decimal128 numerics cast to float64.
- `assessor_history.dta` (~690 MB) — Stata 14+ (`version=118`) file.
  Open in Stata with `use code/data/attom_sample/assessor_history.dta, clear`.
- `merge_to_stata.py` — the conversion script. Source parquet shards live at
  `/Users/tommasodesanto/Desktop/019e6b85-e6d0-79f3-ac21-a5dc416c8dfa/`.

## Coverage (this sample)

- **Rows:** 368,751 parcel-year observations.
- **Geography:** essentially New York County, NY (Manhattan); 9 stray Michigan
  rows. Likely an ATTOM order scoped to NYC.
- **Year:** `assessorhistoryyear == 2025` for all rows.
- **Unit:** one row per `attomid` (ATTOM parcel ID).

## Schema notes

Column names lowercased and truncated to Stata's 32-char limit. Full original
names are preserved as Stata variable labels — view with `describe` or
`describe varname`. Nine columns were truncated; no collisions.

Useful starting variables:

| Stata name | What it is |
| --- | --- |
| `attomid` | ATTOM parcel identifier (key). |
| `assessorhistoryyear` | Assessment year. |
| `situsstatecode`, `situscounty`, `propertyaddresscity`, `propertyaddresszip` | Geography. |
| `latitude`, `longitude` | Parcel centroid (often missing for NYC condos). |
| `propertyusestandardized` | ATTOM standardized property-use code. |
| `taxassessedvaluetotal`, `taxassessedvalueland`, `taxassessedvalueimprovements` | Assessed values. |
| `taxmarketvaluetotal`, `taxbilledamount` | Market value, tax bill. |
| `deedlastsaleprice`, `deedlastsaledate` | Last recorded transfer. |
| `assessorlastsaleamount`, `assessorlastsaledate` | Assessor's last sale (often differs from deed). |
| `yearbuilt`, `areabuilding`, `arealotsf` | Structure characteristics. |
| `bedroomscount`, `bathcount`, `roomscount` | Interior counts. |
| `unitscount`, `storiescount` | Multifamily/structure detail. |

## Known limitations of this sample

- **Interior counts are sparse:** `bedroomscount`, `bathcount`, `roomscount`,
  `areabuilding` are zero/missing for most NYC condos and co-ops. ATTOM
  generally lacks unit-level interior data for high-rise residential in NYC.
  If you need interior counts, use suburban/single-family geographies.
- **62 columns are entirely null** in this sample (mostly rural/agricultural
  flags: grainery, silo, kennel, milkhouse, etc.). They survive in the file as
  empty strings so the column list matches the full ATTOM schema.
- **Single year, single county.** Don't read longitudinal patterns into this.

## Reproducing

```bash
python code/data/attom_sample/merge_to_stata.py
```

Requires the source parquet directory at the hard-coded `SRC_DIR` path. Edit
that constant if the shards move.
