# Lead correction

The calculation docstring and family-link comment now identify `MOMLOC`/`POPLOC` links to the head as the primary child count, with `RELATE=3` described only as a secondary check. Documentation clarifies that the all-age validation is a match of counts capped at 3; it does not establish exact counts above 3 because `NCHILD` is top-coded.

The script now asserts unique head household keys, all four existing target gaps below `1e-9`, zero disagreement between parent-linked child bins and `NCHILD` capped at 3, and complete coverage of the family sample by the two candidate groups. Its `PASS` print occurs after these gates. The stored aggregate JSON satisfies the three numerical/sample assertions; the script compiles in memory. No raw ACS rerun was needed for this correction. The unique-key check is now enforced by the script on a future run; the stored aggregate does not retain household keys for a separate post hoc uniqueness test.

The lead subsequently reran the bounded memory-mapped calculation with all four new gates, including unique head keys. All gates passed and the same four target values and candidate shift were reproduced. This supersedes the earlier no-rerun limitation for the unique-key assertion; the raw source still was not independently rehashed.
