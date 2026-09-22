# National ACS comparison package

The primary readout is
[NATIONAL_ACS_PSID_FIRST_BIRTH_COMPARISON.md](NATIONAL_ACS_PSID_FIRST_BIRTH_COMPARISON.md).
The scientific comparison figure is
[national_acs_primary_vs_psid.png](national_acs_primary_vs_psid.png).

Small, independently reviewable fit receipts and tables are in
[national_continuation_20260921b/](national_continuation_20260921b/). The
directory contains no national panel or RDS checkpoint. The figure is
reproducible with:

    /opt/anaconda3/bin/python code/empirical/acs/kleven_pseudo/build_national_acs_psid_comparison.py

The package compares ACS full and reduced specifications with the saved PSID
rooms reference, PSID original ownership arm, and separately labeled aligned
ownership sensitivity. It preserves the distinction between event-cell
support/ESS and full-regression N, and makes no causal-validation claim.
