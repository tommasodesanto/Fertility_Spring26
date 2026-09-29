# Codex worker task

## Goal

Fix concrete lead-reviewed issues in the prepared isolated packet, before any
test submission. This is a changed targeted scope, not an unchanged retry.

## Scope

Exclusive write ownership same credit_no_taper_v1 folder. No shared edits.
Read existing driver/tests/README and correct only following issues:
1. Torch compute hostname is `cs603` etc, so guard requiring substring'torch'
   will fail. Require Linux+SLURM and authenticated exact frozen source; don't
   invent hostname assumptions.
2. verify_source currently compares manifest physical `/scratch/.../project`
   to container-bound `/Users/...` root. Permit exactly known physical OR the
   established container mount, while validating manifest physicalrootconstant
   and exact source SHA. Never weaken source/reference validation.
3. Native owner floor=min(b,collateral), so testexpected for b=+1,collateral=-1
   is -1, not+1. Test actual native buyer owner_floor with native_purchase_income
   plus incumbent rule, both newflag on/off, not merely an unrelatedformula.
4. README incorrectly writes b'= rather than b'>=; fixboundsnotation.
5. Explicitly test serialized reference survival array and mortality timing
   from manifest rather than only syntheticearlydeath. Default-off identity
   at nonzero lambda_d must remain unchanged; newflagpositivecredit rejects.
6. Add test that changed isolated function is actualmodule imported bynative
   floor path and warn prepared overlay is NOT production runtime acceptance.
   Current one-file overlay lacks full import package; give safe activation
   recipe for fresh process / registeredfunction before runtime imports or a
   future isolated packagecopyonTorch. Don't claim futureproduction integration
   alreadyvalidated. No entire repositorycopyor checkpointcopy.
7. Pin source/test/launcher hashes onTorch launcher before executing. Receipts
   include reference label. Source hashes checks before runningtests unchanged.

## Context

Full AGENTS startup applies, but prior worker alreadyreadcontext. Mathematical
specunchanged: noNEWunsecuredcredit; rolloverfloor beforedeathrisk; b'>=0 at
positive mortality/terminal, exact j-nativeindex. No native naturalcredit mode.
Frozen label/source/checkpoint unchanged. User wants removal, positivecreditopen.

## Do not touch

No sharedsource, imports/numerictests onMac, productionjobs, Git, Googlewrites.
Do not submitTorchjob. Five-minute cap, no furtherdiagnosis/scope expansion.

## Required output

Corrected minimal driver/tests/launcher/README, short finaldiffsummary viawrapper
review_fix_findings.md. Mark tests notrun. Lead reviews then runs one five-minute
Torchzero-solve job, separatelyauthorized notoldfailedaccountingjob.

## Verification

LightweightAST/text only. Lead line-by-line mathreviewandTorch targetedchecks.

## Stop and report if

Cannotfixwithout assumption/newpathsor cap reached. Never claimrun passed.
