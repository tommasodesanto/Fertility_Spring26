# Codex worker task — urgent GE addon, 8-minute implementation cap

Author now explicitly requests BOTH same-price and GE for the two debt rules.
Absolute lead deadline2026-09-30T01:00:25UTC; do not extend. Fullcontextrequired
but boundedmandatoryread latestprefix/daily/topstatus only, no broadscans/dumps.
Own ONLY credit_rule_ge_quick_v1 under this packet. NO SUBMISSION, git, cleanup,
extraagents or modification of the other worker's credit_rule_quick_v1 sources.
Other worker61540 prepares 2PE cases/strict_tenure.py; its promptquick_compare_prompt.md
fully defines both rules. Read that prompt and reuse its functions/code once
present. It will be ready byabout00:40UTC. No duplicated solver/kernelimplementation.

GE contract EXACTLY existingcompleted credit_ge_v1:
fixed allestimatedpreferencesincludingpsi,entry,income,age/birthtiming,fiscalrule,
absoluteHs(q) curve,H0, native16/20queue. q solvesACTUALB/(2.1E)-1=0 using
sol.adult_entry_adjusted_birth_children/sol.entry_rate, neverpsi adjustment.
N=Hs(q)/d(q) scalesactualdemographicstationaryhouseholds to clearphysical supply.
Retentionofcurve meansH0constant; supplymayrespondvia unchanged .63exponent.
Do not confuse this withfixed physicalhousingstock sensitivity, whichis separate.
Retain actualPAYGO/nonnegativeestates/nativequeue checks, no closure substitutions.
Original160grid, not262naturalcredit. CreditrulesONLYdifference inprompt; owner
LTV/grandfathering staysfrozen. No sourceidentitywaivers or parameterfilesoverlays.

Use existingroot/accountinglogic from parentcredit_ge_v1/run_ge.py, especially
closed_accounting, choose_bracket, next_log_factor, scaled_native_step; maycopy
these reviewedsmallfunctions verbatim, pinthesource. Don'tcallitsnaturalcredit
adapter,261grid,psi initializer, hard-codedupward-onlybracket oroldcasebudget.
Root must supportlowerorhigherq: q0 thenperturb±.03 firstaccordingto residual
sign(positivebirthresidual→raiseq,negative→lowerq). Ifnotbracketed expand±.08,
then±.20max, stopwithnofallback. Safeguardedlogsecant insidebracket.
Exactly2GEcases, parallelarray2each1CPU24GiB, atmost8NEWlifecycle solves perGE
case(includingq0ifcan'treusepassedPEreceipt actualrenewalfields). No exactrepeat
requiredforPRELIMINARY readout; labelunverifiedrepeat. 180sec/perpricesolve,
stopnewsolve if cannotfitabsolute deadline minus60sec readout reserve.
No automaticretries/renormalization/newcredit/entry truncation. RawPE q0 canseed
rootonlyifauthenticatedsource/params/strictmodule identical andactualrenewalinputs
present;otherwiseoneq0solve permissiblewithin8 cap. Reportstartdeadlineonce.

Write eachcandidate summaryatonce, latest/best/progress, currentbracket/failure.
Ifselected actualrenewal|residual|<=1e-6, report price,rent, N, tenure,birthtiming,
wealth,C,actualfiscal/solvency/pop-scaledclearing. All14fit/31paramtables and17
standardplots remainrequired; send preliminarysmallnumbers BEFORE slowrendering.
Do NOTpresent lastcandidate asclearedGE iffiscal/entry/solvency/renewalgatefails.
Strictzero initialnegativeentrant wealth mustremain untouched; report quantified
infeasibilityifany rather thanchangeinput ortry priceswhenincome debt cannotclear.
CheckpointstayTorchandatomicsaveifneeded; no largeMac download.

Deliver driver/launcher/plan/shortREADME, exactcommand/interface/dependencies,
smallsyntheticrootchecks ifpossible (localpurealgebraallowed), no submission.
Ifotherworker's interface unavailable, define narrow documented callableAPI
around its authenticate/applyrule/solve andstop ratherthaninvent economicfields.
Leadwillverifyroot math/contracts beforelaunch. Remote root
/scratch/td2248/projects/fixed_reference_credit_rule_ge_quick_20260929.
