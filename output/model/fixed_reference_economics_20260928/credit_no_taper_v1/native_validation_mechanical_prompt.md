# Codex worker task

## Goal / context
Context-minimal MECHANICAL harness repair, no economics/model interpretation.
The lead has already specified and reviewed the scientific contract. Do not
read huge memory/status/model sources. Profile mechanic_fast, five-minute cap.

## Scope / ownership
Own ONLY native_validation_v6 within this packet. Clone exactly driver/launcher/
plan/README from native_validation_v5, preserve v5. No cleanup, git, other edits,
SSH, numerical submission, scientific imports/tests or extra workers.

## Two precise fixes
1. run_native_validation.py function mock_tests nested run(name,...,**kw) passes
sleeper=lambda... explicitly and forwards **kw. Timeout caller supplies sleeper
again. Fix this ordinary Python keyword duplication by selecting the default
sleeper only when kw does not supply one, then forward exactly once. This fixes
Torch preflight18836944 TypeError, before any model work.
2. install_identity_guard_extension currently adds only ROOT/code/model/tools
to sys.path. Also add ROOT/code/model BEFORE its model-module import, so its
existing actual solver import in self-test is resolvable. Do not change any
guard, source SHA constant, equations, function contents beyond these fixes.

Update driver SHA in plan and its pinned-files entry. Launcher source/result
directory versions become v6, but immutable scientific schema MUST REMAIN
block0506_renter_no_taper_native_validation_v3 in BOTH plan and checks. Do not
globally replace that schema while changing folder paths. Keep all budgets,
cases, source pins and economics unchanged. Original driver will change only
the two narrowly identified sites. README briefly records cause; no report.

## Verification / deliverable
Static syntax/linkage only on Mac; no NumPy/model import. Print concise diff,
new SHA and files. Confirm plan-schema literal matches driver and launcher,
updated SHA matches file, two lifecycle maximum/360case/1200total unchanged.
Do not submit. Stop if more than mechanical source edits are needed.
