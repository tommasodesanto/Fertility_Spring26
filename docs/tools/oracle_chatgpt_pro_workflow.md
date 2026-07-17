# Oracle and ChatGPT Pro Workflow

This note records the local workflow for using Oracle to package repository
context and relaying that package to a ChatGPT Pro conversation through the
Codex Chrome plugin / Brave extension path.

## Current Setup

- Oracle skill: `~/.codex/skills/oracle`
- Oracle CLI: `oracle` (`@steipete/oracle`)
- Browser route: Codex Chrome plugin, using the Brave/Chromium extension
  backend when available
- Preferred ChatGPT target for the toy model relay:
  `OLG Model for Fertility`
  (`https://chatgpt.com/c/6a1f50c9-203c-83ea-8e5d-6878dd32e815`)

The Chrome plugin may be cached locally even when the extension backend is not
active in a new Codex session. Before relying on browser automation, check that
the extension backend is available and can see the relevant ChatGPT tab.

## New Codex Chat Bootstrap Prompt

Use this when starting a new Codex chat that should relay with ChatGPT Pro:

```text
We are in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

First, follow the repository startup protocol:
1. Read `memory/AGENT_MEMORY.md`
2. Read the latest `memory/daily/YYYY-MM-DD.md`
3. Read `CALIBRATION_STATUS.md`
4. Only read `SESSION_DIARY.md` if historical context is needed

Task: become a local Codex relay to my ChatGPT Pro conversation about a toy OLG
fertility model.

Use the Codex Chrome plugin / extension backend through Brave if available. The
target ChatGPT conversation is pinned and titled:

`OLG Model for Fertility`
URL: `https://chatgpt.com/c/6a1f50c9-203c-83ea-8e5d-6878dd32e815`

Important:
- Do not inspect unrelated ChatGPT conversations.
- Open or claim only the target pinned OLG conversation.
- The ChatGPT conversation should use Pro / Pro Extended if available.
- First retrieve the latest completed response in that conversation if I have
  already asked it for an orientation summary.
- Report the summary back to me here.
- Then wait for my instructions. Going forward, relay carefully packaged
  prompts from this local repo/session to that ChatGPT Pro conversation, bring
  back its answers, and help me decide what to implement or ask next.
- Do not send broad repository dumps unless I explicitly approve. Use narrow
  context bundles, preferably via Oracle when code/file context is needed.
- Treat ChatGPT Pro's answers as advisory: verify claims against local files,
  equations, and tests before implementing anything.
```

## Oracle Usage

For any file-context request, first preview the bundle:

```bash
oracle --dry-run summary --files-report \
  -p "<prompt for ChatGPT Pro>" \
  --file "AGENTS.md" \
  --file "CALIBRATION_STATUS.md" \
  --file "path/to/relevant/file.py"
```

If the file set is tight and contains no sensitive material, render the bundle:

```bash
oracle --render --copy \
  -p "<prompt for ChatGPT Pro>" \
  --file "AGENTS.md" \
  --file "CALIBRATION_STATUS.md" \
  --file "path/to/relevant/file.py"
```

If `--copy` fails because `pbcopy` is unavailable inside the sandbox, render to
a reviewed temporary file and paste that through the browser bridge or manually:

```bash
oracle --render --render-plain \
  -p "<prompt for ChatGPT Pro>" \
  --file "AGENTS.md" \
  --file "CALIBRATION_STATUS.md" \
  --file "path/to/relevant/file.py" \
  > /private/tmp/oracle_bundle.md
```

## Safety Rules

- Use the smallest file set that contains the relevant truth.
- Do not send `.env` files, keys, tokens, private logs, raw data, browser
  history, or broad generated output.
- Do not send the whole repository.
- Prefer canonical status files and active model files over archived outputs.
- For calibration/model claims, include target values next to model moments
  when available.
- Treat Pro's output as an external review, not as source of truth.
- Verify any proposed code, equation, calibration interpretation, or empirical
  claim locally before acting on it.

## Quick Browser Check

In a fresh Codex chat, the browser automation path is healthy only if the
Chrome/Brave extension backend appears in the available browser list and can
see open ChatGPT tabs. A successful minimal test is to open a fresh ChatGPT tab,
send `hello world`, and confirm ChatGPT replies.

If the backend is unavailable:

1. Confirm the Codex Chrome plugin is enabled in Codex.
2. Confirm the Codex Chrome Extension is installed and enabled in Brave.
3. Restart Codex if the plugin was just enabled.
4. Retry the extension-backend check before attempting Oracle handoff.
