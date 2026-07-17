# Codex workers

`scripts/codex-worker.sh` is a small, explicit launcher for delegated
`codex exec` tasks. It gives the caller a reproducible profile choice and
keeps the default sandbox read-only. `--write` is required before a worker can
modify the workspace; the wrapper also uses `--ephemeral`, rejects workdirs
outside this repository, and appends a no-commit/no-push instruction.

The built-in `spawn_agent` facility does not expose model selection, so it
cannot guarantee a lower-usage worker. This wrapper supplies an explicit model
and reasoning effort to `codex exec` for every named profile in
`config/models.env`.

## Profiles

| Profile | Model | Reasoning | Context | Default access / limit |
| --- | --- | --- | --- | --- |
| `explorer_fast` | GPT-5.6-Terra | low | minimal | read-only, 10 min |
| `mechanic_fast` | GPT-5.3-Codex-Spark | low | minimal | read-only, 10 min |
| `worker_fast` | GPT-5.6-Terra | medium | full | read-only, 30 min |
| `reviewer_strong` | GPT-5.6-Terra | high | full | read-only, 20 min |

Minimal profiles start from an empty temporary directory and receive the
selected repository workdir through `--add-dir`. This prevents Codex from
automatically loading the root instructions and calibration memory stack for
purely mechanical tasks. The wrapper still supplies the absolute repository
and workdir paths in its safety instructions. Minimal context is not permitted
for economics, calibration, targets, model results, or model-code
interpretation. Use `--full-context` to override a minimal profile when a
seemingly simple task needs that context.

## Usage

Start from the task template and fill every section:

```bash
cp ops/codex-workers/prompts/task-template.md /tmp/my-task.md
ops/codex-workers/scripts/codex-worker.sh \
  --profile explorer_fast \
  --prompt-file /tmp/my-task.md
```

Use `--write` only for a task explicitly authorized to edit the repository.
Optionally save the worker's final message inside the repository:

```bash
ops/codex-workers/scripts/codex-worker.sh \
  --profile mechanic_fast \
  --prompt-file /tmp/my-task.md \
  --write \
  --workdir code/model \
  --final-message output/model/worker-final.md
```

Run `ops/codex-workers/scripts/codex-worker.sh --help` for the full interface.
Unknown profiles and missing prompt files fail before Codex is invoked. The
wrapper never invokes a Git commit or push. Each profile has a wall-clock limit;
use `--time-limit-min N` only when the task itself justifies an override. Use `--dry-run` to inspect the
resolved model, reasoning effort, context and launch modes, sandbox, and
workdir without spending a model call.
