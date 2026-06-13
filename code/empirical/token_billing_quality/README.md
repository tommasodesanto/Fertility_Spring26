# Token Billing and Quality Proof of Concept

This folder contains a small public-data proof of concept for the note
`latex/token_billing_quality_minipaper.tex`.

The script uses the Hugging Face dataset
`lmarena-ai/arena-human-preference-55k`. Each row is a matched Chatbot Arena
battle with the same prompt, two model responses, model names, and an observed
preference label. The analysis counts visible output tokens under one common
`tiktoken` tokenizer. It does not attempt to reconstruct provider-specific
billing meters.

Run from the repository root:

```bash
python3 -m venv /tmp/token_billing_poc_venv
/tmp/token_billing_poc_venv/bin/python -m pip install -r code/empirical/token_billing_quality/requirements.txt
/tmp/token_billing_poc_venv/bin/python code/empirical/token_billing_quality/arena_token_billing_poc.py
```

Main outputs are written to `output/token_billing_quality/`:

- `paper_table.tex`: compact results table included in the mini-paper.
- `loser_token_premium_hist.pdf`: figure included in the mini-paper.
- `arena_summary_table.csv` and `arena_regression_table.csv`: machine-readable
  results.
- `interpretation.md`: short generated interpretation of the estimates.
