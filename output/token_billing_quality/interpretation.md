# Interpretation

In 39,687 non-tie Arena battles, the lower-preferred answer uses more visible output tokens in 37.7% of matched tasks. The mean loser-minus-winner token difference is -55.1, while the median is -34.0.

The fixed-effect LPM regresses an indicator that response A wins on A's log token advantage over response B, prompt length, and model A/B fixed effects. The coefficient is 0.087 with standard error 0.003. This is descriptive, not causal: the exercise establishes that visible billable usage and observed preference are not mechanically aligned in matched tasks.
