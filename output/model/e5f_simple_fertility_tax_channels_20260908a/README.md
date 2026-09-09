# Rebated property-tax channel decomposition

Eight fresh fixed-price household solves at the selected calibration. Bits follow tax, asset price, equal rebate. Cells 000 and 111 reproduce the verified rebated 1% and 2% endpoints, including all policy arrays. The six mixed cells require both replays to pass; they retain household accounting gates but intentionally need not clear housing or balance government budgets.

Reuse the existing exact eight-cell Shapley decomposition. Report births, rooms, ownership and family-group allocation. No calibration, future transition, new figures or production promotion. Full calibration tables remain in the overnight morning review.

Budget: eight solves, two endpoint cells then six mixed cells, at roughly 20–60 seconds each based on observed replay timing; allow several minutes per wave plus queue time. Each case has a 30-minute watchdog and 35-minute Slurm cap, one CPU and 24GiB. No retries or gate relaxation. Stop descendants if smoke fails. Heartbeats every 30 seconds; completed-case summaries, fixed selected-best summary and checkpoints saved per cell. Collector verifies all receipts and exact component add-up. No ongoing monitoring automation.

Source: isolated codex/fertility-nest-computation. Scientific bundle unchanged. Contract SHA256: `de767353697c706dec4a7d960d1f3314dfb5cf2debad4c17e8438ce021966288`.

Invocation is frozen in cells.sbatch and collect.sbatch. Submit cells 0,7 first, cells 1-6 with afterok dependency, then collector afterok. Submission receipt records actual IDs.

Completed: all eight cells and collector pass. See [channel results](results/READOUT.md), with all 17 metric decompositions in the companion CSV and JSON. Lead independently verified report hashes and component add-up.
