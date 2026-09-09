# Verified terminal-root receipts

Both no-rebate terminal stationary roots pass the independent receipt review. This verifies endpoints, not a completed historical perfect-foresight path or recalibration.

The review passed 143 local artifact/arithmetic/gate checks and 3034 remote source/input/data hash checks (1535 distinct files). Only JSON/CSV outputs were collected. Large checkpoint/state files remained remote and were checked by streaming hashes; their arrays were not independently deserialized or recomputed in this review.

| Quantity | Sequential | Fertility nests |
|---|---:|---:|
| Stationary asset price | 0.464205470709 | 0.463910356505 |
| Signed housing-market residual | 4.61891104375e-06 | -1.99117584978e-05 |
| Household heads, model units | 0.679257111332 | 0.679365385888 |
| Resident persons, model units | 1.6487991616 | 1.64908925197 |
| Annual births per head | 0.0172545658279 | 0.0172566795054 |
| Tax revenue less transfers, period model units | 0.0997786880761 | 0.0995658920475 |
| Root runtime, seconds | 384.322599289 | 365.376163472 |

Each arm completed eight evaluations. Evaluation eight is an uncached repeat of evaluation seven at the same price, with an exactly identical housing residual. Both satisfy the 2e-4 market gate and every recorded endpoint gate. The largest person one-step relative error is 1.298e-11; household/person head-mass discrepancies are below 1.289e-12.

The fiscal contract is retained 1% annual property tax (0.04 per four-year period), zero transfer. The positive fiscal surplus is intentional under fixed transfer; these results do not establish an equal-rebate fiscal equilibrium.

The roots preserve each normalized 2007 supply anchor with elasticity 0.63. Each arm derives its old fertility intercept to match 2.1, then holds its corresponding normalized 2023 intercept fixed in the tail. Terminal survival, migration and headship inputs are frozen as declared. These are conditional future assumptions. Differences across arms include their separate derived fertility normalization and supply anchoring; they are not a decomposition of the choice formula alone.

## Verified remote checkpoint and packet hashes

Hashes below refer to the saved output contract.json, not its separately referenced launch contract.

### sequential

- `summary.json`: `2cb954e2c66479140f0050a5f0c17666fef4f6423d0e025cb2a46fc9ba37778b`
- `contract.json`: `99a6478b042157e5681fe62f9a6b0b0b3c0a99a4c6bf65d589e0978ca911fda0`
- `terminal.pkl.gz`: `4bc7d631eefc92713ecea119e0d79ee050578191a40b4bcc55212cfbe54b24b5`

### nested

- `summary.json`: `645ead721206ddec344a3b53a15ef76313b07a982ff7babdee6176a79a0cf7c5`
- `contract.json`: `8e8b740322933448fbbc53a95dbe4ab75abd02fa1629e4386d17dd5c093633f8`
- `terminal.pkl.gz`: `a3305c655d948c4f0016e03dbe9829125f17121b3b8857118027f9b8cd6e65ef`

Full checks: `verification.json`. Remote artifact and explicit pin manifest: `remote_hash_verification.json`. Stage receipts: `endpoint_01/`, `initial_01/`, `root_01/`.
