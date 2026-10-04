"""Record actual scratch free space and the cumulative retained-output budget."""
import json
import os
import time
from pathlib import Path

STAGE = Path('/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1')


def main():
    space = os.statvfs(STAGE)
    free = space.f_bavail * space.f_frsize
    receipt = dict(status='staging_storage_observed', observed_epoch=time.time(),
                   free_bytes=free, free_GiB=free / 1024**3,
                   prior_ten_chains_max_GiB=5000, new_five_chains_max_GiB=2500,
                   combined_retained_case_planning_max_GiB=7500,
                   controller_submission_free_floor_GiB=7800,
                   runtime_shared_free_floor_GiB=350,
                   no_pruning=True)
    (STAGE / 'staging_storage.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
