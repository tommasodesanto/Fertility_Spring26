#!/bin/bash
#SBATCH --job-name=credit-ge-contact-sheets
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:02:00
#SBATCH --output=/scratch/td2248/projects/fixed_reference_credit_ge_20260929/verification_v1/contact_sheets/render_%j.log
set -euo pipefail
module load anaconda3/2025.06
python /scratch/td2248/projects/fixed_reference_credit_ge_20260929/verification_v1/render_contact_sheets.py
