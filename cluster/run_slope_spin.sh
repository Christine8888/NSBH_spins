#!/bin/bash
echo "Running TOVs"
conda activate igwn-py39
python mTOV_slope_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2 --folder mTOV_convergence_99pc --slope 0.2 --freespin
python mTOV_slope_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2  --spin_slope 2 --folder lowspin_99pc --slope 0.2 --freespin
python mTOV_slope_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2  --spin_slope 0 --max_jjkep 0.5 --folder medspin_99pc --slope 0.2 --freespin
