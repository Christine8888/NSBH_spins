#!/bin/bash
echo "Running TOVs"
conda activate igwn-py39
python mTOV_bias_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2 --free --folder mTOV_convergence_99pc
python mTOV_bias_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --free --folder mTOV_convergence_99pc
python mTOV_bias_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --free --folder mTOV_convergence_99pc
python mTOV_bias_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2 --free --folder mTOV_convergence_99pc
