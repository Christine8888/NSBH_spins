#!/bin/bash
echo "Running TOVs"
conda activate igwn-py39
python mTOV_slope_2c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --slope 0.2 --free --folder mTOV_convergence
python mTOV_slope_2c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --slope 0.2 --free --max_jjkep 0.5 --folder medspin
python mTOV_slope_2c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --slope 0.4 --free --max_jjkep 0.5 --folder medspin
python mTOV_slope_2c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2 --slope 0.4 --free --max_jjkep 0.5 --folder medspin
python mTOV_slope_2c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2 --slope 0.4 --free --spin_slope 2 --folder lowspin
