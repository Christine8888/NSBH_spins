#!/bin/bash
python mTOV_bias_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --free
python mTOV_bias_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --free
python mTOV_bias_u.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --free
python mTOV_bias_u.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --free
