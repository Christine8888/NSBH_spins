#!/bin/bash
echo "Running TOVs"
conda activate igwn-py39
python mTOV_slope_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2 --slope 0.2
python mTOV_slope_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --slope 0.2
python mTOV_slope_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --slope 0.2
python mTOV_slope_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2 --slope 0.2

python mTOV_slope_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2  --spin_slope 2 --folder lowspin --slope 0.2
python mTOV_slope_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2  --spin_slope 2 --folder lowspin --slope 0.2
python mTOV_slope_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2  --spin_slope 2 --folder lowspin --slope 0.2
python mTOV_slope_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2  --spin_slope 2 --folder lowspin --slope 0.2

python mTOV_slope_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --slope 0.2
python mTOV_slope_1c.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --slope 0.2
python mTOV_slope_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --slope 0.2
python mTOV_slope_1c.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2  --spin_slope 0 --max_jjkep 0.5 --folder medspin --slope 0.2
