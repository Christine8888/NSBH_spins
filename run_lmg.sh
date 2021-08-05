#!/bin/bash
echo "Running TOVs"
conda activate igwn-py39
python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 0
python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 3 --bh_min 3  --spin_slope 0
python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.5 --bh_min 3.01  --spin_slope 0
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 0
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 3 --bh_min 3 --spin_slope 0
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.5 --bh_min 3.01 --spin_slope 0

#python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 2
#python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 3 --bh_min 3  --spin_slope 2
#python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.5 --bh_min 3.01  --spin_slope 2
#python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 2
#python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 3 --bh_min 3 --spin_slope 2
#python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.5 --bh_min 3.01 --spin_slope 2
