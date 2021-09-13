#!/bin/bash
echo "Running LMG"
conda activate igwn-py39
#python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 5 --spin_slope 0
#python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.2 --spin_slope 0
#python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.65 --spin_slope 0
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 5 --spin_slope 2 --folder LMG_convergence_lowspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.2 --spin_slope 2 --folder LMG_convergence_lowspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.65 --spin_slope 2 --folder LMG_convergence_lowspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 5 --spin_slope 0 --max_jjkep 0.5 --folder LMG_convergence_medspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.2 --spin_slope 0 --max_jjkep 0.5 --folder LMG_convergence_medspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.65 --spin_slope 0 --max_jjkep 0.5 --folder LMG_convergence_medspin
