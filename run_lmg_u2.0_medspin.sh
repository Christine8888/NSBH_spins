#!/bin/bash

python LMG_u.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin
#python LMG_u.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.5 --bh_min 2.5  --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin
python LMG_u.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.0 --bh_min 2.0  --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin
python LMG_u.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.0 --bh_min 2.41  --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin
#python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.5 --bh_min 2.5 --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.0 --bh_min 2.0 --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin
python LMG_u.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.0 --bh_min 2.41 --spin_slope 0 --max-jjkep 0.5 --folder LMG_convergence_medspin