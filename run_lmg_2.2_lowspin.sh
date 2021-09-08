#!/bin/bash
echo "Running LMG"
conda activate igwn-py39

#python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2 --bh_min 5 --spin_slope 2 --folder LMG_convergence_lowspin
#python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2 --bh_min 2.2  --spin_slope 2 --folder LMG_convergence_lowspin
#python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.2 --bh_min 2.65  --spin_slope 2 --folder LMG_convergence_lowspin
#python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 5 --spin_slope 2 --folder LMG_convergence_lowspin
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.2 --spin_slope 2 --folder LMG_convergence_lowspin
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 --bh_min 2.65 --spin_slope 2 --folder LMG_convergence_lowspin
