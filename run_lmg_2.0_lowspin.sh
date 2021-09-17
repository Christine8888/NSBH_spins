#!/bin/bash
echo "Running LMG"
conda activate igwn-py39

python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 2 --folder LMG_convergence_lowspin_99pc --free
python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.0 --bh_min 2.0  --spin_slope 2 --folder LMG_convergence_lowspin_99pc --free
python LMG.py --detector "Design" --event_min 10 --event_max 50 --n_events 5 --mtov_true 2.0 --bh_min 2.41  --spin_slope 2 --folder LMG_convergence_lowspin_99pc --free
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2 --bh_min 5 --spin_slope 2 --folder LMG_convergence_lowspin_99pc --free
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.0 --bh_min 2.0 --spin_slope 2 --folder LMG_convergence_lowspin_99pc --free
python LMG.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.0 --bh_min 2.41 --spin_slope 2 --folder LMG_convergence_lowspin_99pc --free
