#!/bin/bash
echo "Running TOVs"
conda activate igwn-py39
python mTOV_convergence.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2
python mTOV_convergence.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2
#python mTOV_convergence_2component.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2
#python mTOV_convergence_2component.py --detector "APlus" --event_min 30 --event_max 150 --n_events 5 --mtov_true 2.2 
