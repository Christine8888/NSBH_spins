#!/bin/bash
echo "Running real data"
conda activate igwn-py39

python run_real_data.py --gw190426 --gw190814 --type direct --name default
python run_real_data.py --gw190426 --gw190814 --type chieff --name default
python run_real_data.py --gw190426 --gw190814 --type pos_chieff --name default
