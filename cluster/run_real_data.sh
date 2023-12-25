#!/bin/bash
echo "Running real data"
conda activate igwn-py39

python run_real_data.py --gw190426 --gw190814 --type direct --name widespin --freespin --beta 0
python run_real_data.py --gw190426 --gw190814 --type chieff --name widespin --freespin --beta 0
python run_real_data.py --gw190426 --gw190814 --type pos_chieff --name widespin --freespin --beta 0

python run_real_data.py --gw190426 --type direct --name widespin --freespin --beta 0
python run_real_data.py --gw190426 --type chieff --name widespin --freespin --beta 0
python run_real_data.py --gw190426 --type pos_chieff --name widespin --freespin --beta 0

python run_real_data.py --gw190814 --type direct --name widespin --freespin --beta 0
python run_real_data.py --gw190814 --type chieff --name widespin --freespin --beta 0
python run_real_data.py --gw190814 --type pos_chieff --name widespin --freespin --beta 0

python run_real_data.py --type direct --name widespin --freespin --beta 0
python run_real_data.py --type chieff --name widespin --freespin --beta 0
python run_real_data.py --type pos_chieff --name widespin --freespin --beta 0