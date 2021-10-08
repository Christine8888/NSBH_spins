#!/bin/bash
echo "Running real data"
conda activate igwn-py39

python run_real_data.py --gw190426 --gw190814 --type direct --name widespin --poptype 2c --freespin
python run_real_data.py --gw190426 --gw190814 --type chieff --name widespin --poptype 2c --freespin
python run_real_data.py --gw190426 --gw190814 --type pos_chieff --name widespin --poptype 2c --freespin

python run_real_data.py --gw190426 --type direct --name widespin --poptype 2c --freespin
python run_real_data.py --gw190426 --type chieff --name widespin --poptype 2c --freespin
python run_real_data.py --gw190426 --type pos_chieff --name widespin --poptype 2c --freespin

python run_real_data.py --gw190814 --type direct --name widespin --poptype 2c --freespin
python run_real_data.py --gw190814 --type chieff --name widespin --poptype 2c --freespin
python run_real_data.py --gw190814 --type pos_chieff --name widespin --poptype 2c --freespin

python run_real_data.py --type direct --name widespin --poptype 2c --freespin
python run_real_data.py --type chieff --name widespin --poptype 2c --freespin
python run_real_data.py --type pos_chieff --name widespin --poptype 2c --freespin


#python run_real_data.py --gw190426 --gw190814 --type direct --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --gw190426 --gw190814 --type chieff --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --gw190426 --gw190814 --type pos_chieff --name widespin --beta 0 --poptype 2c --freespin

#python run_real_data.py --gw190426 --type direct --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --gw190426 --type chieff --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --gw190426 --type pos_chieff --name widespin --beta 0 --poptype 2c --freespin

#python run_real_data.py --gw190814 --type direct --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --gw190814 --type chieff --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --gw190814 --type pos_chieff --name widespin --beta 0 --poptype 2c --freespin

#python run_real_data.py --type direct --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --type chieff --name widespin --beta 0 --poptype 2c --freespin
#python run_real_data.py --type pos_chieff --name widespin --beta 0 --poptype 2c --freespin
