#!/bin/bash
echo "Running real data"
conda activate igwn-py39

python run_real_data.py --gw190426 --gw190814 --type direct --name wide --poptype 2c
python run_real_data.py --gw190426 --gw190814 --type chieff --name wide --poptype 2c
python run_real_data.py --gw190426 --gw190814 --type pos_chieff --name wide --poptype 2c

python run_real_data.py --gw190426 --type direct --name wide --poptype 2c
python run_real_data.py --gw190426 --type chieff --name wide --poptype 2c
python run_real_data.py --gw190426 --type pos_chieff --name wide --poptype 2c

python run_real_data.py --gw190814 --type direct --name wide --poptype 2c
python run_real_data.py --gw190814 --type chieff --name wide --poptype 2c
python run_real_data.py --gw190814 --type pos_chieff --name wide --poptype 2c

python run_real_data.py --type direct --name wide --poptype 2c
python run_real_data.py --type chieff --name wide --poptype 2c
python run_real_data.py --type pos_chieff --name wide --poptype 2c


#python run_real_data.py --gw190426 --gw190814 --type direct --name wide --beta 0 --poptype 2c
#python run_real_data.py --gw190426 --gw190814 --type chieff --name wide --beta 0 --poptype 2c
#python run_real_data.py --gw190426 --gw190814 --type pos_chieff --name wide --beta 0 --poptype 2c

#python run_real_data.py --gw190426 --type direct --name wide --beta 0 --poptype 2c
#python run_real_data.py --gw190426 --type chieff --name wide --beta 0 --poptype 2c
#python run_real_data.py --gw190426 --type pos_chieff --name wide --beta 0 --poptype 2c

#python run_real_data.py --gw190814 --type direct --name wide --beta 0 --poptype 2c
#python run_real_data.py --gw190814 --type chieff --name wide --beta 0 --poptype 2c
#python run_real_data.py --gw190814 --type pos_chieff --name wide --beta 0 --poptype 2c

#python run_real_data.py --type direct --name wide --beta 0 --poptype 2c
#python run_real_data.py --type chieff --name wide --beta 0 --poptype 2c
#python run_real_data.py --type pos_chieff --name wide --beta 0 --poptype 2c
