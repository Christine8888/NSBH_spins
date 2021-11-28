#!/bin/bash

python run_real_data.py --gw190426 --gw190814 --type fake --name fake --poptype 1c --freespin --beta 3
python run_real_data.py --gw190426 --type fake --name fake --poptype 1c --freespin --beta 3
python run_real_data.py --gw190814 --type fake --name fake --poptype 1c --freespin --beta 3
#python run_real_data.py --type fake --name fake --poptype 1c --freespin --beta 3

python run_real_data.py --gw190426 --gw190814 --type fake --name fake --poptype u --freespin --beta 3
python run_real_data.py --gw190426 --type fake --name fake --poptype u --freespin --beta 3
python run_real_data.py --gw190814 --type fake --name fake --poptype u --freespin --beta 3
#python run_real_data.py --type fake --name fake --poptype u --freespin --beta 3

python run_real_data.py --gw190426 --gw190814 --type fake --name fake --poptype 2c --freespin --beta 3
python run_real_data.py --gw190426 --type fake --name fake --poptype 2c --freespin --beta 3
python run_real_data.py --gw190814 --type fake --name fake --poptype 2c --freespin --beta 3
#python run_real_data.py --type fake --name fake --poptype 2c --freespin --beta 3
