#!/bin/bash

python run_real_data.py --gw190426 --gw190814 --type fake --name slope --poptype 1c --freespin --beta 0 --slope
python run_real_data.py --gw190426 --type fake --name slope --poptype 1c --freespin --beta 0 --slope
python run_real_data.py --gw190814 --type fake --name slope --poptype 1c --freespin --beta 0 --slope
python run_real_data.py --type fake --name slope --poptype 1c --freespin --beta 0 --slope

python run_real_data.py --gw190426 --gw190814 --type fake --name slope --poptype u --freespin --beta 0 --slope
python run_real_data.py --gw190426 --type fake --name slope --poptype u --freespin --beta 0 --slope
python run_real_data.py --gw190814 --type fake --name slope --poptype u --freespin --beta 0 --slope
python run_real_data.py --type fake --name slope --poptype u --freespin --beta 0 --slope

python run_real_data.py --gw190426 --gw190814 --type fake --name slope --poptype 2c --freespin --beta 0 --slope
python run_real_data.py --gw190426 --type fake --name slope --poptype 2c --freespin --beta 0 --slope
python run_real_data.py --gw190814 --type fake --name slope --poptype 2c --freespin --beta 0 --slope
python run_real_data.py --type fake --name slope --poptype 2c --freespin --beta 0 --slope
