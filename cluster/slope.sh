#!/bin/bash

python run_real_data.py --slope --gw190426 --gw190814 --type chieff --name slope --poptype 1c --freespin --beta 0
python run_real_data.py --slope --gw190426 --type chieff --name slope --poptype 1c --freespin --beta 0
python run_real_data.py --slope --gw190814 --type chieff --name slope --poptype 1c --freespin --beta 0
python run_real_data.py --slope --type chieff --name slope --poptype 1c --freespin --beta 0

python run_real_data.py --slope --gw190426 --gw190814 --type chieff --name slope --poptype u --freespin --beta 0
python run_real_data.py --slope --gw190426 --type chieff --name slope --poptype u --freespin --beta 0
python run_real_data.py --slope --gw190814 --type chieff --name slope --poptype u --freespin --beta 0
python run_real_data.py --slope --type chieff --name slope --poptype u --freespin --beta 0

python run_real_data.py --slope --gw190426 --gw190814 --type chieff --name slope --poptype 2c --freespin --beta 0
python run_real_data.py --slope --gw190426 --type chieff --name slope --poptype 2c --freespin --beta 0
python run_real_data.py --slope --gw190814 --type chieff --name slope --poptype 2c --freespin --beta 0
python run_real_data.py --slope --type chieff --name slope --poptype 2c --freespin --beta 0
