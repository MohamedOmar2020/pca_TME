#! /bin/bash -l
#SBATCH --partition=panda-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name=jup
#SBATCH --time=2-00:00:00
#SBATCH --mem=350G
#SBATCH --gres=gpu:1


source ~/.bashrc
conda activate /home/mao4005/.conda/envs/pathml

jupyter notebook --no-browser --ip 0.0.0.0 --port=8958
