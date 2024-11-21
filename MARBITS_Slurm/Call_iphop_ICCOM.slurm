#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=15-00:00:00
#SBATCH --job-name=iphop_ICCOM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/users/fcoutinho/ICCOM/

# module load perl
# perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/split_fasta.pl /mnt/smart/users/fcoutinho/ICCOM/IMGVR_Seqs_in_ICCOM_Models.fasta 5000

module load iphop/iphop_env

conda activate /mnt/smart/apps/condaEnvs/iphop_env

for i in Part*fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; mkdir iPhop_Results_${smp} ; iphop predict --min_score 75 --num_threads 23 --fa_file $i --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir iPhop_Results_${smp}/; rm -fr iPhop_Results_${smp}/Wdir/ ; done

# mkdir iPhop_Results

# iphop predict --num_threads 47 --fa_file /mnt/smart/users/fcoutinho/ICCOM/IMGVR_Seqs_in_ICCOM_Models.fasta --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir /mnt/smart/users/fcoutinho/ICCOM/iPhop_Results/

#Above command never finishes and gives error:
# Process ForkPoolWorker-559:
# Traceback (most recent call last):
#   File "/mnt/smart/apps/condaEnvs/iphop_env/lib/python3.8/multiprocessing/process.py", line 315, in _bootstrap
#     self.run()
#   File "/mnt/smart/apps/condaEnvs/iphop_env/lib/python3.8/multiprocessing/process.py", line 108, in run
#     self._target(*self._args, **self._kwargs)
#   File "/mnt/smart/apps/condaEnvs/iphop_env/lib/python3.8/multiprocessing/pool.py", line 114, in worker
#     task = get()
#   File "/mnt/smart/apps/condaEnvs/iphop_env/lib/python3.8/multiprocessing/queues.py", line 356, in get
#     res = self._reader.recv_bytes()
#   File "/mnt/smart/apps/condaEnvs/iphop_env/lib/python3.8/multiprocessing/connection.py", line 216, in recv_bytes
#     buf = self._recv_bytes(maxlength)
#   File "/mnt/smart/apps/condaEnvs/iphop_env/lib/python3.8/multiprocessing/connection.py", line 421, in _recv_bytes
#     return self._recv(size)
#   File "/mnt/smart/apps/condaEnvs/iphop_env/lib/python3.8/multiprocessing/connection.py", line 386, in _recv
#     buf.write(chunk)
# MemoryError


#iphop predict --num_threads 23 --fa_file /mnt/smart/users/fcoutinho/ICCOM/IMGVR_Seqs_in_ICCOM_Models.fasta --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir /mnt/smart/users/fcoutinho/ICCOM/iPhop_Results/

