#Script that splits the files in the folder into batches of n files each as specified in the command args
import glob
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--batch_size", help="number of files in each sub directory", type =int, default=1000)
parser.add_argument("--input_dir", help="input directory with files to split", type =str)
parser.add_argument("--extension", help="extension of files in the input directory to split", type =str, default="")
args = parser.parse_args()

in_files = glob.glob(f'{args.input_dir}/*{args.extension}')

print(f'Number of files in the input directory: {len(in_files)}')

batch_num = 0
for (i, file) in enumerate(in_files):
    if (i % args.batch_size == 0):
        print(f'Batch {batch_num} done')
        batch_num += 1
        os.makedirs(f'Batch_{batch_num}', exist_ok=True)
    os.rename(file, f'Batch_{batch_num}/{os.path.basename(file)}')
