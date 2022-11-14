#!/usr/bin/env python

# For benchmarking metaSPAdes
import os
import pandas as pd
import psutil as ps
import time
from subprocess import PIPE, check_output

# For plotting the results of the benchmark
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

##########################
# Run metaSPAdes benchmark
##########################

# List of 21 samples for the HTP dataset
samples = ['SRR14530767','SRR14530766','SRR14530765','SRR14530764','SRR14530763','SRR14530762','SRR14530881','SRR14530880','SRR14530772','SRR14530771','SRR14530770','SRR14530769','SRR14530888','SRR14530890','SRR14530889','SRR14530887','SRR14530886','SRR14530885','SRR14530884','SRR14530882','SRR14530891']

# List of number of samples to iterate through
n_samples = [1, 5, 10, 15, 20]

# Create a dataframe to store the results
df = pd.DataFrame(columns=['n_samples', 'runtime', 'peak_mem', 'peak_cpu', 'avg_mem', 'avg_cpu', 'sum_gb_size', 'n_reads'])

# Define paths
workdir = '/nobackup1/phil_p/STRONG/STRONG_Runs/RothmanUnenrichedHTPBenchmarkResults'
fastq_dir = '/nobackup1c/users/phil_p/STRONG/STRONG_Runs/RothmanUnenriched/'
bind_dirs = '/nobackup1/phil_p/STRONG,/nobackup1c/users/phil_p/STRONG'
img = '/home/phil_p/.singularity/strong-1.0.simg'
out_path = f'{workdir}/spades_benchmark.csv'
out_fig = f'{workdir}/spades_benchmark.png'

# Run SPAdes - log the time and resources it takes to run
for n in n_samples:

    samples_subset = samples[:n]

    mems = []
    cpus = []
    start_time = time.time()

    fastq1_paths = [f'{fastq_dir}/{sample}/{sample}_1.fastq.gz' for sample in samples_subset]
    fastq2_paths = [f'{fastq_dir}/{sample}/{sample}_2.fastq.gz' for sample in samples_subset]
    fastq1_inputs = ' -1 ' + ' -1 '.join(fastq1_paths)
    fastq2_inputs = ' -2 ' + ' -2 '.join(fastq2_paths)
    spades_cmd = f'/STRONG/SPAdes/assembler/bin/spades.py --meta --only-assembler --save-gp -t 20 -k 77 -m 200 {fastq1_inputs} {fastq2_inputs} -o assembly/spades/ &> assembly/spades_{n}.log'
    singularity_cmd = f'singularity exec --home {workdir} --bind {bind_dirs} {img} bash -c \'set -euo pipefail; {spades_cmd}\''
    # singularity_cmd = f'echo "singularity exec --home {workdir} --bind /nobackup1/phil_p/STRONG,/nobackup1c/users/phil_p/STRONG /home/phil_p/.singularity/strong-1.0.simg bash -c \'set -euo pipefail; {spades_cmd}\'" && sleep 10'

    # Create the process
    print(f'Running command {singularity_cmd}')
    process = ps.Popen(singularity_cmd, stdout=PIPE, shell=True)

    # while the process is running calculate resource utilization.
    while(process.is_running()):
        # set the sleep time to monitor at an interval of every second.
        time.sleep(5)

        # capture the memory and cpu utilization at an instance
        mems.append((ps.virtual_memory().total - ps.virtual_memory().available) / (float)(2**30))
        cpus.append(ps.cpu_percent())

        if float(process.memory_info().rss) == 0.0:
            break
    
    # Track the peak utilization of the process
    peak_mem = max(mems)
    peak_cpu = max(cpus)
    avg_mem = sum(mems)/len(mems)
    avg_cpu = sum(cpus)/len(cpus)
    total_mem = ps.virtual_memory().total / (float)(2**30)
    # Record time for the process to finish
    run_time = time.time() - start_time

    # Print the results to the monitor for each subprocess run.
    print(f'Time to run SPAdes with {n} samples: {run_time} seconds')
    print(f'Peak memory usage is {peak_mem} of {total_mem} GB')
    print(f'Peak CPU utilization is {peak_cpu}')

    # Get the size of the sum of the fastq paths and in GB and the total number of reads
    sum_gb_size = 0
    n_reads = 0
    for sample in samples_subset:
        sum_gb_size += sum(os.path.getsize(f'{fastq_dir}/{sample}/{sample}_{i}.fastq.gz') for i in [1, 2]) / 1000000000
        print(f'Total size of fastq files: {sum_gb_size} GB')
        n_reads += sum(int(check_output(f'gunzip -c {fastq_dir}/{sample}/{sample}_{i}.fastq.gz | wc -l', shell=True)) for i in [1, 2]) / 4
        print(f'Total number of reads: {n_reads}')
    
    # Aggregate all of the logged data into a single dataframe
    df = df.append({'n_samples': n, 'runtime': run_time, 'peak_mem': peak_mem, 'peak_cpu': peak_cpu, 'avg_mem': avg_mem, 'avg_cpu':avg_cpu, 'sum_gb_size': sum_gb_size, 'n_reads': n_reads}, ignore_index=True)

    # Write the dataframe to CSV
    df.to_csv(out_path, index=False)

############################
# Plot the benchmark results
############################

df = pd.read_csv(out_path)

# Fix the data types
df['n_samples'] = df['n_samples'].astype(int)
df['runtime'] = df['runtime'].astype(float)
df['peak_mem'] = df['peak_mem'].astype(float)
df['peak_cpu'] = df['peak_cpu'].astype(int)
df['sum_gb_size'] = df['sum_gb_size'].astype(float)
df['n_reads'] = df['n_reads'].astype(int)
# Convert runtime to minutes
df['runtime'] = df['runtime'] / 60

# Define vars for plotting
xs = ['n_samples', 'sum_gb_size', 'n_reads']
x_labels = ['Number of Samples', 'Total Size of FASTQ Files (GB)', 'Total Number of Reads']
ys = ['runtime', 'peak_mem', 'peak_cpu'] # , 'avg_mem', 'avg_cpu'
y_labels = ['Runtime (minutes)', 'Peak Memory Usage (GB)', 'Peak CPU Utilization (%)'] # , 'Average Memory Usage (GB)', 'Average CPU Utilization (%)'

sns.set(font_scale=1.5)
sns.set_style('whitegrid')

# Create a grid of subplots for all combinations of X and Y
fig, axs = plt.subplots(len(xs), len(ys), figsize=(20, 20))
# Create scatter plots for each combination of X and Y with a linear regression line
for i, x in enumerate(xs):
    for j, y in enumerate(ys):
        sns.regplot(x=x, y=y, data=df, ax=axs[i, j])
        axs[i, j].scatter(df[x], df[y])
        axs[i, j].set(ylim=(0, max(df[y]+10)))
        axs[i, j].set_xlabel(x_labels[i])
        axs[i, j].set_ylabel(y_labels[j])
        # Calculate and annotate Rp
        r, p = stats.pearsonr(df[x], df[y])
        axs[i, j].annotate(f'$R_\\rho = {r:.2f}$', xy=(.5, .4), xycoords='axes fraction', fontsize=20) # horizontalalignment='right', verticalalignment='bottom'
        fig.tight_layout()

# Save the figure
fig.savefig(out_fig, dpi=300)