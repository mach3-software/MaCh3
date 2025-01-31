#!/bin/bash

# Remove any existing log files
rm -f memory_usage.log memory_usage_data.txt cpu_usage_data.txt

# Install necessary tools if not already installed
dnf install -y procps
pip install matplotlib

# Source environment setup scripts
source bin/setup.MaCh3.sh
source bin/setup.MaCh3Tutorial.sh
source bin/setup.NuOscillator.sh

# Get the number of CPU cores to have normalised CPU usage
NUM_CORES=$(nproc)

# Start the executable in the background
./bin/MCMCTutorial Inputs/FitterConfig.yaml General:MCMC:NSteps:100000 &

# Get the PID of the process
PID=$!

# Monitor memory and CPU usage while the process is running
while kill -0 $PID 2>/dev/null; do
    # Log memory usage (RSS in KB) and CPU usage (%CPU)
    ps -o pid,rss,%cpu,command -p $PID --no-headers >> memory_usage.log
    sleep 0.1
done

# Extract memory usage (RSS values) and CPU usage from the log file
awk '{print $2}' memory_usage.log > memory_usage_data.txt
awk '{print $3}' memory_usage.log > cpu_usage_data.txt

# Inform that the process is completed
echo "Process completed" >> memory_usage.log

# Run the Python script to plot RAM and CPU usage
python3 - <<EOF
import matplotlib.pyplot as plt

# Number of CPU cores
num_cores = $NUM_CORES

# Read the memory usage data
memory_data = []
cpu_data = []
with open('memory_usage_data.txt', 'r') as mem_file:
    for line in mem_file:
        try:
            memory_data.append(int(line.strip()) / 1024)  # Convert KB to MB
        except ValueError:
            continue

with open('cpu_usage_data.txt', 'r') as cpu_file:
    for line in cpu_file:
        try:
            cpu_data.append(float(line.strip()) / num_cores)  # Normalize CPU usage
        except ValueError:
            continue

# Calculate mean and max for RAM and CPU usage
mean_ram = sum(memory_data) / len(memory_data) if memory_data else 0
max_ram = max(memory_data) if memory_data else 0
mean_cpu = sum(cpu_data) / len(cpu_data) if cpu_data else 0
max_cpu = max(cpu_data) if cpu_data else 0

# Save the summary table in markdown format
with open('summary_table.txt', 'w') as summary_file:
    summary_file.write(f"### Memory and CPU Usage Summary\n\n")
    summary_file.write(f"| Metric  |     Mean    |      Max    |\n")
    summary_file.write(f"|---------|-------------|-------------|\n")
    summary_file.write(f"| **RAM** | {mean_ram:8.2f} MB | {max_ram:8.2f} MB |\n")
    summary_file.write(f"| **CPU** | {mean_cpu:8.2f} %  | {max_cpu:8.2f} %  |\n")

# Create a time axis (in seconds) corresponding to each memory reading
time_data = [i * 0.1 for i in range(len(memory_data))]

# Set the figure size and resolution (larger figure, higher dpi)
plt.figure(figsize=(12, 8), dpi=150)

# Plot RAM usage
plt.subplot(2, 1, 1)
plt.plot(time_data, memory_data, label='RAM Usage (MB)', color='blue')
plt.xlabel('Time (seconds)')
plt.ylabel('RAM (MB)')
plt.title('Memory Usage Over Time')
plt.grid(True)
plt.legend()

# Plot CPU usage
plt.subplot(2, 1, 2)
plt.plot(time_data, cpu_data, label='CPU Usage (%)', color='red')
plt.xlabel('Time (seconds)')
plt.ylabel('CPU (%)')
plt.title('CPU Usage Over Time')
plt.grid(True)
plt.legend()

# Save the plot as a PDF file
plt.tight_layout()
plt.savefig('ram_cpu_usage.pdf')
plt.close()

print("Summary table and plots saved.")
EOF

echo "Summary table and plots saved."
