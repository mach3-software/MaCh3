#!/bin/bash

# Remove any existing log files
rm -f memory_usage.log memory_usage_data.txt

dnf install -y procps

pip install matplotlib

# Source environment setup scripts
source bin/setup.MaCh3.sh
source bin/setup.MaCh3Tutorial.sh
source bin/setup.NuOscillator.sh

# Start the executable in the background
./bin/MCMCTutorial Inputs/FitterConfig.yaml &

# Get the PID of the process
PID=$!

# Monitor memory usage while the process is running
while kill -0 $PID 2>/dev/null; do
    ps -o pid,rss,command -p $PID >> memory_usage.log
    sleep 0.1
done

# Extract the memory usage (RSS values) from the log file
awk 'NR>1 {print $2}' memory_usage.log | tail -n +2 > memory_usage_data.txt

# Inform that the process is completed
echo "Process completed" >> memory_usage.log

# Run the Python script to plot the memory usage
python3 - <<EOF
import matplotlib.pyplot as plt

# Read the memory usage data
memory_data = []
with open('memory_usage_data.txt', 'r') as file:
    for line in file.readlines():
        try:
            # Convert RSS from KB to MB (1 MB = 1024 KB)
            memory_data.append(int(line.strip()) / 1024)
        except ValueError:
            continue

# Create a time axis (in seconds) corresponding to each memory reading
time_data = [i * 0.1 for i in range(len(memory_data))]

# Set the figure size and resolution (larger figure, higher dpi)
plt.figure(figsize=(12, 8), dpi=150)  # figsize: (width, height) in inches, dpi for resolution

# Plot the data
plt.plot(time_data, memory_data)
plt.xlabel('Time (seconds)')
plt.ylabel('RAM Usage (RSS in MB)')
plt.title('Memory Usage Over Time')
plt.grid(True)

# Save the plot as a PNG file
plt.savefig('memory_usage_plot.png')

# Close the plot
plt.close()
EOF

echo "Memory usage plot saved as memory_usage_plot.png"
