import matplotlib.pyplot as plt

# Sample data: replace these values with your actual recorded runtimes
cores = [1, 2, 4, 8, 16, 32]  # Number of cores used in each test
runtimes = [143, 73, 72, 59, 29, 32]  # Corresponding runtimes for each core count

# Calculate the speedup (compared to single-core runtime)
single_core_runtime = runtimes[0]
speedup = [single_core_runtime / runtime for runtime in runtimes]

# Calculate parallel efficiency
efficiency = [s / c for s, c in zip(speedup, cores)]

# Plot the runtime, speedup, and efficiency
plt.figure(figsize=(15, 7))

# Subplot 1: Runtime vs. Number of Cores
plt.subplot(1, 3, 1)
plt.plot(cores, runtimes, marker='o', color='b', label="Runtime")
for x, y in zip(cores, runtimes):
    plt.text(x, y, f"{y:.2f}", ha='center', va='bottom')  # Label each point
plt.xlabel("Number of Cores")
plt.ylabel("Runtime (seconds)")
plt.title("Runtime vs. Number of Cores")
plt.grid(True)
plt.legend()

# Subplot 2: Speedup vs. Number of Cores
plt.subplot(1, 3, 2)
plt.plot(cores, speedup, marker='o', color='g', label="Speedup")
for x, y in zip(cores, speedup):
    plt.text(x, y, f"{y:.2f}", ha='center', va='bottom')  # Label each point
plt.xlabel("Number of Cores")
plt.ylabel("Speedup (relative to 1 core)")
plt.title("Speedup vs. Number of Cores")
plt.grid(True)
plt.legend()

# Subplot 3: Efficiency vs. Number of Cores
plt.subplot(1, 3, 3)
plt.plot(cores, efficiency, marker='o', color='r', label="Efficiency")
for x, y in zip(cores, efficiency):
    plt.text(x, y, f"{y:.2f}", ha='center', va='bottom')  # Label each point
plt.xlabel("Number of Cores")
plt.ylabel("Efficiency")
plt.title("Parallel Efficiency vs. Number of Cores")
plt.grid(True)
plt.legend()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
