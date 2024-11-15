import matplotlib.pyplot as plt

# Sample data: replace these values with your actual recorded runtimes
cores = [1, 2, 4, 8, 16, 32]  # Number of cores used in each test
runtimes = [26.95 ,13.88, 7.24, 4.09, 2.79, 1.75]  # Corresponding runtimes for each core count

# Calculate the speedup (compared to single-core runtime)
single_core_runtime = runtimes[0]
speedup = [single_core_runtime / runtime for runtime in runtimes]

# Plot the runtime vs. number of cores
plt.figure(figsize=(12, 5))

# Subplot 1: Runtime vs. Number of Cores
plt.subplot(1, 2, 1)
plt.plot(cores, runtimes, marker='o', color='b', label="Runtime")
for x, y in zip(cores, runtimes):
    plt.text(x, y, f"{y:.2f}", ha='center', va='bottom')  # Label each point
plt.xlabel("Number of Cores")
plt.ylabel("Runtime (seconds)")
plt.title("Runtime vs. Number of Cores")
plt.grid(True)
plt.legend()

# Subplot 2: Speedup vs. Number of Cores
plt.subplot(1, 2, 2)
plt.plot(cores, speedup, marker='o', color='g', label="Speedup")
for x, y in zip(cores, speedup):
    plt.text(x, y, f"{y:.2f}", ha='center', va='bottom')  # Label each point
plt.xlabel("Number of Cores")
plt.ylabel("Speedup (relative to 1 core)")
plt.title("Speedup vs. Number of Cores")
plt.grid(True)
plt.legend()

# Show the plot
plt.tight_layout()
plt.show()