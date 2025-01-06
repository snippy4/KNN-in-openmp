import matplotlib.pyplot as plt

# Sample data: replace these values with your actual recorded runtimes
cores = [1, 2, 4, 8, 16, 32]  # Number of cores used in each test
runtimes = [19.385, 10.972, 10.704, 14.388, 22.522, 33.540]  # Corresponding runtimes for each core count

# Calculate the speedup (compared to single-core runtime)
single_core_runtime = runtimes[0]
speedup = [single_core_runtime / runtime for runtime in runtimes]

# Ideal speedup is equal to the number of cores
ideal_speedup = cores

# Plot the speedup against the ideal speedup
plt.figure(figsize=(10, 6))

# Actual Speedup
plt.plot(cores, speedup, marker='o', color='g', label="Actual Speedup")
for x, y in zip(cores, speedup):
    plt.text(x, y, f"{y:.2f}", ha='center', va='bottom')  # Label each point

# Ideal Speedup
plt.plot(cores, ideal_speedup, linestyle='--', color='r', label="Ideal Speedup")
for x, y in zip(cores, ideal_speedup):
    plt.text(x, y, f"{y}", ha='center', va='bottom')  # Label each point

# Labels, Title, and Legend
plt.xlabel("Number of Cores")
plt.ylabel("Speedup")
plt.title("Speedup vs. Ideal Speedup")
plt.grid(True)
plt.legend()

# Show the plot
plt.tight_layout()
plt.show()
