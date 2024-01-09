import matplotlib.pyplot as plt
import numpy as np
import os.path as p
import json

thread_numbers = [1, 2, 5, 10, 20, 40, 60]

json_path = "times.json"
time_json = {thread: 0.0 for thread in thread_numbers}


def get_array_for_size(size: int):
    out = []
    for thread in thread_numbers:
        out.append(time_json[str(size)][str(thread)])
    return out


def plot_all():
    pass
    """for size in sizes_to_time:
        times = get_array_for_size(size)
        plt.plot(thread_numbers, times[0]/np.array(times), label=f"{size}")

    plt.title("Speedup for different sizes and threads")
    plt.xlabel("Number of threads")
    plt.ylabel("Speedup")
    plt.legend(loc="upper left")
    plt.show()"""


if __name__ == '__main__':
    if p.exists(json_path):
        with open(json_path, 'r') as file:
            time_json = json.load(file)
    else:
        print("no json you idiot run timing.py first")

    plot_all()
