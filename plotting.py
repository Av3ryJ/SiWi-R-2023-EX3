import matplotlib.pyplot as plt
import numpy as np
import os.path as p
import json

thread_numbers = [1, 2, 5, 10, 20, 40, 60, 80]

json_path = "times.json"
time_json = {thread: 0.0 for thread in thread_numbers}


def get_array():
    out = []
    for thread in thread_numbers:
        out.append(time_json[str(thread)])
    return out


def plot_all():
    times = get_array()
    plt.plot(thread_numbers, times[0]/(np.array(times)*np.array(thread_numbers)))

    plt.title("Parallel efficiency for n processes")
    plt.xlabel("Number of processes")
    plt.ylabel("parallel efficiency")
    plt.show()

if __name__ == '__main__':
    if p.exists(json_path):
        with open(json_path, 'r') as file:
            time_json = json.load(file)
    else:
        print("no json you idiot run timing.py first")

    plot_all()
