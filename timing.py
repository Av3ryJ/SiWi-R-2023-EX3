import subprocess
import os.path as p
import json

binary = "./cg"
thread_numbers = [1, 2, 5, 10, 20, 40, 60]
number_of_iterations = 100
eps = -1
nxy = 4096

json_path = "times.json"
time_json = {thread: 0.0 for thread in thread_numbers}


def run_bin(threads):
    result = subprocess.run(["mpirun", "-np", threads, binary, nxy, nxy, number_of_iterations, eps], capture_output=True, text=True)
    return result.stdout


def time_all():
    # iterate over sizes and threads
    for thread_num in thread_numbers:
        print(f"running: {thread_num}")
        returned = run_bin(str(thread_num))
        time = float(returned)
        print(f"time was: {time}")
        time_json[thread_num] = time
    with open(json_path, 'w') as f:
        json.dump(time_json, f)


if __name__ == '__main__':
    if p.exists(json_path):
        with open(json_path, 'r') as file:
            time_json = json.load(file)
    else:
        time_all()
