import os
import sys


def run_perf(file_out):
    perf_command = "perf stat"
    options = ['cycles', 'cpu-clock', 'task-clock', 'branches', 'branch-misses', 
               'L1-dcache-loads', 'L1-dcache-load-misses', 'L1-dcache-stores',
               'L1-dcache-store-misses', 'branch-loads', 'branch-load-misses', 'cache-misses']
    func = "../kepler"
    func_args = "-f ../eos/eosA -e 1e15"
    command = perf_command + " -e " + " -e ".join(options) + ' ' + func + ' ' + func_args  + "> /dev/null 2>> " + str(file_out)
    os.system(command)

def get_perf_output():
    perf_lines = []
    with open('perf.out', 'r') as f:
        recording = False
        for line in f:
            if "Performance counter stats for" in line:
                recording = True
            if recording:
                perf_lines += line
    print "Huurrrr"
    print perf_lines
    print "done"
    return perf_lines

def print_lines_to_file(lines, file_out):
    with open(file_out, 'a') as f:
        for line in lines:
            f.write(line)

if __name__ == "__main__":
    run_perf(sys.argv[1]) 
   # lines = get_perf_output()
   # print_lines_to_file(lines, sys.argv[1])

