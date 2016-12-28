#!/usr/bin/env python

import os
import sys
import timers

def multi_exec(iterations):
    base_command = '../kepler -f ../eos/eosA -e 1e15 '
    command = base_command + " > out.out "
    for x in range(iterations-1):
        command += ' && ' + base_command + " >> out.out"
    print "Executing " + str(iterations) +  " iterations of kepler"
    os.system(command)
    return timers.read_timer_data()

def write_to_file(lines, num, file_out):
    with open(file_out, 'w') as f:
        f.write("Timing for %d iterations of kepler\n" % num)
        for line in lines:
            f.write("%s\n" % line) 

def exec_and_write(num, file_out, print_out=False):
    spin, lines = multi_exec(num)
    write_to_file(lines, num,  file_out)
    if print_out:
        for line in lines:
            print line

if __name__ == "__main__":
    num = int(sys.argv[1])
    spin, lines = multi_exec(num)
    if len(sys.argv) > 2:
        write_to_file(lines, num, sys.argv[2])
    for line in lines:
        print line
