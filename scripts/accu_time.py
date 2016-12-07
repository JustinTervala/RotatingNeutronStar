#!/usr/bin/env python

import os
import sys
import timers

if __name__ == "__main__":
    num = int(sys.argv[1])
    base_command = '../kepler -f ../eos/eosA -e 1e15 '
    command = base_command + " > out.out "
    for x in range(num-1):
        command += ' && ' + base_command + " >> out.out"
    print "Executing " + str(num) +  " iterations of kepler"
    os.system(command)
    timers.main()
