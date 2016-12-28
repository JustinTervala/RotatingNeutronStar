#!/usr/bin/env python

import run_perf
import accu_time
import os
import sys

def both(iterations, file_out):
    accu_time.exec_and_write(int(iterations), file_out, True)
    run_perf.run_perf(file_out)


if __name__ == "__main__":
    both(sys.argv[1], sys.argv[2])
