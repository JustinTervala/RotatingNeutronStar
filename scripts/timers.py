import numpy as np
import sys

def read_data():
    with open('out.out', 'r') as f:
        spin_time = {x:[] for x in ['total', 'vep', 'metric', 'ang', 'rad', 'coeff', 'alpha']}
        mass_radius_time = []
        for line in f:
            if "spin():" in line:
                spin_time['total'].append(float(line.split(':')[1]))
            elif "spin()" in line:	
                timer = line.split(',')[-1]	
                if "vep:" in timer:
                    spin_time['vep'].append(float(timer.split(':')[1]))
                elif "metric:" in timer:
                    spin_time['metric'].append(float(timer.split(':')[1]))
                elif "ang:" in timer:
                    spin_time['ang'].append(float(timer.split(':')[1]))
                elif "rad:" in timer:
                    spin_time['rad'].append(float(timer.split(':')[1]))
                elif "coeff:" in timer:
                    spin_time['coeff'].append(float(timer.split(':')[1]))
                elif "alpha:" in timer:
                    spin_time['alpha'].append(float(timer.split(':')[1]))
            elif "mass_radius()" in line:
                mass_radius_time.append(float(line.split(':')[1]))

    return spin_time, mass_radius_time

def print_spin_header():
    spin_header = []
    spin_header += ["spin():"]
    spin_header += ["Func".center(9) + "Avg".center(9) + "Med".center(9) + "Min".center(9) + "Max".center(9) + "Std".center(9) + "Calls".ljust(7) + "Percent".ljust(7)]
    return spin_header

def print_stats_line(name, data, total):
    stddev = round(np.std(data)/1.e6, 2)
    avg = round(np.average(data)/1.e6,2)
    med = round(np.median(data)/1.e6,2)
    min_x = round(min(data)/1.e6, 2)
    max_x = round(max(data)/1.e6, 2)
    return [name.rjust(9) + str(avg).center(9) + str(med).center(9) + str(min_x).center(9) \
        + str(max_x).center(9) + str(stddev).center(9) + str(len(data)).ljust(7) + str(round(np.sum(data)/total*100.0,3)).ljust(7)]

def print_spin_table(data):
    spin_table_lines = []
    total = np.sum(data['total'])
    for x in sorted(data, key=lambda k: np.sum(data[k]), reverse=True):
        spin_table_lines += print_stats_line(x, data[x], total) 
    return spin_table_lines

def printdata(name, data):
    return ["  " + name + ": " + str(data)]

def print_stats(name, data):
    stddev = np.std(data)/1.e6
    avg = np.average(data)/1.e6
    lines = []
    lines += [name + "():"]
    lines += printdata("min", min(data)/1.e6)
    lines += printdata("max", max(data)/1.e6)
    lines += printdata("avg", avg)
    lines += printdata("std", stddev)
    return lines

def read_timer_data():
    spin, mr = read_data()
    lines = []
    lines += ["Execution Time in Milliseconds"]
    lines += print_spin_header()
    lines += print_spin_table(spin)
    lines += print_stats("mass_radius", mr)
    return spin, lines

if __name__ == "__main__":
    spin, ll = read_timer_data()
    for line in lines:
        print line
