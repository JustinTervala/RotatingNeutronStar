import numpy as np
import sys
import matplotlib.pyplot as plt

def read_data():
    with open('out.out', 'r') as f:
        spin_time = {x:[] for x in ['total', 'rho_gamma', 'pn', 'vep', 'metric', 'ang', 'rad', 'coeff', 'alpha']}
        mass_radius_time = []
        for line in f:
            if "spin():" in line:
                spin_time['total'].append(float(line.split(':')[1]))
            elif "spin()" in line:	
                timer = line.split(',')[-1]	
                if "rho_gamma:" in timer:
                    spin_time['rho_gamma'].append(float(timer.split(':')[1]))
                elif "pn:" in timer:
                    spin_time['pn'].append(float(timer.split(':')[1]))
                elif "vep:" in timer:
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
    print "spin():"
    print "Func".center(9) + "Avg".center(9) + "Min".center(9) + "Max".center(9) + "Std".center(9) + "Calls".ljust(7) + "Percent".ljust(7)

def print_stats_line(name, data, total):
    stddev = round(np.std(data)/1.e6, 2)
    avg = round(np.average(data)/1.e6,2)
    min_x = round(min(data)/1.e6, 2)
    max_x = round(max(data)/1.e6, 2)
    print name.rjust(9) + str(avg).center(9) + str(min_x).center(9) \
        + str(max_x).center(9) + str(stddev).center(9) + str(len(data)).ljust(7) + str(round(np.sum(data)/total*100.0,3)).ljust(7)

def print_spin_table(data):
    total = np.sum(data['total'])
    for x in sorted(data, key=lambda k: np.sum(data[k]), reverse=True):
        print_stats_line(x, data[x], total) 

def printdata(name, data):
    print "  " + name + ": " + str(data) 

def print_stats(name, data):
    stddev = np.std(data)/1.e6
    avg = np.average(data)/1.e6
    print name + "():"
    printdata("min", min(data)/1.e6)
    printdata("max", max(data)/1.e6)
    printdata("avg", avg)
    printdata("std", stddev)

def main():
    spin, mr = read_data()
    print "Execution Time in Milliseconds"
    print_spin_header()
    print_spin_table(spin)
    print_stats("mass_radius", mr)
    return spin

if __name__ == "__main__":
    spin = main()
    #if len(sys.argv) == 2:
        #plt.plot(spin[sys.argv[1]]
        #plt.show()    
