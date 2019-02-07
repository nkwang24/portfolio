# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 09:41:06 2018

@author: nkwang

Process coverage gap reports from targeted panel sequencing
to determine common gaps for additional bait design and rare
gaps to identify potential homozygous deletions
"""

import os
import numpy as np
import csv
import time

training_directory = os.fsencode('./Training Set/')
report_directory = os.fsencode('./Reports to Flag/')

def process_gaps(chr, report):
    with open(report, 'r') as f:
        gaps = list(csv.reader(f,delimiter='\t'))[5:]
    
    report_gaps = np.zeros(shape=(1,249250621),dtype=np.uint8)
    
    for gap in gaps:
        if gap[2] == str(chr):
            report_gaps[0][int(gap[3])-1:int(gap[4])] = 1
    
    return(report_gaps)


def integrate_gaps(chr):
    global training_directory
    
    print('Chr:',chr)
    cumulative_gaps = np.zeros(shape=(1,249250621),dtype=np.uint8)
    training_n = 0
    
    for report in os.listdir(training_directory):
        if report.endswith(b'.txt'):
    #        print(report)
            training_n += 1
            cumulative_gaps += process_gaps(chr, training_directory + report)
    
    print('Cumulative gap length (N=' + str(training_n) + ') = ' + str(np.sum(cumulative_gaps)) + '\n')
    
    return(cumulative_gaps) #number of gaps by chr and chr pos
    

def identify_common(cumulative_gaps):
    
    
    
    return


def flag_uncommon(report, cumulative_gaps, training_n, chr):
    global report_directory
    
    print('\t' + os.fsdecode(report))
    with open(report_directory + report, 'r') as f:
        readCSV = list(csv.reader(f))
    
    total = 0
    flagged = 0
    
    if chr == 1:
        filteredCSV = readCSV[:5]
    else:
        filteredCSV = []
    
    for i in range(5,len(readCSV)):
        gap_start = int(readCSV[i][3])-1
        gap_end = int(readCSV[i][4])
        gap_length = gap_end - gap_start
        
        avg_gap_freq = float(sum(cumulative_gaps[0][gap_start:gap_end]))/gap_length/training_n
#        min_gap_freq = float(min(cumulative_gaps[0][gap_start:gap_end]))/training_n
        
        if readCSV[i][2] == str(chr):
            total += 1
#            print(cumulative_gaps[0][gap_start:gap_end])
#            print(avg_gap_freq)
#            print(min_gap_freq)
        
        if readCSV[i][2] == str(chr) and avg_gap_freq <= 0.05 and gap_length > 100:# or min_gap_freq <= 0.1):
            flagged += 1
            readCSV[i] += ['', 'Rare', round(avg_gap_freq, 2)]#, round(min_gap_freq, 2)]
            filteredCSV.append(readCSV[i])
        
    with open(report_directory + report, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(readCSV)
        
    with open(report_directory + report[:-4] + b'_Filtered.csv', 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(filteredCSV)
         
    print('\tFlagged gaps in Chr' + str(chr) + '= ' + str(flagged) + '/' + str(total) + '\n')
    
    return


def main():
    start_time = time.time()
    
    global training_directory
    global report_directory
    training_n = len(os.listdir(training_directory))
    
    for report in os.listdir(report_directory):
        if report.endswith(b'.txt'):
            with open(report_directory + report, 'r') as f:
                readCSV = list(csv.reader(f,delimiter='\t'))
            readCSV[2] += ['avg_gap_freq <= 0.05 and gap_length > 100']
            readCSV[4] += ['Gap Flag', 'Average Gap Frequency']#, 'Minimum Gap Frequency']
            with open(report_directory + report[:-18] + b'.csv', 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(readCSV)
        if report.endswith(b'Coverage.csv'):
                with open(report_directory + report[:-4] + b'_Filtered.csv', 'w', newline='') as f:
                    pass
    
    for chr in list(range(1,23)) + ['X']:
        cumulative_gaps = integrate_gaps(chr)
        identify_common(cumulative_gaps)
        for report in os.listdir(report_directory):
            if report.endswith(b'Coverage.csv'):
                flag_uncommon(report, cumulative_gaps, training_n, chr)
    
    print('Complete')
    
    end_time = time.time()
    print('Run time = ' + str(round(end_time - start_time, 2)) + ' seconds')
    
    return

    
if __name__ == '__main__':
    main()
    
    