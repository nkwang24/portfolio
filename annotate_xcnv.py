# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:16:15 2018

@author: nkwang

Annotate XHMM CNV caller output with gene and transcript information
"""

from pyensembl import Genome
import csv
import re
import time

def annotate_xcnv():

    with open(r'Z:\gatk\data\DATA.xcnv','r') as f1:
        DATA_xcnv = list(csv.reader(f1, delimiter='\t'))
        
    data = Genome(
            reference_name='GRCh37',
            annotation_name='GRCh37',
            gtf_path_or_url=r'Z:\gatk\reference\Homo_sapiens.GRCh37.75.gtf')
    data.index()
    
    with open(r'./VPv1_exons_positive','r') as f2:
        exons_positive = list(csv.reader(f2, delimiter='\t'))
    with open(r'./VPv1_exons_negative_flipped','r') as f3:
        exons_negative = list(csv.reader(f3, delimiter='\t'))
    
    for index, row in enumerate(DATA_xcnv[1:]):
        chr = str(row[4])
        chr_start = int(re.split(':|-', row[2])[1])+1
        CNV_start = data.gene_names_at_locus(contig=chr, position=int(chr_start))
        try:
            chr_end = int(re.split(':|-', row[2])[2])
        except IndexError:
            chr_end = chr_start
        CNV_end = data.gene_names_at_locus(contig=chr, position=chr_end)
        
        for row in exons_positive:
            if chr == row[0][3:]:
                if chr_start >= int(row[1]) and chr_start <= int(row[2]):
                    match = row[3].split('_')
                    CNV_start.append(match[0] + '_' + match[1] + ' ' + match[2] + ' ' + str(int(match[3])+1))
                if chr_end >= int(row[1]) and chr_end <= int(row[2]):
                    match = row[3].split('_')
                    CNV_end.append(match[0] + '_' + match[1] + ' ' + match[2] + ' ' + str(int(match[3])+1))
                
        for row in exons_negative:
            if chr == row[0][3:]:
                if chr_start >= int(row[1]) and chr_start <= int(row[2]):
                    match = row[3].split('_')
                    CNV_start.append(match[0] + '_' + match[1] + ' ' + match[2] + ' ' + match[3])
                if chr_end >= int(row[1]) and chr_end <= int(row[2]):
                    match = row[3].split('_')
                    CNV_end.append(match[0] + '_' + match[1] + ' ' + match[2] + ' ' + match[3])
        
        DATA_xcnv[index+1].append(CNV_start)
        DATA_xcnv[index+1].append(CNV_end)
   
    DATA_xcnv[0] += ['CNV_START', 'CNV_END']
    
    with open(r'Z:\gatk\data\DATA_annotated.xcnv','w',newline='') as fo:
        writer = csv.writer(fo, delimiter='\t')
        for row in DATA_xcnv:
            writer.writerow(row)
    
    return
    
    
def main():
    
    start_time = time.time()
    annotate_xcnv()
    end_time = time.time()    
    print('Completed ' + str(round(end_time - start_time, 2)))
    
    return
    
    
if __name__ == '__main__':
    main()
    
    
    