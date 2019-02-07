# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 10:18:38 2017

@author: nkwang

From sample variant list, select variants for Sanger confirmation
and generate primer map for downstream workflows
"""

import requests
import re
import os
import csv
from bs4 import BeautifulSoup
import time
import xlrd

sample_list = []

def process_report(filename, directory):
    sample_number = filename[0:7]
    print('\n' + sample_number)
    
    report_path = './' + os.fsdecode(directory) + '/' + filename
    sheet = xlrd.open_workbook(report_path).sheet_by_index(0)
    readCSV = []
    
    for row in range(sheet.nrows):
        readCSV.append(sheet.row_values(row))
    
#    with open('./' + os.fsdecode(directory) + '/' + filename, 'r') as f:
#        readCSV = list(csv.reader(f))
    
    row = 0
    file_header = []
    variant_data = []
    
    while readCSV[row][0][0] == '#':
        file_header.append(readCSV[row])
#        row_variant_start = row
        if row < len(readCSV)-1:
            row += 1
        else:
            break
    
    while readCSV[row][0]:
        variant_data.append(readCSV[row])
#        row_variant_end = row
        if row < len(readCSV)-1:
            row += 1
        else:
            break
      
    column_headers = file_header[-1]
    gene_index = 0
    HGVSc_index = column_headers.index('HGVSCoding')
    HGVSp_index = column_headers.index('HGVSProtein')
    zygosity_index = column_headers.index('Zygosity')
    try:
        patho_index = column_headers.index('PredictedPathogenicity')
    except ValueError:
        patho_index = column_headers.index('Pathogenicity')
    chr_index = column_headers.index('Chr:ChrPos')
    CDS_index = column_headers.index('ExonNumber')
#    dbSNPAF_index = column_headers.index('dbSNPAF')
#    variant_frequency_index = column_headers.index('VariantFrequency')
    
    global sample_list
    
#    flag_confirm = 0
    
    for variant in variant_data:
        if variant[patho_index + 1].lower().strip() == 'confirm':
#            flag_confirm = 1
            
            variant_info = [sample_number, variant[gene_index], variant[HGVSc_index],
                            variant[HGVSp_index], variant[zygosity_index]]
            print(', '.join(variant_info[1:3]))
            primer_info = find_primer(variant[gene_index], variant[chr_index])
            if primer_info:
                with open('./Output/sanger_map.csv', 'a', newline='') as sanger_map:
                    writer = csv.writer(sanger_map)
                    sample_name = variant_info[0] + ' ' + primer_info[0]
                    sample_list.append(sample_name)
                    writer.writerow(variant_info + primer_info)
            else:
                with open('./Output/sanger_map.csv', 'a', newline='') as sanger_map:
                    writer = csv.writer(sanger_map)
                    sample_name = variant_info[0]
                    sample_list.append(sample_name)
                    writer.writerow(variant_info)
                with open('./Output/primers_to_design.csv', 'a', newline='') as regions:
                    writer = csv.writer(regions)
                    writer.writerow([filename, variant[gene_index], variant[HGVSc_index],
                                     variant[HGVSp_index], variant[chr_index], variant[CDS_index]])
                
#    select_2 = 0
#    enumerate_variants = 0
    
#    if flag_confirm == 0:
#        print('Negative')
#        while select_2 < 2:
#            try:
#                variant = variant_data[enumerate_variants]
#            except IndexError:
#                print("Couldn't find 2 to confirm")
#                return
#            if variant[patho_index].lower() == 'benign' and float(variant[variant_frequency_index]) < 0.65\
#                                                        and float(variant[variant_frequency_index]) > 0.4:
#                if variant[dbSNPAF_index] == 'N/A':
#                    variant_info = [sample_number, variant[gene_index], variant[HGVSc_index],
#                                    variant[HGVSp_index], variant[zygosity_index]]
#                    print(', '.join(variant_info[1:3]))
#                    primer_info = find_primer(variant[gene_index], variant[chr_index])
#                    if primer_info:
#                        select_2 += 1
#                        well_counter += 1
#                        with open('./Output/sanger_map.csv', 'a', newline='') as sanger_map:
#                            writer = csv.writer(sanger_map)
#                            sample_name = variant_info[0] + ' ' + primer_info[0]
#                            sample_list.append(sample_name)
#                            writer.writerow([enumerate_well(well_counter)] + variant_info
#                                             + primer_info + [sample_name])                            
#                else:
#                    try:
#                        if float(variant[dbSNPAF_index]) < 0.05:
#                            variant_info = [sample_number, variant[gene_index], variant[HGVSc_index],
#                                            variant[HGVSp_index], variant[zygosity_index]]
#                            print(', '.join(variant_info[1:3]))
#                            primer_info = find_primer(variant[gene_index], variant[chr_index])
#                            if primer_info:
#                                select_2 += 1
#                                well_counter += 1
#                                with open('./Output/sanger_map.csv', 'a', newline='') as sanger_map:
#                                    writer = csv.writer(sanger_map)
#                                    sample_name = variant_info[0] + ' ' + primer_info[0]
#                                    sample_list.append(sample_name)
#                                    writer.writerow([enumerate_well(well_counter)] + variant_info
#                                                     + primer_info + [sample_name]) 
#                    except ValueError:
#                        pass
#                    
#            enumerate_variants += 1
    return
    

def enumerate_well(well_counter):
    if well_counter > 96:
        return('~~~')
    column = well_counter%12
    row_index = 'ABCDEFGH'
    row = row_index[well_counter//8]
    well = row + str(column).zfill(2)
    return(well)


def find_primer(gene, chr_location):
    primer_match = []
    
    with open('./primer_database.csv', 'r') as f:
        invCSV = list(csv.reader(f))
    
#    for file in os.listdir('./'):
#        if file.endswith('.csv'):
#            with open('./' + file, 'r') as f:
#                invCSV = list(csv.reader(f))
            
    
    gene_index = invCSV[0].index('Gene')
    forward_index = invCSV[0].index('Primer Forward')
    reverse_index = invCSV[0].index('Primer Reverse')
    
    pattern1 = re.compile(r':[\d]+')
    pattern2 = re.compile(r':[\d]+-')
    pattern3 = re.compile(r'-[\d]+')
    
    for primer in invCSV:
        has_position = 0
        if primer[gene_index] == gene and primer[forward_index] and primer[reverse_index]:
            try:
                chr_range = primer[14]
                has_position = 1
            except IndexError:
#                print('\t' + primer[4])
                chr_range = in_silico_pcr(primer[forward_index], primer[reverse_index])
            
            if has_position == 0:
                primer.append(chr_range)
                with open('./primer_database.csv', 'w', newline='') as f:
                    writeCSV = csv.writer(f)
                    writeCSV.writerows(invCSV)
                    
            if 'chr' in chr_range:
                target = pattern1.search(chr_location).group()[1:]
                range_low = pattern2.search(chr_range).group()[1:-1]
                range_high = pattern3.search(chr_range).group()[1:]

                if target > range_low and target < range_high:
                    print('\t', primer[4], primer[6], primer[7], primer[8])
                    primer_match = [primer[4], primer[6], primer[7], primer[8]]
                    break
                
    return(primer_match)

    
def in_silico_pcr(forward, reverse):

    session = requests.Session()
    
    if forward[0:30] == 'CCATCTCATCCCTGCGTGTCTCCGACTCAG':
        forward = forward[30:]
    if reverse[0:23] == 'CCTCTCTATGGGCAGTCGGTGAT':
        reverse = reverse[23:]

    product_size = '50000'
    
    url = ('https://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=650326111_OUAsAbu6RHjC6yoRfPNso6DtsaIy' '&'
           'org=Human' '&'
           'db=hg19' '&'
           'wp_target=genome' '&'
           'wp_f=' + forward + '&'
           'wp_r=' + reverse + '&'
           'Submit=submit' '&'
           'wp_size=' + product_size + '&'
           'wp_perfect=15' '&'
           'wp_good=15' '&'
           'boolshad.wp_flipReverse=0')
    
    r = session.get(url, stream=True)
    soup = BeautifulSoup(r.text,'html.parser')
    
    tags = soup.find_all('pre')
    if tags:
        chr_pattern = re.compile(r'chr[XxMm\d:+-]+')
        chr_range = chr_pattern.search(str(tags[0])).group()
    else:
        chr_range = 'No Match'
    
    return(chr_range)


#def generate_map():
#    global sample_list
#    with open('./Output/sanger_map.csv', 'a', newline='') as sanger_map:
#        writer = csv.writer(sanger_map)
#        writer.writerows([[],[],[]])
#        row_index = 'ABCDEFGH'
#        sample_counter = 0
#        for row in row_index:
#            row_label = []
#            row_sample = sample_list[sample_counter:(sample_counter + 12)]
#            sample_counter += 12
#            for column in range(12):
#                row_label.append(row + str(column + 1).zfill(2))
#            writer.writerow(row_label)
#            writer.writerow(row_sample)
#    
#    return
    
    
def main():
    start_time = time.time()
    
    directory = os.fsencode('./Reports to Process')
    
    if not os.path.exists('./Output'):
        os.makedirs('./Output')
    with open('./Output/sanger_map.csv', 'w', newline='') as sanger_map:
        writer = csv.writer(sanger_map)
        writer.writerow(['Sample #', 'Gene', 'HGVS Coding', 'HGVS Protein', 'Zygosity',
                         'Primer Name', 'Rack', 'Box', 'Well'])
    
    with open('./Output/primers_to_design.csv', 'w', newline='') as regions:
        writer = csv.writer(regions)
        writer.writerow(['Sample #', 'Gene', 'HGVScoding', 'HGVSprotein', 'Chr Pos', 'CDS #'])
    
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('.xlsx'): 
            process_report(filename, directory)
    
    print('\nComplete')
    end_time = time.time()
    print('Run time = ' + str(round(end_time - start_time, 2)) + ' seconds')
    
    input()
    
    return


if __name__ == '__main__':
    main()
    

