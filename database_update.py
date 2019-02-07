# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 09:47:36 2017

@author: nkwang

Generate SQL file for updating database variant pathogenicity
from analyzed variant tables
"""

import re
import os
import xlrd

def process_report(report_filename, directory):        
    report_path = './' + os.fsdecode(directory) + '/' + report_filename
    sheet = xlrd.open_workbook(report_path).sheet_by_index(0)
    readCSV = []
    
    for row in range(sheet.nrows):
        readCSV.append(sheet.row_values(row))
    
    print(report_filename)
    
    row = 0
    file_header = []
    variant_data = []
    
    while readCSV[row][0][0] == '#':
        file_header.append(readCSV[row])
        row_variant_start = row
        if row < len(readCSV)-1:
            row += 1
        else:
            break
    
    while readCSV[row][0]:
        variant_data.append(readCSV[row])
        row_variant_end = row
        if row < len(readCSV)-1:
            row += 1
        else:
            break
    
    num_variants = row_variant_end - row_variant_start
    print('Number of variants = ' + str(num_variants))
    
    column_headers = file_header[-1]
    variant_ID = column_headers.index('ID')
    pathogenicity = column_headers.index('Pathogenicity')
                    
    patho_dict = {'Benign':'6', 'Likely Benign':'4', 'Unknown':'3', 'Possibly damaging':'8',
                  'Reported mutation':'9', 'Predicted benign':'12', 'Splice site mutation':'13',
                  'Frameshift':'14', 'False Positive':'15', 'Nonsense Mutation':'16', 'In-Frame':'17',
                  'Probably Damaging':'21', 'Unassigned':'25', 'Pathogenic':'26', 'Hypomorphic allele':'22',
                  'VOUS Clinvar':'27', 'VOUS':'27', 'VUS Clinvar':'27', 'VUS':'27', 'Conflicting Interpretations':'28',
                  'Co-segregated':'29', 'Likely pathogenic':'30', 'Start codon mutation':'31', 'Stop codon mutation':'32'}
    patho_dict_lower = {k.lower():v for k,v in patho_dict.items()}
    
    patho_pattern = re.compile(r'[a-zA-Z -]+')
    
    SQL_list = []
    unmatched = 0
    NA = 0
    unassigned = 0
    for variant in variant_data:
        if variant[1] != 'N/A' and variant[1]:
            patho_match = patho_pattern.match(variant[pathogenicity])
            formatted_patho = patho_match.group().strip().lower()
            if formatted_patho != 'unassigned':
                try:
                    SQL_list.append("UPDATE softgenetics.histdb_panelgroupvariant SET latest_pathogenicity_id = '"
                                + patho_dict_lower[formatted_patho] + "' WHERE variant_id = '" + str(int(variant[variant_ID])) + "';")
                except IndexError:
                    unmatched += 1
#                    print('Unmatched variant ' + variant[1] + ' ' + formatted_patho)
            else:
                unassigned += 1
        else:
#            print(variant[0:2])
            NA += 1
    
    print('Unassigned = ' + str(unassigned) + ', Unmatched = ' + str(unmatched) + ', N/A = ' + str(NA))
    print('Length of SQL = ' + str(len(SQL_list)))
    if num_variants - unassigned - unmatched - NA == len(SQL_list):
        print('QC Pass')
    else:
        print('QC Fail')
        
#    with open('./Generated SQL Files/' + report_filename[:-4] + '.sql', 'w') as fo:
#        for SQL in SQL_list:
#            fo.write(SQL + '\n')
    
    with open('./Generated SQL Files/batch.sql', 'a') as batch:
        for SQL in SQL_list:
            batch.write(SQL + '\n')
    print()
      
    return


def main():
    directory = os.fsencode('./Reports to Update')
    
    with open('./Generated SQL Files/batch.sql', 'w') as batch:
        pass
    
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('Report.xlsx'): 
            process_report(filename, directory)
    
    batch_length = 0
    with open('./Generated SQL Files/batch.sql', 'r') as batch:
        for count, line in enumerate(batch):
            batch_length += 1
        print('Batch SQL length = ' + str(batch_length))
        
    return


if __name__ == '__main__':
    main()
    
    