# -*- coding: utf-8 -*-
"""
Created on Thu Nov 9 12:49:23 2017

@author: nkwang

Batch primer design from list of target region coordinates
"""

import csv
import requests
import os
from bs4 import BeautifulSoup
import bs4
import time
import urllib.parse

headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_5) AppleWebKit 537.36 (KHTML, like Gecko) Chrome',
           'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8'}
ptHeaders = { 'Origin': 'https://www.ncbi.nlm.nih.gov',
              'Referer': 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/',
              'Cache-Control': 'max-age=0'}
jobHeaders = { 'Origin': 'https://www.ncbi.nlm.nih.gov',
              'Referer': 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'}

referer = {'Referer': 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1510766788&job_key=UVuOAlNUXvx5xs7Dw6Pq8bm4-8OUq-DelQ'}
cookies = {}

PBBase = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/'
PTcgi = 'primertool.cgi'

rtf_dir = './rtfs'
logPath = './logs'
rtf_flank_bp = '2000' # how many bp flanking gene to include in rtf download, default=400
rtf_overwrite = 1 # overwrite rtf if previous version exists, default 0
region_pad_bp = 600 # how many bp flanking region to input into primer3, default=300
max_product_length = '450' # maximum PCR product length, default=500, up to 1000 for sanger, 20kb for LR-PCR
delay = 20  # how often to check job completion in seconds, default=20
deletion_analysis = 0 # 1 for large deletion analysis


def download_file(session, url, dirPath, geneName):
    r = session.get(url, cookies=cookies, headers=headers, stream=True)
    filename = dirPath + '/' + geneName + '.rtf'
    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=4096):
            if chunk:
                f.write(chunk)
    return filename


def download_rtf(geneName, dirPath):
    if not os.path.exists(dirPath):
        os.makedirs(dirPath)

    session = requests.Session()

    searchUrl = 'http://grch37.ensembl.org/Multi/Ajax/search?q=(+' + geneName + \
                '%5E190+AND+species%3A"Human"+)'
    r = session.get(searchUrl, cookies=cookies, headers=headers)
    geneInfo = r.json()['result']['response']['docs'][0]
    geneRange = geneInfo['location']
    geneID = geneInfo['id']
    print(geneName + '\t' + geneID + '\t' + geneRange)
    formData = ('compression=uncompressed' '&'
                'data_type=Gene' '&'
                'data_action=' '&'
                'component=GeneSeq' '&'
                'export_action=GeneSeq' '&'
                'adorn=both' '&'
                'r=' + geneRange + '&'
                'g=' + geneID + '&'
                'db=core' '&'
                'name=Homo_sapiens_' + geneName + '_sequence' '&'
                'format=RTF' '&'
                'extra_FASTA=sequence' '&'
                'flank5_display_FASTA=600' '&'
                'flank3_display_FASTA=600' '&'
                'flank5_display_RTF=' + rtf_flank_bp + '&'
                'flank3_display_RTF=' + rtf_flank_bp + '&'
                'display_width_RTF=60' '&'
                'exon_display_RTF=core' '&'
                'exon_ori_RTF=all' '&'
                'snp_display_RTF=on' '&'
                'hide_long_snps_RTF=on' '&'
                'hide_rare_snps_RTF=0.01' '&' # default = 0.001
                'consequence_filter_RTF=off' '&'
                'line_numbering_RTF=slice' '&'
                'variants_as_n_RTF=on')

    url = 'http://grch37.ensembl.org/Homo_sapiens/DataExport/Output?' + formData
    print('Downloaded file: ' + download_file(session, url, dirPath, geneName))


def format_rtf(rtf_filename):
    rtf_file = open(rtf_filename, 'r')

    formatted_list = []
    line_list = rtf_file.readlines()
    chr_start = []
    chr_end = []

    for index in list(range(1,len(line_list[11]))):
        if line_list[11][index].isnumeric():
            chr_start.append(line_list[11][index])
        else:
            break

    flag = 0
    for index in list(range(-1,len(line_list[-1]))):
        if line_list[-1][index].isnumeric():
            chr_end.append(line_list[-1][index])
            flag = 1
        else:
            if flag == 1:
                break

    chr_start = int(''.join(chr_start))
    chr_end = int(''.join(chr_end))
#    chr_len = abs(chr_end - chr_start) + 1

    for line in line_list[9:]:
        formatted_list.append(line.strip('\n'))

    rtf_string = ''.join(formatted_list).replace(' ', '')

    sequence = []
    for char in rtf_string:
        if char in "ATCGNul}":
            sequence.append(char)

    for index in list(range(0,len(sequence))):
        if sequence[index] == 'u' and sequence[index+1] == 'l':
            bracket = 0
            index_enum = 2
            while bracket == 0:
                if sequence[index+index_enum] != '}':
                    sequence[index+index_enum] = 'N'
                    index_enum += 1
                else:
                    bracket = 1

    sequence2 = []
    for char in sequence:
        if char in "ATCGN":
            sequence2.append(char)

    sequence_string = ''.join(sequence2)

#    print(chr_start, chr_end)
#    print('Sequence length = ' + str(len(sequence_string)))
#    if chr_len == len(sequence_string):
#        print('Length OK')
#    else:
#        print('LENGTH MISMATCH')

    return([sequence_string, chr_start, chr_end])


def outputTag(tag):
    output = ''
    if tag is not None:
        if 'class' in tag.parent.attrs.keys():
            #print(tag.parent['class'])
            if len(tag.parent['class']) > 1:
                if tag.parent['class'][1] != 'shown':
                    return output
        for item in tag.contents:
            if isinstance(item, bs4.element.NavigableString):
                output += str(item)
            elif isinstance(item, bs4.element.Tag):
                if item.name == 'a':
                    output += str(item.string) + ' ' + item['href']
                else:
                    output += outputTag(item)
    return output


def format_prPairInfo(tags):
    fStr = ''
    fp = ''
    rp = ''
    tagIndex = 0
    for tag in tags:
        h2Tag = tag.h2
        fStr += outputTag(h2Tag) + '\n'
        trs = tag.table.find_all('tr',{'class':None})
        rowIndex = 0
        colWidth = [20, 25, 16, 7, 6, 6, 6, 6, 21, 24]
        for tr in trs:
            row = ''
            strings = list(tr.strings)
            if rowIndex == 0:
                fStr += '-' * sum(colWidth) + '\n'
                row += ' ' * colWidth[0]
                for index in range(1, len(strings) + 1):
                    row +=  strings[index-1] + ' ' * (colWidth[index] - len(strings[index-1]))
            else:
                for index in range(0, len(strings)):
                    row +=  strings[index] + ' ' * (colWidth[index] - len(strings[index]))
            fStr += row + '\n'

            # extract first primer pair
            if tagIndex == 0:
                if rowIndex == 1:
                    fp = strings[1]
                if rowIndex == 2:
                    rp = strings[1]
                if rowIndex == 3:
                    product_len = strings[1]

            rowIndex += 1
        fStr += '\n'
        prPairTl = tag.find_all('div', {'class':'prPairTl'})
        for ppt in prPairTl:
            fStr += outputTag(ppt)
        prPairDtl = tag.find_all('div', {'class':'prPairDtl'})
        for ppd in prPairDtl:
            fStr += outputTag(ppd)
        tagIndex += 1
    return [fStr, fp, rp, product_len]


def generate_primer(sequence, region_coord_left, region_coord_right,
                    ext_coord_left, ext_coord_right):
    global deletion_analysis
    
    if sequence[2] > sequence[1]:  # forward sequence
        input_length = ext_coord_right - ext_coord_left + 1
        start_index = ext_coord_left - sequence[1]
        end_index = input_length + start_index
        
        if deletion_analysis == 0:
            input_sequence = sequence[0][start_index:end_index] # default
        else:
            input_sequence = sequence[0][start_index:region_coord_left-sequence[1]] + sequence[0][region_coord_right-sequence[1]:end_index] # for large deletion analysis
        region_left_index = region_coord_left - ext_coord_left + 1  # - input_sequence.count('N',0,region_coord_left)
        region_right_index = region_coord_right - ext_coord_left + 1  # - input_sequence.count('N',0,region_coord_right)

        #        excluded_regions = []
        #        for index, char in enumerate(input_sequence):
        #            if char == 'N':
        #                excluded_regions.append(str(index+1) + '%2C1%20')
        #        excluded_regions = ''.join(excluded_regions)

    else:  # reverse sequence
        input_length = int(ext_coord_right) - int(ext_coord_left) + 1
        start_index = sequence[1] - ext_coord_right + 1
        end_index = input_length + start_index
        input_sequence = sequence[0][start_index:end_index]
        region_left_index = ext_coord_right - region_coord_right + 1
        region_right_index = ext_coord_right - region_coord_left + 1

    inputURL = ('https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=bookmark' '&'
                    'INPUT_SEQUENCE=' + input_sequence + '&'
                    'PRIMER5_START=' + '1' + '&'
                    'PRIMER5_END=' + str(region_left_index) + '&'
                    'PRIMER3_START=' + str(region_right_index) + '&'
                    'PRIMER3_END=' + str(input_length) + '&'
#                    'PRIMER3_END=' + str(input_length - input_sequence.count('N',0,input_length-1)) + '&'
#                    'EXCLUDED_REGIONS=' + excluded_regions + '&'
                    'OVERLAP_5END=7' '&'
                    'OVERLAP_3END=4' '&'
                    'PRIMER_PRODUCT_MIN=50' '&'
                    'PRIMER_PRODUCT_MAX=' + max_product_length + '&' ###
                    'PRIMER_NUM_RETURN=5' '&'###
                    'PRIMER_MIN_TM=57.0' '&'
                    'PRIMER_OPT_TM=60.0' '&'
                    'PRIMER_MAX_TM=63.0' '&'
                    'PRIMER_MAX_DIFF_TM=3' '&'
                    'PRIMER_ON_SPLICE_SITE=0' '&'
                    'SEARCHMODE=0' '&'
                    'SPLICE_SITE_OVERLAP_5END=7' '&'
                    'SPLICE_SITE_OVERLAP_3END=4' '&'
                    'SPAN_INTRON=off' '&'
                    'MIN_INTRON_SIZE=1000' '&'
                    'MAX_INTRON_SIZE=1000000' '&'
                    'SEARCH_SPECIFIC_PRIMER=off' '&' #default = on
                    'EXCLUDE_ENV=off' '&'
                    'EXCLUDE_XM=off' '&'
                    'TH_OLOGO_ALIGNMENT=off' '&'
                    'TH_TEMPLATE_ALIGNMENT=off' '&'
                    'ORGANISM=Homo%20sapiens' '&'
                    'PRIMER_SPECIFICITY_DATABASE=refseq_representative_genomes' '&'
                    'TOTAL_PRIMER_SPECIFICITY_MISMATCH=1' '&'
                    'PRIMER_3END_SPECIFICITY_MISMATCH=0' '&'
                    'MISMATCH_REGION_LENGTH=10' '&'
                    'TOTAL_MISMATCH_IGNORE=6' '&'
                    'MAX_TARGET_SIZE=4000' '&'
                    'ALLOW_TRANSCRIPT_VARIANTS=off' '&'
                    'HITSIZE=50000' '&'
                    'EVALUE=30000' '&'
                    'WORD_SIZE=7' '&'
                    'MAX_CANDIDATE_PRIMER=500' '&'
                    'PRIMER_MIN_SIZE=15' '&'
                    'PRIMER_OPT_SIZE=20' '&'
                    'PRIMER_MAX_SIZE=25' '&'
                    'PRIMER_MIN_GC=20.0' '&'
                    'PRIMER_MAX_GC=80.0' '&'
                    'GC_CLAMP=0' '&'
                    'NUM_TARGETS_WITH_PRIMERS=1000' '&'
                    'NUM_TARGETS=20' '&'
                    'MAX_TARGET_PER_TEMPLATE=100' '&'
                    'POLYX=5' '&'
                    'SELF_ANY=8.00' '&'
                    'SELF_END=3.00' '&'
                    'PRIMER_MAX_END_STABILITY=9' '&'
                    'PRIMER_MAX_END_GC=5' '&'
                    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00' '&'
                    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00' '&'
                    'PRIMER_MAX_SELF_ANY_TH=45.0' '&'
                    'PRIMER_MAX_SELF_END_TH=35.0' '&'
                    'PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0' '&'
                    'PRIMER_PAIR_MAX_COMPL_END_TH=35.0' '&'
                    'PRIMER_MAX_HAIRPIN_TH=24.0' '&'
                    'PRIMER_MAX_TEMPLATE_MISPRIMING=12.00' '&'
                    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00' '&'
                    'PRIMER_PAIR_MAX_COMPL_ANY=8.00' '&'
                    'PRIMER_PAIR_MAX_COMPL_END=3.00' '&'
                    'PRIMER_MISPRIMING_LIBRARY=AUTO' '&'
                    'NO_SNP=off' '&'
                    'LOW_COMPLEXITY_FILTER=on' '&'
                    'MONO_CATIONS=50.0' '&'
                    'DIVA_CATIONS=1.5' '&'
                    'CON_ANEAL_OLIGO=50.0' '&'
                    'CON_DNTPS=0.6' '&'
                    'SALT_FORMULAR=1' '&'
                    'TM_METHOD=1' '&'
                    'PRIMER_INTERNAL_OLIGO_MIN_SIZE=18' '&'
                    'PRIMER_INTERNAL_OLIGO_OPT_SIZE=20' '&'
                    'PRIMER_INTERNAL_OLIGO_MAX_SIZE=27' '&'
                    'PRIMER_INTERNAL_OLIGO_MIN_TM=57.0' '&'
                    'PRIMER_INTERNAL_OLIGO_OPT_TM=60.0' '&'
                    'PRIMER_INTERNAL_OLIGO_MAX_TM=63.0' '&'
                    'PRIMER_INTERNAL_OLIGO_MAX_GC=80.0' '&'
                    'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50' '&'
                    'PRIMER_INTERNAL_OLIGO_MIN_GC=20.0' '&'
                    'PICK_HYB_PROBE=off' '&'
                    'NEWWIN=off' '&'
                    'NEWWIN=off' '&'
                    'SHOW_SVIEWER=true')

    session = requests.Session()
    
    if deletion_analysis == 0:
        primer3_start = str(region_right_index)
        primer3_end = str(input_length)
    else:
        primer3_start = '20001'
        primer3_end = '40000'
        
    mpData = \
        {
            'INPUT_SEQUENCE':(None,input_sequence),
            'SEQFILE': ('', '', 'application/octet-stream'),
            'PRIMER5_START': (None, '1'),
            'PRIMER5_END': (None, str(region_left_index)),
            'PRIMER3_START': (None, primer3_start),
            'PRIMER3_END': (None, primer3_end),
            'PRIMER_LEFT_INPUT': (None, ''),
            'PRIMER_RIGHT_INPUT': (None, ''),
            'PRIMER_PRODUCT_MIN': (None, '50'),
            'PRIMER_PRODUCT_MAX': (None, max_product_length),
            'PRIMER_NUM_RETURN': (None, '5'),
            'PRIMER_MIN_TM': (None, '55.0'), # default = 57
            'PRIMER_OPT_TM': (None, '60.0'), # default = 60
            'PRIMER_MAX_TM': (None, '65.0'), # default = 63
            'PRIMER_MAX_DIFF_TM': (None, '5'), # default = 3
            'PRIMER_ON_SPLICE_SITE': (None, '0'),
            'SPLICE_SITE_OVERLAP_5END': (None, '7'),
            'SPLICE_SITE_OVERLAP_3END': (None, '4'),
            'MIN_INTRON_SIZE': (None, '1000'),
            'MAX_INTRON_SIZE': (None, '1000000'),
            'SEARCH_SPECIFIC_PRIMER': (None, 'off'), #default = on
            'SEARCHMODE': (None, '0'),
            'PRIMER_SPECIFICITY_DATABASE': (None, 'refseq_mrna'), #default = refseq_representative_genomes
            'CUSTOM_DB': (None, ''),
            'CUSTOMSEQFILE': ('', '', 'application/octet-stream'),
            'ORGANISM': (None, 'Homo sapiens'),
            'ORG_DBS': (None, ''),
            'slctOrg': (None, ''),
            'ENTREZ_QUERY': (None, ''),
            'TOTAL_PRIMER_SPECIFICITY_MISMATCH': (None, '1'),
            'PRIMER_3END_SPECIFICITY_MISMATCH': (None, '0'),
            'MISMATCH_REGION_LENGTH': (None, '10'),
            'TOTAL_MISMATCH_IGNORE': (None, '6'),
            'MAX_TARGET_SIZE': (None, '4000'),
            'SHOW_SVIEWER': (None, 'on'),
            'HITSIZE': (None, '50000'),
            'UNGAPPED_BLAST': (None, 'on'),
            'EVALUE': (None, '30000'),
            'WORD_SIZE': (None, '7'),
            'MAX_CANDIDATE_PRIMER': (None, '500'),
            'NUM_TARGETS': (None, '20'),
            'NUM_TARGETS_WITH_PRIMERS': (None, '1000'),
            'MAX_TARGET_PER_TEMPLATE': (None, '100'),
            'PRODUCT_MIN_TM': (None, ''),
            'PRODUCT_OPT_TM': (None, ''),
            'PRODUCT_MAX_TM': (None, ''),
            'PRIMER_MIN_SIZE': (None, '15'),
            'PRIMER_OPT_SIZE': (None, '20'),
            'PRIMER_MAX_SIZE': (None, '25'),
            'PRIMER_MIN_GC': (None, '20.0'),
            'PRIMER_MAX_GC': (None, '80.0'),
            'GC_CLAMP': (None, '0'),
            'POLYX': (None, '5'),
            'PRIMER_MAX_END_STABILITY': (None, '9'),
            'PRIMER_MAX_END_GC': (None, '5'),
            'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': (None, '40.00'),
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': (None, '70.00'),
            'PRIMER_MAX_SELF_ANY_TH': (None, '45.0'),
            'PRIMER_MAX_SELF_END_TH': (None, '35.0'),
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': (None, '45.0'),
            'PRIMER_PAIR_MAX_COMPL_END_TH': (None, '35.0'),
            'PRIMER_MAX_HAIRPIN_TH': (None, '24.0'),
            'PRIMER_MAX_TEMPLATE_MISPRIMING': (None, '12.00'),
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': (None, '24.00'),
            'SELF_ANY': (None, '8.00'),
            'SELF_END': (None, '3.00'),
            'PRIMER_PAIR_MAX_COMPL_ANY': (None, '8.00'),
            'PRIMER_PAIR_MAX_COMPL_END': (None, '3.00'),
            'EXCLUDED_REGIONS': (None, ''),
            'OVERLAP': (None, ''),
            'OVERLAP_5END': (None, '7'),
            'OVERLAP_3END': (None, '4'),
            'MONO_CATIONS': (None, '50.0'),
            'DIVA_CATIONS': (None, '1.5'),
            'CON_DNTPS': (None, '0.6'),
            'SALT_FORMULAR': (None, '1'),
            'TM_METHOD': (None, '1'),
            'CON_ANEAL_OLIGO': (None, '50.0'),
#            'PRIMER_MISPRIMING_LIBRARY': (None, 'AUTO'),
            'LOW_COMPLEXITY_FILTER': (None, 'on'),
            'PRIMER_INTERNAL_OLIGO_MIN_SIZE': (None, '18'),
            'PRIMER_INTERNAL_OLIGO_OPT_SIZE': (None, '20'),
            'PRIMER_INTERNAL_OLIGO_MAX_SIZE': (None, '27'),
            'PRIMER_INTERNAL_OLIGO_MIN_TM': (None, '57.0'),
            'PRIMER_INTERNAL_OLIGO_OPT_TM': (None, '60.0'),
            'PRIMER_INTERNAL_OLIGO_MAX_TM': (None, '63.0'),
            'PRIMER_INTERNAL_OLIGO_MIN_GC': (None, '20.0'),
            'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT': (None, '50'),
            'PRIMER_INTERNAL_OLIGO_MAX_GC': (None, '80.0'),
            #'SHOW_SVIEWER': (None, 'on'),
            'LINK_LOC': (None, ''),
            'SVIEWER_DATA_KEY': (None, ''),
            'CMD': (None, 'request'),
            'NUM_DIFFS': (None, '5'),
            'NUM_OPTS_DIFFS': (None, '0')
        }
    primerPairList = []
    resultsHtml = ''

    session.headers.update(headers)
    r = session.get(PBBase)
    r = session.post(PBBase+PTcgi, headers=ptHeaders, files=mpData)
    soup = BeautifulSoup(r.text,'html.parser')

    # if the landing page is titled 'Primer-Blast results'
    # get submitted job progress tracking URL and keep ckecking job status until it's done
    firstJobURL = ''
    while soup.title.string == 'Primer-Blast results':
        formTag = soup.form
        if formTag is None:
            break
        else:
            firstJobURL = formTag['action']
#            print('First Job URL:')
#            print(firstJobURL)
            time.sleep(delay)   # wait before checking job status
            r = session.get(firstJobURL, headers=jobHeaders)
            soup = BeautifulSoup(r.text,'html.parser')
            #print("Response HTML:")
            #print(soup.prettify())

    # if we are getting sequence confirmation page, collect all params from the input fields
    # and submit the job again.
    if soup.title.string == 'Primer-Blast hit review':
        params = soup.find_all('input')
        postData = {}
        for param in params:
            if all(x in param.attrs for x in ['type', 'name', 'value']):
                postData[param['name']] = param['value']
        #print(urllib.parse.urlencode(postData))
        r = session.post(PBBase+PTcgi,headers=jobHeaders,data=urllib.parse.urlencode(postData))
        soup = BeautifulSoup(r.text,'html.parser')
        #print(soup.prettify())
        #jobURL = soup.form['action']
        #print(jobURL)

    # if the landing page is titled 'Primer-Blast results'
    # get submitted job progress tracking URL and keep ckecking job status until it's done
    secondJobURL = ''
    while soup.title.string == 'Primer-Blast results':
        formTag = soup.form
        if formTag is None:
            break
        else:
            secondJobURL = formTag['action']
#            if not secondJobURL:
#                print('Second Job URL:')
#                print(secondJobURL)
            time.sleep(delay)   # wait before checking job status
            r = session.get(secondJobURL, headers=jobHeaders)
            soup = BeautifulSoup(r.text,'html.parser')
            #print("Response HTML:")
            #print(soup.prettify())

    # we're getting the real result page
    if soup.title.string == 'Primer-Blast results':
        resultsHtml = r.text
        
        tags = soup.find_all('div', {'class': 'prPairInfo'})
        if len(tags) > 0:         # Primers Found
            print('Primers found')
            primerPairList = format_prPairInfo(tags)
        else:
            tags = soup.find_all('p', {'class': 'info'})
            if len(tags) > 0: # No primers found
                print('Primers not found')
                info = ''
                for tag in tags:
                    info += tag.string
                    primerPairList = [info, 'Not Found', 'Not Found']
            else:
                print('Unexpected page content:')
                print(r.text)
    else: # error info
        print('Unexpected page:')
        print(r.text)
		
    if secondJobURL:
        resultsURL = secondJobURL
    else:
        resultsURL = firstJobURL
    print(resultsURL)    
    
    return([resultsHtml, primerPairList, resultsURL, inputURL])


def main():
    with open('regions.csv', 'r') as f:
        regions = list(csv.reader(f))[1:]

    gene_list = []

    for region in regions:
        gene_list.append(region[0])
        region[1] = int(region[1])
        region[2] = int(region[2])
        ext_coord_left = region[1] - region_pad_bp
        ext_coord_right = region[2] + region_pad_bp
        region.append(ext_coord_left)
        region.append(ext_coord_right)

    sequence_info = {}
    gene_set = set(gene_list)
    
    for gene in gene_set:
        geneRtf = rtf_dir + '/' + gene + '.rtf'
        # Gene RTF file does not exist
        if not os.path.isfile(geneRtf):
            download_rtf(gene, rtf_dir)
        elif rtf_overwrite:
            print(geneRtf + ' exists. Overwrite file.')
            download_rtf(gene, rtf_dir)
        else:	# Gene RTF file exists. Skip download
            print(geneRtf + ' exists. Skip download.')
        sequence_info.update({gene: format_rtf(geneRtf)})

    # Make log directory
    if not os.path.exists(logPath):
        os.makedirs(logPath)

    with open('output.csv', 'w', newline='') as output:
        writer = csv.writer(output)
        writer.writerow(['Gene','Chr Start','Chr End', 'Forward Primer', 'Reverse Primer', 'Product Length', 'Results URL', 'Input URL'])
        for region in regions:
            print('\n' + region[0])
            (html, *results) = generate_primer(sequence_info[region[0]], region[1], region[2], region[3], region[4])
            primers = results[0]
            if len(primers) == 4:
                writer.writerow(region[0:3] + primers[1:] + list(results[1:]))
                print(region[0:3] + primers[1:])
                logName = logPath + '/' + region[0] + '-' + str(region[1]) + '-' + str(region[2])
                with open(logName + '.log', 'w') as log:
                    log.write(primers[0])
                with open(logName + '.html', 'w') as log:
                    log.write(html)
    return



if __name__ == '__main__':
    main()

