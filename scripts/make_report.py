#!/usr/bin/env python

import sys
from time import gmtime, strftime
import os
import re
from cyvcf2 import VCF
import tempfile
import csv
import binascii
import argparse

def parse_csq_header(vcf_file):
    for header in vcf_file.header_iter():
        info = header.info(extra=True)
        if b'ID' in info.keys() and info[b'ID'] == b'CSQ':
            format_pattern = re.compile('Format: (.*)"')
            match = format_pattern.search(info[b'Description'].decode())
            return match.group(1).split('|')

def parse_csq_entries(csq_entries, csq_fields):
    transcripts = {}
    for entry in csq_entries:
        values = entry.split('|')
        transcript = {}
        for key, value in zip(csq_fields, values):
            transcript[key] = value
        if transcript['Allele'] not in transcripts.keys():
            transcripts[transcript['Allele']] = []
        transcripts[transcript['Allele']].append(transcript)
    return transcripts

def get_csq_entries_bygene(csq_entries):
    genes = {}
    for entry in csq_entries:
        genes[entry['SYMBOL']] = entry
        
    return genes

def resolve_alleles(entry, csq_alleles):
    alleles = {}
    if entry.is_indel:
        for alt in entry.ALT:
            alt = str(alt)
            if alt[0:1] != entry.REF[0:1]:
                csq_allele = alt
            elif alt[1:] == "":
                csq_allele = '-'
            else:
                csq_allele = alt[1:]
            alleles[alt] = csq_allele
    elif entry.is_sv:
        for alt in alts:
            if len(alt) > len(entry.REF) and 'insertion' in csq_alleles:
                alleles[alt] = 'insertion'
            elif len(alt) < len(entry.REF) and 'deletion' in csq_alleles:
                alleles[alt] = 'deletion'
            elif len(csq_alleles) == 1:
                alleles[alt] = list(csq_alleles)[0]
    else:
        for alt in entry.ALT:
            alt = str(alt)
            alleles[alt] = alt
    return alleles

def get_genes(transcripts):
    genes = set()
    for transcript in transcripts:
        genes.add(transcript['SYMBOL'])

    return genes

def decode_hex(string):
    hex_string = string.group(0).replace('%', '')
    return binascii.unhexlify(hex_string).decode('utf-8')

def convert_aa(codon):
    three = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter"]
    one  = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
    
    for i in range(0,len(three)):
        p = re.compile(three[i])
        codon = p.sub(one[i],codon)

    return codon


#
# Script
#

MinGeneCov = 20
MinRegionCov = 20
MinFracCov = 90
MinReads = 3
MinVAF = 5.0

parser = argparse.ArgumentParser(description='Make ChromoSeq report.')
parser.add_argument('name',help='Sample name')
parser.add_argument('genevcffile',help='Gene VCF')
parser.add_argument('svvcffile',help='SV VCF')
parser.add_argument('genelist',help='Gene list')
parser.add_argument('mapsum',help='Map summary')
parser.add_argument('genecov',help='Gene coverage')
parser.add_argument('svcov',help='SV coverage')
parser.add_argument('haplotect',help='Haplotect output')
parser.add_argument('exception',help='Exception')
parser.add_argument('-v',"--minvaf",help='Minimum validated VAF')
parser.add_argument('-r',"--minreads",type=int,help='Minimum validated variant supporting reads')
parser.add_argument('-g',"--mingenecov",type=int,help='Min gene coverage')
parser.add_argument('-s',"--minsvcov",type=int,help='Min SV region coverage')
parser.add_argument('-f',"--fraccovqc",type=int,help='Min fraction at coverage')

args = parser.parse_args()

Name = args.name
genevcf_file = args.genevcffile
svvcf_file = args.svvcffile
genelist = args.genelist
mapsum = args.mapsum
genecov = args.genecov
svcov = args.svcov
haplotect = args.haplotect
exception = args.exception

if args.mingenecov:
    MinGeneCov = args.mingenecov

if args.minsvcov:
    MinRegionCov = args.minsvcov

if args.fraccovqc:
    MinFracCov = args.fraccovqc

if args.minreads:
    MinReads = args.minreads

if args.minvaf:
    MinVAF = args.minvaf
                
# variants to print out
vars = {}
vars['knownsv'] = []
vars['cna'] = []
vars['genelevel'] = []
vars['filteredgenelevel'] = []
vars['novelsv'] = []

mapdata = {}
rgmapdata = {}

mapdata['Average sequenced coverage over genome'] = ''
mapdata['Average alignment coverage over genome'] = ''
mapdata['Total input reads'] = ''
mapdata['Mapped reads'] = ''
mapdata['Number of duplicate marked reads'] = '';
mapdata['Number of unique reads (excl. duplicate marked reads)'] = ''
mapdata['Properly paired reads'] = ''
mapdata['Insert length: mean'] = ''
      
# read in mapping summary data
with open(mapsum, 'r') as m:
    reader = csv.reader(m, delimiter=',')
    for row in reader:
        mapdata[row[0]] = row[1]
    
# read genelist
knowngenelist = set()

with open(genelist, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        l = row[0]
        p = re.compile("_exon_\d+")
        l = p.sub("",l)
        knowngenelist.add(l)

# get gene-level variants

genevcf = VCF(genevcf_file)
genecsq_fields = parse_csq_header(genevcf)

for variant in genevcf:

    csq = variant.INFO.get('CSQ')
    
    if csq is None:
        sys.exit("No VEP fields")
                    
    transcripts = list(parse_csq_entries(csq.split(','), genecsq_fields).items())[0][1] # just get the first allele in the list.
    genes = get_csq_entries_bygene(transcripts)

    vartype = ''
    if len(variant.REF) == len(variant.ALT[0]):
        vartype = 'SNV'
    else:
        vartype = 'INDEL'

    filter = 'PASS'
    if variant.FILTER is not None:
        filter = variant.FILTER

    chr1 = str(variant.CHROM)
    pos1 = variant.POS
    svlen = len(variant.ALT) - len(variant.REF)

    gene = None
    if len(knowngenelist.intersection(get_genes(transcripts))) > 0:
        gene = knowngenelist.intersection(get_genes(transcripts)).pop()

    # otherwise dont report synonymous or filtered variants
    if gene is None or re.match("^synonymous|^5_prime_UTR|^3_prime_UTR|^upstream|^downstream|^intron",genes[gene]['Consequence']) is not None:
        continue

    abundance = variant.format("VAF")[0][0] * 100
        
    csyntax = 'NA'
    if genes[gene]['HGVSc'] is not None and genes[gene]['HGVSc'] is not '':
       csyntax = genes[gene]['HGVSc'].split(":")[1]
       
    psyntax = 'NA'
    if genes[gene]['HGVSp'] is not None and genes[gene]['HGVSp'] is not '':
        psyntax = genes[gene]['HGVSp'].split(":")[1]

    if psyntax is 'NA':
        psyntax = csyntax
       
    pmaf = genes[gene]['MAX_AF']
    if pmaf is None or pmaf == '':
        pmaf = 'NA'
    else:
        pmaf = str(float(pmaf) * 100) + '%'
        
    psyntax = convert_aa(psyntax)

    if abundance < float(MinVAF)*100 or int(variant.format("NV")[0][0]) < MinReads:
        filter = 'LowReads'

    if filter != 'PASS':
        vars['filteredgenelevel'].append([vartype,chr1,str(pos1),variant.REF,variant.ALT[0],gene,genes[gene]['Consequence'],csyntax,psyntax,str(genes[gene]['EXON']),filter,
                                          str(variant.ID),str(round(abundance,1))+"%",str(variant.format("NV")[0][0]),str(variant.format("NR")[0][0]),pmaf])
    else:
        vars['genelevel'].append([vartype,chr1,str(pos1),variant.REF,variant.ALT[0],gene,genes[gene]['Consequence'],csyntax,psyntax,str(genes[gene]['EXON']),filter,
                                  str(variant.ID),str(round(abundance,1))+"%",str(variant.format("NV")[0][0]),str(variant.format("NR")[0][0]),pmaf])
    
    
# done getting gene variants

# get SVs. have to go through twice because BNDs have 2 entries per

passedvars = {} # the ones that passed all filters
alreadydone = set() # so that BND mates can be skipped appropriately 

svvcf = VCF(svvcf_file)
svcsq_fields = parse_csq_header(svvcf)

for variant in svvcf:

    passedvars[variant.ID] = variant

# done with first pass

for v in passedvars.items():

    out = []
    
    variant = v[1]
    
    isknown = 0
    knowngenes = []
    mate = '';

    if variant.INFO.get('KNOWNSV') is not None:
        isknown = 1
        
    vartype = variant.INFO.get('SVTYPE')

    csq = variant.INFO.get('CSQ')

    if csq is None:
        sys.exit("No VEP fields")

    transcripts = list(parse_csq_entries(csq.split(','), svcsq_fields).items())[0][1] # just get the first allele in the list.
    genes = get_csq_entries_bygene(transcripts)
                                    
    filter = 'PASS'
    if variant.FILTER is not None:
        filter = variant.FILTER

    if vartype in ('DEL','DUP','INV'):
        
        chr1 = str(variant.CHROM)
        pos1 = variant.POS
        
        chr2 = chr1
        pos2 = int(variant.INFO.get('END'))
        svlen = pos2 - pos1 + 1
        
        gs = get_genes(transcripts)
        if '' in gs:
            gs.remove('')

        if len(gs) == 0:
            gs.add('None')
            
        bands = transcripts[0]['cytobands'].split("&")
        bandstr = bands[0]
        if len(bands) > 1:
            bandstr = bands[0] + bands[-1]

        knowngenes = list(gs.intersection(knowngenelist))
        knowngenes.sort()
        knowngenestring = ",".join(knowngenes)
        if len(knowngenes) == 0:
            knowngenestring = 'None'
            
        GS = list(gs)
        GS.sort()
        genestring = ",".join(GS)
        if len(genestring) > 10:
            genestring = str(len(gs)) + " genes"
            
        # For example:  seq[GRCh38] del(X)(q21.31q22.1)
        #          chrX:g.89555676_100352080del

        csyntax = '.'
        psyntax = '.'
        if vartype == 'DEL':
            csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "del"
            if bands[0].find('p') > -1 and bands[-1].find('q') > -1: # if the CNA spans the centromere then the whole chromosome is lost/gained
                psyntax = "seq[GRCh38] -" + chr1.replace('chr','')
                
            elif bands[0].find('p') > -1:
                bands.reverse()
                psyntax = "seq[GRCh38] del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            else:
                psyntax = "seq[GRCh38] del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
            
        elif vartype == 'DUP':
            csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "dup"
            if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
                psyntax = "seq[GRCh38] +" + chr1.replace('chr','')
                
            elif bands[0].find('p') > -1:
                bands.reverse()
                psyntax = "seq[GRCh38] dup(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            else:
                psyntax = "seq[GRCh38] dup(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"

        elif vartype == 'INV':
            csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "inv"
            if bands[0].find('p') > -1:
                bands.reverse()
                psyntax = "seq[GRCh38] inv(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            else:
                psyntax = "seq[GRCh38] inv(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
            
        # abundance
        abundance = 0.0
        pr = (0,0)
        sr = (0,0)
        if variant.INFO.get('LOG2RATIO') is not None:
            
            abundance = variant.INFO.get('ABUNDANCE') #((2**variant.INFO.get('LOG2RATIO') - 1.0) / ((CN/2.0 - 1.0)))*100;
            if abundance > 100:
                abundance = ">95%"
            else:   
                abdunance = str(round(abundance,1))+"%"
                
            infostring = 'CN=' + str(variant.INFO.get('CN')) + ';LOG2RATIO=' + str(round(variant.INFO.get('LOG2RATIO'),3))
            out = [vartype,chr1,str(pos1),chr2,str(pos2),str(svlen),bandstr,knowngenestring,csyntax,psyntax,genestring,filter,str(variant.ID),abundance,infostring]
                
        elif variant.format("SR") is not None and variant.format("PR")[0] is not None:
            sr = variant.format("SR")[0]
            pr =  variant.format("PR")[0]
                
            abundance = (sr[1] + pr[1]) / (pr[0] + pr[1] + sr[0] + sr[1]) * 100
            
            numhits = '.'
            if (variant.INFO.get('CONTIGHITS') is not None):
                numhits = variant.INFO.get('CONTIGHITS')

            infostring = 'PR_READS=' + str(pr[1]) + '/' + str(pr[0]+pr[1]) + ';SR_READS=' + str(sr[1]) + '/' + str(sr[0]+sr[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG'))
            out = [vartype,chr1,str(pos1),chr2,str(pos2),str(svlen),bandstr,knowngenestring,csyntax,psyntax,genestring,filter,str(variant.ID),str(round(abundance,1))+"%",infostring]
                
    elif vartype == 'BND':

        chr1 = str(variant.CHROM)
        pos1 = variant.POS
                        
        # skip if this is the mate
        if variant.INFO.get('MATEID') in alreadydone:
            continue
            
        # get the mate
        mate = passedvars[variant.INFO.get('MATEID')]

        if mate.INFO.get('KNOWNSV') is not None:
            isknown = 1
            
        matecsq = mate.INFO.get('CSQ')

        if matecsq is None:
            sys.exit("No VEP fields")
                
        matetranscripts = list(parse_csq_entries(matecsq.split(','), svcsq_fields).items())[0][1] # just get the first allele in the list.
        mategenes = get_csq_entries_bygene(matetranscripts)
            
        chr2 = mate.CHROM
        pos2 = mate.POS
        
        gs1 = get_genes(transcripts)
        if '' in gs1:
            gs1.remove('')
            
        if len(gs1) == 0:
            gs1.add('INTERGENIC')
        
        bands1 = transcripts[0]['cytobands'].split("&")
        
        gs2 = get_genes(matetranscripts)
        if '' in gs2:
            gs2.remove('')
            
        if len(gs2) == 0:
            gs2.add('INTERGENIC')
            
        bands2 = matetranscripts[0]['cytobands'].split("&")
        bandstr = bands1[0] + bands2[0]
        
        GS1 = list(gs1)
        GS2 = list(gs2)
        GS1.sort()
        GS2.sort()
        genestring = ",".join(GS1) + "--" + ",".join(GS2)
                                
        # abundance
        abundance = 0.0
        pr = (0,0)
        sr = (0,0)            
        if variant.format("SR") is not None:
            sr = variant.format("SR")[0]
            
        if variant.format("PR")[0] is not None:                
            pr =  variant.format("PR")[0]

        abundance = (sr[1] + pr[1]) / (pr[0] + pr[1] + sr[0] + sr[1]) * 100

        numhits = '.'
        if (variant.INFO.get('CONTIGHITS') is not None):
            numhits = variant.INFO.get('CONTIGHITS')

        knowngenes1 = gs1.intersection(knowngenelist)
        knowngenes2 = gs2.intersection(knowngenelist)
        knowngenes = list(knowngenes1.union(knowngenes2))
        knowngenes.sort() 

        knowngenestring = ",".join(knowngenes)
        if len(knowngenes) == 0:
            knowngenestring = 'None'
                                
        alt = variant.ALT[0]
        strand = '+'
        if alt.find("[") == 0 or alt.find("]") > 0:
            strand = '-'
            
        csyntax = '';
        psyntax = '';
        if (chr1.find('X') == -1 and chr2.find('X') == -1 and chr1.find('Y') == -1 and chr2.find('Y') == -1 and int(chr1.replace('chr','')) < int(chr2.replace('chr',''))) or chr1.find('X') > -1 or chr1.find('Y') > -1: # this isnt working. Want to list lower chromosome first in these strings. If X is involved, then X first.
            csyntax = chr1 + ":g." + str(pos1) + "(+)::" + chr2 + ":g." + str(pos2) + "(" + strand + ")"
            psyntax = 'seq[GRCh38] t(' + chr1.replace('chr','') + ';' + chr2.replace('chr','') + ')(' + bands1[0] + ';' + bands2[0] + ')'
        else:
            csyntax = chr2 + ":g." + str(pos2) + "(+)::" + chr1 + ":g." + str(pos1) + "(" + strand + ")"
            psyntax = 'seq[GRCh38] t(' + chr2.replace('chr','') + ';' + chr1.replace('chr','') + ')(' + bands2[0] + ';' + bands1[0] + ')'

        infostring = 'PR_READS=' + str(pr[1]) + '/' + str(pr[0]+pr[1]) + ';SR_READS=' + str(sr[1]) + '/' + str(sr[0]+sr[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG')) + ';'
        out = [vartype,chr1,str(pos1),chr2,str(pos2),"NA",bandstr,knowngenestring,csyntax,psyntax,genestring,filter,str(variant.ID) + ";" + str(mate.ID),str(round(abundance,1))+"%",infostring]

        alreadydone.add(variant.ID)
        
    if isknown == 1:
        vars['knownsv'].append(out)
        
    elif variant.INFO.get('LOG2RATIO') is not None and filter is 'PASS':
        vars['cna'].append(out)

    else:
        vars['novelsv'].append(out)


print("ChromoSeq Report for " + Name + " ---- Generated on: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\n")

print("*** COPY NUMBER ALTERATIONS ***\n")
if len(vars['cna']) > 0:
    print("\t".join(('TYPE','CHR1','POS1','CHR2','POS2','LENGTH','BANDS','KNOWN_GENES','HGVS-LIKE','ISCN-LIKE','TOTAL_GENES','FILTERS','ID','ABUNDANCE','INFO')))
    for i in vars['cna']:
        print("\t".join(i))
    print()        
else:
    print("\tNONE DETECTED\n")    

print("*** RECURRENT TRANSLOCATIONS ***\n")
if len(vars['knownsv']) > 0:
    print("\t".join(('TYPE','CHR1','POS1','CHR2','POS2','LENGTH','BANDS','KNOWN_GENES','HGVS-LIKE','ISCN-LIKE','TOTAL_GENES','FILTERS','ID','ABUNDANCE','INFO')))
    for i in vars['knownsv']:
        print("\t".join(i))
    print()
else:
    print("\tNONE DETECTED\n")    

print("*** GENE MUTATIONS ***\n")
if len(vars['genelevel']) > 0:
    print("\t".join(('TYPE','CHR','POS','REF','ALT','GENE','CONSEQUENCE','HGVSc','HGVSp','EXON','FILTERS','ID','VAF','VariantReads','TotalReads','PopulationMAF')))
    for r in vars['genelevel']:
        print("\t".join(r))
    print()
else:
    print("\t"+"NONE DETECTED\n")

print("*** FILTERED GENE MUTATIONS ***\n")
if len(vars['filteredgenelevel']) > 0:
    print("\t".join(('TYPE','CHR','POS','REF','ALT','GENE','CONSEQUENCE','HGVSc','HGVSp','EXON','FILTERS','ID','VAF','VariantReads','TotalReads','PopulationMAF')))
    for r in vars['filteredgenelevel']:
        print("\t".join(r))
    print()
else:
    print("\t"+"NONE DETECTED\n")

print("*** OTHER STRUCTURAL VARIANTS ***\n")
if len(vars['novelsv']) > 0:
    print("\t".join(('TYPE','CHR1','POS1','CHR2','POS2','LENGTH','BANDS','KNOWN_GENES','HGVS-LIKE','ISCN-LIKE','TOTAL_GENES','FILTERS','ID','ABUNDANCE','INFO')))
    for r in vars['novelsv']:
        print("\t".join(r))
    print()
else:
    print("\t"+"NONE DETECTED\n")

## now print QC stuff

print("*** CHROMOSEQ QC ***\n")

print("AVERAGE COVERAGE:\t",mapdata['Average alignment coverage over genome'])
print("TOTAL READS:\t",mapdata['Total input reads'])
print("MAPPED READS:\t",mapdata['Mapped reads']," ",str(round(int(mapdata['Mapped reads'])/int(mapdata['Total input reads'])*100,1))+"%")
print("DUPLICATES:\t",mapdata['Number of duplicate marked reads']," ",str(round(int(mapdata['Number of duplicate marked reads'])/int(mapdata['Total input reads'])*100,1))+"%")
print("UNIQUE READS:\t",mapdata['Number of unique reads (excl. duplicate marked reads)']," ",str(round(int(mapdata['Number of unique reads (excl. duplicate marked reads)'])/int(mapdata['Total input reads'])*100,1))+"%")
print("PROPERLY PAIRED READS:\t"+str(mapdata['Properly paired reads'])+" ",str(round(int(mapdata['Properly paired reads'])/int(mapdata['Total input reads'])*100,1))+"%")
print("MEAN INSERT SIZE:\t"+mapdata['Insert length: mean']+"\n")

print("*** Haplotect Contamination Estimate ***\n")
with open(haplotect, 'r') as hap:
    reader = csv.reader(hap, delimiter='\t')
    for row in reader:
        print("\t".join(row))

print("\n*** Exon Coverage Metrics: Exons with <"+str(MinFracCov)+"% at " + str(MinGeneCov) + "x ***\n")
# gene cov
gc = {}
gcheader = []
covlevelindex = -1
with open(genecov, 'r') as g:
    reader = csv.reader(g, delimiter='\t')
    line_count = 0
    for row in reader:
        if line_count == 0:
            print("\t".join(row))
            gcheader = row[4:]
            covlevelindex = row.index("%"+str(MinGeneCov)+"X")
            line_count += 1
            continue
        
        elif float(row[covlevelindex]) < MinFracCov: # this is for 
            print("\t".join(row))
        
        # get gene name
        gn = row[3].split("_")[0]
        if gn not in gc.keys():
            gc[gn] = {'l':0, 'c':0.0, '10x':0, '20x':0, '30x':0, '40x':0}
            
        # sum coverages for each exon
        gc[gn]['l'] += int(row[2])-int(row[1])
        gc[gn]['c'] += float(row[12]) * (float(row[2])-float(row[1]))
        gc[gn]['10x'] += int(row[4])
        gc[gn]['20x'] += int(row[5])
        gc[gn]['30x'] += int(row[6])
        gc[gn]['40x'] += int(row[7])

        line_count += 1
        
print("\n*** Gene Coverage Metrics: Genes with <"+str(MinFracCov)+"% at " + str(MinGeneCov) + "x ***\n")
# print gene coverage metrics if lower than MinFracGene or MinGeneCov
print("#gene\tlength\t" + "\t".join(gcheader))
for g in gc.keys():
    if gc[g][str(MinGeneCov)+"x"] / gc[g]['l'] * 100 < MinFracCov:
        print('{}\t{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}'.format(g,gc[g]['l'],gc[g]['10x'],gc[g]['20x'],gc[g]['30x'],gc[g]['40x'],gc[g]['10x']/gc[g]['l']*100,gc[g]['20x']/gc[g]['l']*100,gc[g]['30x']/gc[g]['l']*100,gc[g]['40x']/gc[g]['l']*100,gc[g]['c']/gc[g]['l']))

print("\n*** SV Region Coverage Metrics: Exons with <"+str(MinFracCov)+"% at " + str(MinRegionCov) + "x ***\n")
with open(svcov, 'r') as sv:
    reader = csv.reader(sv, delimiter='\t')
    line_count = 0
    for row in reader:
        if line_count == 0:
            print("\t".join(row))
        elif float(row[8]) < MinFracCov or float(row[12]) < MinRegionCov:
            print("\t".join(row))
            
        line_count += 1

print("\n*** EXCEPTIONS ***\n")
print(exception + "\n")

print("\nThis laboratory developed test (LDT) was developed and its performance characteristics determined by the CLIA Licensed Environment laboratory at the McDonnell Genome Institute at Washington University (MGI-CLE, CLIA #26D2092546, CAP #9047655), Dr. David H. Spencer MD, PhD, FCAP, Medical Director. 4444 Forest Park Avenue, Rm 4111 St. Louis, Missouri 63108 (314) 286-1460 Fax: (314) 286-1810. The MGI-CLE laboratory is regulated under CLIA as certified to perform high-complexity testing. This test has not been cleared or approved by the FDA.")
