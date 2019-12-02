#!/usr/bin/env python

import csv, sys,re,os

def mean(list):
    sum = 0.0
    ct = 0.0

    for n in list:
        sum += n
        ct += 1
    return "%.4f"%(sum / ct)


pfam2gofile="pfam2go"
samples="samples.csv"
pfamtbl="M_tuberculosis.domtbl_simple.tab"

resultsdir="kallisto_results"
# read in pfam2go
pfam2go = {}

with open(pfam2gofile,"rU") as pf2go:
    for line in pf2go:
        line = line.strip()
        if line.startswith("!"):
            continue
        (pfam,go) = line.split(" > ")
        (pfamacc,pfamname)  = pfam.split(" ")
        pfamacc = re.sub("Pfam:","",pfamacc)
#        pfamacc2name[pfamacc] = pfamname
        (goname,goid) = go.split(" ; ")
        if pfamacc in pfam2go:
            pfam2go[pfamacc].append(goid)
        else:
            pfam2go[pfamacc] = [goid]

# read in expression
gene_expression = {}
gene_expression_avg = {} # bonus points to compute the mean for this experiment
experiments = set()
experiments_avg =  set()

with open(samples,"rU") as samp:
    sampcsv = csv.reader(samp,delimiter=",")
    header = next(sampcsv)
    for line in sampcsv:
        # make the experiment a full string matches folder name
        expname="%s_%s_r%s"%(line[1],line[2],line[3])
        # make experiment and ph a variable
        condition_ph = "%s_%s"%(line[1],line[2])
        experiments.add(expname)
        experiments_avg.add(condition_ph)

        filepath = os.path.join(resultsdir,expname,"abundance.tsv")
        if os.path.exists(filepath):
            with open(filepath,"rU") as kallistodata:
                kallisto = csv.reader(kallistodata,delimiter="\t")
                # treat first row differently since it is a header
                header = next(kallisto)
                tpm_col = -1
                i = 0
                # how to programmatically figure which column is the 'tpm' one
                for n in header:
                    if n == "tpm":
                        tpm_col = i
                    i += 1
                if tpm_col < 0:
                    print("expected to find column name 'tpm' in %s",filepath)
                    exit()
                # process each other line in the report
                for row in kallisto:
                    gene_name = row[0] # gene name is the first column
                    tpm = float(row[tpm_col]) # get the tpm
                    if gene_name not in gene_expression:
                        gene_expression[gene_name] = {}
                        gene_expression_avg[gene_name] = {}
                    gene_expression[gene_name][expname] = tpm
                    if condition_ph in gene_expression_avg[gene_name]:
                        gene_expression_avg[gene_name][condition_ph].append(tpm)
                    else:
                        gene_expression_avg[gene_name][condition_ph] = [tpm]
        else:
            print("no ",filepath)

# pfam parsing
pep2pfam = {}
pfamacc2name = {}
with open(pfamtbl,"rU") as pfam:
    l = csv.reader(pfam,delimiter=" ")
    for row in l:
        domain = row[0]
        domainacc = re.sub(r"\.\d+$","",row[1])
        gene   = row[2]
        if domainacc not in pfamacc2name:
            pfamacc2name[domainacc] = domain
        if gene not in pep2pfam:
            pep2pfam[gene] = []
        pep2pfam[gene].append(domainacc)

with open("Mtub.summary_table.tab","w") as rpt:
    resultout = csv.writer(rpt,delimiter="\t")
    header = ['GENE']
    header.extend(sorted(experiments))
    header.extend(["Pfam",'GO'])
    resultout.writerow(header)
    for gene in sorted(gene_expression):
        row = [gene]
        for n in sorted(experiments):
            row.append(gene_expression[gene][n])
        pfams = set()
        go    = set()
        if gene in pep2pfam:
            for pfam_acc in pep2pfam[gene]:
                pfams.add(pfamacc2name[pfam_acc])
                if pfam_acc in pfam2go:
                    for g in pfam2go[pfam_acc]:
                        go.add(g)
        row.append(",".join(sorted(pfams)))
        row.append(",".join(sorted(go)))

        resultout.writerow(row)

with open("Mtub.summary_mean_table.tab","w") as rpt:
    resultout = csv.writer(rpt,delimiter="\t")
    header = ['GENE']
    header.extend(sorted(experiments_avg))
    header.extend(["Pfam",'GO'])
    resultout.writerow(header)
    for gene in sorted(gene_expression_avg):
        row = [gene]
        for n in sorted(experiments_avg):
            row.append(mean(gene_expression_avg[gene][n]))
        pfams = set()
        go    = set()
        if gene in pep2pfam:
            for pfam_acc in pep2pfam[gene]:
                pfams.add(pfamacc2name[pfam_acc])
                if pfam_acc in pfam2go:
                    for g in pfam2go[pfam_acc]:
                        go.add(g)
        row.append(",".join(sorted(pfams)))
        row.append(",".join(sorted(go)))
        resultout.writerow(row)
