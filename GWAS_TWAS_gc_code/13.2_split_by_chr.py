import os, sys, collections
from collections import defaultdict


hapsname = sys.argv[1] # name of the whole genome "haps" file (e.g. FINAL.haps)
snpname = sys.argv[2] # name of the chrpos file (e.g FINAL.info)

n = 1 # this counts the number of SNPs per chrm
m = 0 # this counts where we are in the file
prev_chrm = 'notstartedyet'
#sumcheckm = []
for line in open(snpname):
    if len(line) <= 1:
        continue
    a = line.strip().split('\t')
    curr_chrm = a[0]
    pos = a[1]
    #
    if prev_chrm != curr_chrm:
        # when we change chromosome, we subset the amount of SNPs we have read
        if m > 0: # that should not be true when opening the first chromosome 
            cmd = 'cut -f1-2,' + str(m+2) + '-' + str(m+n+1) + ' ' + hapsname + ' > ' + prev_chrm + '.haps'     
            os.system(cmd)
            print(prev_chrm, m, n)
        m += n # remember how many snp we have reached
        n = 0 # and reset n
        #sumcheckm.append(m)
        if 'outputinfo' in globals():
            if not outputinfo.closed:
                outputinfo.close()
        outputinfo = open(curr_chrm + '.info', 'w')
    #
    #outputinfo.write('marker' + str(n) + '\t' + str(pos) + '\n')
    outputinfo.write(str(curr_chrm)+"_"+str(pos) + '\t' + str(pos) + '\n')
    n += 1
    prev_chrm = curr_chrm
#hapsfile.close()
outputinfo.close()
#sumcheckm.append(m+n)
# write the last one (we don't have a change of chr for the last block since we hit the end of the file)
cmd = 'cut -f1-2,' + str(m+2) + '-' + str(m+n+1) + ' ' + hapsname + ' > ' + prev_chrm + '.haps'     
os.system(cmd)
