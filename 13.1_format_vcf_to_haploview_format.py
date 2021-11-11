import os, sys

# FAM1    FAM1M01    0    4    2    2
# FAM1    FAM1M01    0    4    2    2
# FAM1    FAM1F02    3    h    1    2
# FAM1    FAM1F02    3    h    1    2

# with 1=A, 2=C, 3=G, 4=T, 0=Missing, h=phase unsure

filename='Merged_AGPv4_sorted1.vcf'
outname='final.haps'
samplesize=435
filename = sys.argv[1] # input vcf, here FINAL.vcf 
outname = sys.argv[2] # name of the output, here final_.haps
samplesize = int(sys.argv[3]) # how many samples, here 782
dictallele={}
dictallele['A'] = '1'
dictallele['C'] = '2'
dictallele['G'] = '3'
dictallele['T'] = '4'

output = open(outname,'w')   
for indcounter in range(samplesize):
    curr_ind1 = []
    curr_ind2 = []
    for line in open(filename):
        if len(line) <= 1:
            continue


        if line.startswith('#CHROM'):
            a = line.strip().split('\t')
            curr_ind1.append(a[9+indcounter])
            curr_ind2.append(a[9+indcounter])
            continue      

  
        if line.startswith('#'):
            continue


        a = line.strip().split('\t')
        all1 = str(a[3])
        all2 = str(a[4])
        curr_geno = a[9+indcounter]
        if curr_geno.startswith('.'):
            curr_ind1.append('0')
            curr_ind2.append('0')


        elif curr_geno.startswith(('1/0','1|0')):
            curr_ind1.append('h')
            curr_ind2.append('h')


        elif curr_geno.startswith(('0/1','0|1')):
            curr_ind1.append('h')
            curr_ind2.append('h')


        elif curr_geno.startswith(('0/0','0|0')):    
            curr_ind1.append(dictallele[all1])
            curr_ind2.append(dictallele[all1])


        elif curr_geno.startswith(('1/1','1|1')):
            curr_ind1.append(dictallele[all2])
            curr_ind2.append(dictallele[all2])


        else:
            print ('There is a problem at ', a[0], ':', a[1], indcounter+1)





    output.write('FAM1\t'+ '\t'.join(map(str,curr_ind1))+'\n')
    output.write('FAM1\t'+ '\t'.join(map(str,curr_ind2))+'\n')





output.close()





