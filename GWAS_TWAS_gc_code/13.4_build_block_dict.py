import os,sys
#os.chdir("/Users/Meng/Desktop/RNAseq_temp/build_block_test/300K_window")

blockfile='4.res.GABRIELblocks'
outname='test.chr4.hmp.txt'
genofile='../HMP3_geno_chr4_istl1_region_3e+05.hmp.txt'
#outname='chr1_test.txt'

blockfile=sys.argv[1]
outname=sys.argv[2]
genofile=sys.argv[3]

with open(blockfile, 'r') as file:
    content = file.read()

dictallele={}
dictallele['1'] = 'A'
dictallele['2'] = 'C'
dictallele['3'] = 'G'
dictallele['4'] = 'T'


line=content.split("\n")
#line.append('Multiallelic')
#b=1
dicthaps={}
nm_line=len(line)
for a in range(nm_line):
    if line[a].startswith('BLOCK'):
        #b+=1
        info=line[a].split() # split by white space
        snp=info[3:]
        snp=tuple(snp)
        for c in range(1,10000):
            if line[a+c].startswith('Multiallelic') or len(line[a+c].strip())==0:
                break
            else:
                haps=line[a+c].split()
                hap=list(haps[0])
                for p in range(len(hap)):
                    hap[p]=dictallele[hap[p]]

                try: dicthaps[snp]+=[hap]
                except KeyError: dicthaps[snp]=[hap]



#file='test.txt'
#snps=('1','2','4')
#snps=tuple(snps)

snp_grp=dicthaps.keys()


output = open(outname,'w')
f = open(genofile).readline()
f1=f.strip().split('\t')
#f2=f1[:11]+['EndPos','nmSNP']+f1[11:]
f2=f1[:11]+['EndPos','nmSNP','block_number']+f1[11:]
for k in range(len(f2)-1):
    output.write(str(f2[k]) + '\t')
output.write(str(f2[-1]) + '\n')


block_nm=0 # generate block numbers, require the .hmp.txt genotype file has the same snp order as the input of haploview

for snps in snp_grp:
    print(snps)
    block_nm+=1
    print(block_nm)
    f = open(genofile)
    markers=[]
    for i, line in enumerate(f):
        if str(i) in snps:
            markers.append(line.strip().split('\t'))

    nm_col=len(markers[0])
    nm_snp=len(markers)
    dictblock={}
    for m in range(11,nm_col):
        for n in range(nm_snp):
            try: dictblock[m]+= markers[n][m]
            except KeyError: dictblock[m]=[markers[n][m]]

    #geno_hap=markers[0][:11]+[str(markers[-1][3])]+[str(len(markers))]
    geno_hap=markers[0][:11]+[str(markers[-1][3])]+[str(len(markers))]+[str(block_nm)]
    taxa=list(dictblock.keys())
    taxa.sort()
    for l in taxa:
        if dictblock[l] in dicthaps[snps]:
            haptype=dicthaps[snps].index(dictblock[l])
            geno_hap.append(haptype)
        else:
            geno_hap.append('NA')

    #print(geno_hap)
    for p in range(nm_col+2):
        output.write(str(geno_hap[p]) + '\t')
    output.write(str(geno_hap[nm_col+2]) + '\n')



output.close()
