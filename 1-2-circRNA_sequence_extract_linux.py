#!usr/bin/python3
from pathlib import Path
import os, sys, re,shutil,multiprocessing,time
maindir=str(Path(__file__).parent)
os.chdir(maindir)
split_num=88


def task(split_character,split_num):
#    split_character=0
    splitdir=maindir+"/tmp/split{split_character}".format(split_character=split_character)

    print("Step1:   extract circRNA containing exon in CIRI.result")
    with open(splitdir+"/tmp_circ_contain_exon.tmp", "w") as circ_select:
        
        for line in open(splitdir+"/split{split_character}.result".format(split_character=split_character), "r"):
            if re.search("exon", line):
                circ_select.write(line)


    print("Step2:   extract exon location in gtf for each circRNA")
    gtf = open(maindir+"/Homo_sapiens.GRCh37.75.gtf", "r")
    gene = gtf.readlines()
    exon_list = open(splitdir+"/tmp_exon_list.tmp", "w")
    with open(splitdir+"/tmp_circ_contain_exon.tmp", "r") as circ_select:
        for i2 in circ_select:
            circ = i2.strip().split('\t')
            #chr_number
            chromo = circ[1].split('chr')[1]
            exon_list.write(i2)
            for line in gene:
                block = line.split('\t')
                if chromo == block[0] and block[2] == "exon" and int(circ[2]) <= int(
                        block[3]) and int(circ[3]) >= int(block[4]) and str(circ[10]) == str(block[6]):
                    #In the same chr,if CIRI_start <= exon_start and CIRI_end >= exon_end, then write the exon line.
                    infor = line
                    exon_list.write(infor)
    exon_list.close()

    print("Step2-2:   filter out exon boundary in circRNA")
    with open(splitdir+"/exon_filter.result","w")as f:
        dic={}
        for line in open(splitdir+"/tmp_exon_list.tmp","r"):
            j = line.split('\t')
            if j[2]!="exon":
                pass
            else:
                trans_id=re.findall("(?<=transcript_id \")(\w+)(?=\"\;)",line)[0]
                if trans_id not in dic.keys():
                    chrom=j[0]
                    exon_start=j[3]
                    exon_end=j[4]
                    dic[trans_id]="chr"+chrom+":"+exon_start+"|"+exon_end
                else:
                    if j[3]<exon_start:
                        exon_start=j[3]
                    if j[4]>exon_end:
                        exon_end=j[4]
                    dic[trans_id]="chr"+chrom+":"+exon_start+"|"+exon_end


        writenum=0
        for line in open(splitdir+"/tmp_exon_list.tmp","r"):
            j = line.split('\t')
            if j[2]!="exon":
                circ_id=str(j[0])
                circ_info=line
                writenum=0
            else:
                trans_id=re.findall("(?<=transcript_id \")(\w+)(?=\"\;)",line)[0]
                if str(dic[trans_id])==circ_id:
                    if writenum==0:
                        f.write(circ_info)
                        writenum=1
                    f.write(line)
        f.write("let	the	next	step	run	to	the	last	sequence")


    print("Step3:   extract & merge exon sequence from genome")
    exon_list2 = open(splitdir+"/exon_filter.result", "r")
    rna_seq = open(splitdir+"/circRNA_rna_sequence.out", "w")
    for i in exon_list2:
        j=i.split('\t')
        if j[2]!="exon":    ###start to balance the previous circRNA sequence
            if 'seq' in locals().keys():
            #     l={}
            #     for key in seq:
            #         l[len(seq[key])]=seq[key]   ### l is a dirctionary containing all sequence length and corresponding sequence:{sequence_length:sequence}
            #     maxchang=max(l.keys())  ### maxchang is the longest length in sequence_length
            #     results=l[maxchang] ### results is the corresponding longest sequence
            #     rna_seq.write(results+'\n')
            # else:
            #     pass
                sequence_dedup={}
                for key in seq.keys():
                    if seq[key] not in sequence_dedup.keys():
                        sequence_dedup[seq[key]]=key    ### {sequence:transcript_ID}
                    else:
                        sequence_dedup[seq[key]]=sequence_dedup[seq[key]]+";"+key   ### {sequence:transcript_ID1;transcript_ID2}
                for uni_seq in sequence_dedup.keys():
                    trans_id_list=sequence_dedup[uni_seq]
                    results=i_previous.replace(j_previous[0],j_previous[0]+"__"+trans_id_list)
                    rna_seq.write(results)
                    rna_seq.write(uni_seq+"\n")
            trans_id=[]
            seq={}
            i_previous=i
            j_previous=j
        elif j[2]=="exon":  ### start to calculate the next circRNA sequence
            k=re.findall("(?<=transcript_id \")(\w+)(?=\"\;)",i)[0]  ### e.g.: gene_id "ENSMUSG00000025917" gene_version "9" transcript_id "ENSMUST00000027050" transcript_version "9" exon_number "2" gene_name "Cops5" gene_source "ensembl_havana" gene_biotype "protein_coding" transcript_name "Cops5-201" transcript_source "ensembl_havana" transcript_biotype "protein_coding" tag "CCDS" ccds_id "CCDS14818" exon_id "ENSMUSE00001218429" exon_version "1" tag "basic" transcript_support_level "1"
            file=open(maindir+"/Chromosomes/"+j[0]+".fa","r")
            alseq_1=file.readline() ### skip the first annotation line
            alseq = file.read() ### read all sequence in a chromosome
            alseq=alseq.replace("\n","")    ###put all sequence in the same line
            file.close()
            if j[6]=="+":
                circseq=alseq[(int(j[3])-1):(int(j[4]))]
            else:
                circseq=alseq[(int(j[3])-1):(int(j[4]))]
                circseq=circseq[::-1]
                circseq=circseq.replace("A","X")
                circseq=circseq.replace("a","X")
                circseq=circseq.replace("T","A")
                circseq=circseq.replace("t","A")
                circseq=circseq.replace("X","T")
                circseq=circseq.replace("C","Y")
                circseq=circseq.replace("c","Y")
                circseq=circseq.replace("G","C")
                circseq=circseq.replace("g","C")
                circseq=circseq.replace("Y","G")
            if k in trans_id : ### k is the transcript_ID e.g. ENSMUST00000027050
                seq[k]=seq[k]+str(circseq)
            else:
                trans_id.append(k)
                seq[k]=str(circseq)
    exon_list2.close()
    rna_seq.close()

    print("Step4:   circRNA sequence *4 ")
    rna_seq = open(splitdir+"/circRNA_rna_sequence.out","r")
    rna_seq_triple = open(splitdir+"/circRNA_rna_sequence_triple.out","w")
    for i in rna_seq:
        j=i.split('\t')
        if len(j)>3:
            rna_seq_triple.write(i)
        else:
            seq=j[0].strip('\n')
            seq=seq*4
            rna_seq_triple.write(seq+"\n")
    rna_seq.close()
    rna_seq_triple.close()

    print("Step5:   circRNA translate to protein ")
    def translate_dna(sequence):    ### obtain all ORF in a circRNA sequence
        codontable={                                  
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
        proteinsequence={}
        head={}
        end={}
        lenth=len(sequence)
        thirdlenth=lenth/4
        number={}
        for x in range(3):  ### 3 start position options
            proteinsequence[x]={}
            n=0 ### nt position
            i=0 ### ORF count
            head[x]={}  ### ATG position in RNA sequence
            end[x]={}   ### stop position codon in RNA sequence
            while (n+3<=lenth):  ###    n+3<lenth
                proteinsequence[x][i]=''
                if sequence[n:n+3]=='ATG': ###  index: n|n+1|n+2
                        head[x][i]=x+n+1
                        while (sequence[n:n+3] in codontable):
                            proteinsequence[x][i]+=codontable[sequence[n:n+3]]
                            if (codontable[sequence[n:n+3]]=='_'):
                                end[x][i]=x+n+3
                                i=i+1
                                break
                            n=n+3
                n=n+3
            number[x]=i
            sequence=sequence[1:]+sequence[0]   ### move the first nucleotide to the end of sequence
        return head,end,proteinsequence,number,thirdlenth

    def writing(head,end,proteinsequence,number,circ_header,thirdlenth):
        newfile=open(splitdir+'/circ_protein_sequence.out','a')
        newfile.write(str(circ_header)+'\t'+str(thirdlenth)+'\n')
        aa=[]
        for x in range(3):  ### 3 start position options
            for i in range(number[x]):
                if (proteinsequence[x][i][-1]=='_'):
                    if (int(head[x][i]/thirdlenth)<int(end[x][i]/thirdlenth)):  ### ATG position in circle[1] and stop codon position in another circle[2] 
                            l=str(proteinsequence[x][i])
                            if not str(l) in aa:    ### deduplication
                                aa.append(str(l))
                                if int(head[x][i])<int(thirdlenth) and int(end[x][i])>int(thirdlenth):
                                    newfile.write('+'+str(x)+'\t'+str(head[x][i])+'\t'+str(end[x][i])+'\t'+str(len(l))+'\t'+proteinsequence[x][i]+'\n')
                                elif int(head[x][i])>int(thirdlenth) and int(head[x][i])<int(thirdlenth)*2 and int(end[x][i])>int(thirdlenth)*2:
                                    newfile.write('+'+str(x)+'\t'+str(head[x][i])+'\t'+str(end[x][i])+'\t'+str(len(l))+'\t'+proteinsequence[x][i]+'\n')
                                elif int(head[x][i])>int(thirdlenth)*2 and int(head[x][i])<int(thirdlenth)*3 and int(end[x][i])>int(thirdlenth)*3:
                                    newfile.write('+'+str(x)+'\t'+str(head[x][i])+'\t'+str(end[x][i])+'\t'+str(len(l))+'\t'+proteinsequence[x][i]+'\n')
        newfile.close()

    for line in open(splitdir+'/circRNA_rna_sequence_triple.out','r'):
        if re.search("^chr",line):
            circ_header=line.split("\t")[0]
        else:
            sequence=line.strip()
            head,end,proteinsequence,number,thirdlenth=translate_dna(sequence)
            writing(head,end,proteinsequence,number,circ_header,thirdlenth)


if __name__ == "__main__":
    print("cpu_count: "+str(multiprocessing.cpu_count()))
    multiprocessing.freeze_support()  # Windows add,avoid RuntimeError
    pool = multiprocessing.Pool()
    for i in range(split_num):
        pool.apply_async(task, args=(i,split_num))
    pool.close()
    pool.join()











