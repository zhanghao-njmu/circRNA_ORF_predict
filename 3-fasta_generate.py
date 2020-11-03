import os, sys, re
from pathlib import Path
maindir = str(Path(__file__).parent)
os.chdir(maindir)
with open("circ_rna_sequence.tsv", "w") as f:
    f.write("circ_id	gene	strand	#junction_reads	sequence\n")
    for line in open("./circRNA_rna_sequence.result", "r"):
        if re.search("^chr", line):
            block = line.split("\t")
            loc = block[0]
            gene = block[9]
            strand = block[10]
            junctionreads = block[4]
        else:
            seq = line.strip()
            content = loc + "\t" + gene + "\t" + strand + "\t" + junctionreads + "\t" + seq
            f.write(content + "\n")

with open("circ_protein_peptide_design.tsv", "w") as l:
    l.write(
        "location/model	nt_start	nt_end	juction location/circRNA length(nt)	protein length	juction location(aa)	juction aa	protein_sequence	peptide_start	peptide_end	peptide_length	peptide	circ_ID"+"\n"
        )
    with open("circ_protein_sequence.min11.fasta", "w") as h:
        with open("circ_protein_sequence.fasta", "w") as f:
            for line in open("./circ_protein_sequence.result", "r"):
                block = line.strip().split("\t")
                if re.search("^chr", line):
                    coord = block[0]
                    ntlen = int(float(block[1]))
                    num=0
                else:
                    num=num+1
                    loc = block[0]
                    start = block[1]
                    end = block[2]
                    aalen = block[3]
                    seq = block[4]
                    circle = int(start) // int(ntlen)
                    ntjunction = int(ntlen) * (circle + 1) - int(start)
                    aajunction = int(ntjunction / 3) + 1
                    aajunc = seq[aajunction - 1]
                    circid0=coord+loc+"pro"+str(num)
                    header = circid0 + " ntlen:" + str(
                        ntlen
                    ) + " ntloc:" + loc + " ntstart:" + start + " ntend:" + end + " aalen:" + aalen + " junctionaa:" + str(
                        aajunction) + aajunc
                    f.write(">" + header + "\n" + seq + "\n")
                    if len(seq) > 12:
                        seq_reduce = re.sub("_", "", seq)
                        h.write(">" + header + "\n" + seq_reduce + "\n")
                    else:
                        pass
                    
                    seqraw=seq
                    seq="_"*10+seq+"_"*10
                    circid=coord+loc+"pro"+str(num)
                    pepstart="none"
                    pepend="none"
                    for i in range(aajunction+10-2,0,-1):
                        aa=seq[i]
                        aa2=seq[i+1]
                        if (aa=="K" or aa=="R") and (aa2!="P"):
                            pepstart=i+2
                            break
                        if aa=="_":
                            break
                    for i in range(aajunction+10-1,len(seq),1):
                        aa=seq[i]
                        aa2=seq[i+1]
                        if (aa=="K" or aa=="R")and (aa2!="P"):
                            pepend=i+1
                            break
                        if aa2=="_":
                            break
                    pep="none"
                    pep_len=0
                    if pepstart!="none" and pepend!="none":
                        pep=seq[(pepstart-11):(pepstart-1)]+"--"+seq[(pepstart-1):(aajunction+10-1)]+"*"+seq[(aajunction+10-1):(pepend)]+"--"+seq[(pepend):(pepend+10)]
                        pep_len=pepend-pepstart+1
                    elif pepstart=="none" and pepend!="none":
                        pep=seq[0:10]+"--"+seq[10:(aajunction+10-1)]+"*"+seq[(aajunction+10-1):(pepend)]+"--"+seq[(pepend):(pepend+10)]
                        pep_len=pepend-10
                    elif pepstart!="none" and pepend=="none":
                        pep=seq[(pepstart-11):(pepstart-1)]+"--"+seq[(pepstart-1):(aajunction+10-1)]+"*"+seq[(aajunction+10-1):-10]+"--"+seq[-10:]
                        pep_len=len(seq)-10 -pepstart+1
                    if pepstart!="none":
                        pepstart=pepstart-10
                    if pepend!="none":
                        pepend=pepend-10
                    l.write(circid+"\t"+start+"\t"+end+"\t"+str(ntlen)+"\t"+aalen+"\t"+str(aajunction)+"\t"+aajunc+"\t"+seqraw+"\t"+str(pepstart)+"\t"+str(pepend)+"\t"+str(pep_len)+"\t"+pep+"\t"+coord+"\n")
                    

                        






