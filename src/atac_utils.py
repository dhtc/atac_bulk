from __future__ import division
import os
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
def to_log(msg, log_filepath):
	with open(log_filepath, 'a') as log_file:
		log_file.write(msg + '\n')

def ref_size(genome_sizes_filepath):
	genome_sizes={}
	with open(genome_sizes_filepath) as genome_sizes_file:
		for line in genome_sizes_file:
			chrom, size = line.strip().split('\t')
			genome_sizes[chrom] = int(size)
	return genome_sizes


def within_genome_boundary(chrom, check_pos, genome_sizes_dict):
	if check_pos < 1:
		return False
	if check_pos > genome_sizes_dict[chrom]:
		return False
	return True

def Fragdistribution(Inbam,outfig,threads):
	sam=os.popen("samtools view -@ %s %s"%(threads,Inbam))
	fragL=[]
	for line in sam:
		items=line.rstrip("\n").split('\t')
		flag=items[1]
		flag_length=abs(int(items[8]))
		if (int(flag)==99 or int(flag)==163):
			fragL.append(flag_length)
	count=Counter(fragL)
	total_frag=sum(count.values())
	with open(outfig.replace('_fragLen.pdf','_fragLen.txt'),'w')as f:
		for k,c in dict(count).items():
			f.write(str(k)+'\t'+str(c)+'\n')
	fig=plt.figure(figsize=(10.0,4.0))
	plt.plot(count.keys(),np.array(count.values())/float(total_frag),'r')
	plt.xlabel('Read length')
	plt.ylabel('Read counts %')
	fig.savefig(outfig)
	plt.close(fig)
	return



########################
#shift for peak calling#
########################
def shift_reads_bedpe(input_bed_filepath, output_bed_filepath, genome_sizes_filepath,log_filepath, minus_strand_shift, plus_strand_shift):
	genome_sizes=ref_size(genome_sizes_filepath)
	output_bed = open(output_bed_filepath, 'w')
	with open(input_bed_filepath) as input_bed:
		for line in input_bed:
			record=line.strip('\n').split('\t')
			if len(record)<7:
				to_log('Mal-formatted line in file, skipping:', log_filepath)
				to_log(line.strip(), log_filepath)
				continue;
			else:
				chrom, start, end, name, score, strand1, strand2 = record[0], int(record[1]), int(record[2]), record[3], int(record[4]), record[5], record[6].rstrip('\n')
				if strand1 == '+' and strand2 == '-':
					new_start = start + plus_strand_shift
					new_end = end + plus_strand_shift
				elif strand1 == '-' and strand2 == '+':
					new_start = start + minus_strand_shift
					new_end = end + minus_strand_shift
				else:
					to_log('Mal-formatted line in file, skipping:', log_filepath)
					to_log(line.strip(), log_filepath)
					continue
				if (not within_genome_boundary(chrom, new_start, genome_sizes)) or (not within_genome_boundary(chrom,new_end,genome_sizes)):
					to_log('Shift end site beyond chromosome boundary:', log_filepath)
					to_log(line.strip(), log_filepath)
					continue
				output_bed.write('\t'.join([chrom, str(new_start), str(new_end)] + record[3:]) + '\n')
		output_bed.flush()
		output_bed.close()
###############
#visualization#
###############
def extend_bed(inbed,outbed,ref_size_file,extend,log_filepath):
	extend=int(extend)
	out=open(outbed,'w')
	genome_sizes= ref_size(ref_size_file)
	with open(inbed,'r')as inbed:
		for line in inbed:
			record=line.strip('\n').split('\t')
			if len(record)<7:
				to_log('Mal-formatted line in file, skipping:', log_filepath)
				to_log(line.strip(), log_filepath)
				continue;
			else:
				chrom, start, end, name, score, strand1, strand2 = record[0], int(record[1]), int(record[2]), record[3], int(record[4]), record[5], record[6].rstrip('\n')
				if strand1=='+' and strand2=='-':
					new_start=start-extend
					new_end=end+extend
				if strand1=='-' and strand2=='+':
					new_start=start-extend
					new_end=end+extend
				if (not within_genome_boundary(chrom, new_start, genome_sizes)) or (not within_genome_boundary(chrom,new_end,genome_sizes)):
					to_log('extend end site beyond chromosome boundary:', log_filepath)
					to_log(line.strip(), log_filepath)
					continue
				out.write('\t'.join([chrom, str(new_start), str(new_end)] + record[3:]) + '\n')
		out.flush()
		out.close()

def PerbasebedGraph(InbedGraph):
	tmp=InbedGraph+'.tmp'
	tmpf=open(tmp,'w')
	lines=open(InbedGraph,'r')
	for line in lines:
		items=line.split('\t')
		for i in range(int(items[2])-int(items[1])):
			tmpf.write(items[0]+'\t'+str(int(items[1])+i)+'\t'+str(int(items[1])+i+1)+'\t'+str(items[3]))
	os.system('mv %s %s'%(tmp,InbedGraph))
	return

def Normbedgraph(inbdg,outbdg,Readlength=50,Totalcount=10000000):
	infile=open(inbdg,'r')
	sumOfRead=0
	for line in infile:
		items=line.rstrip('\n').split()
		sumOfRead=sumOfRead+int(items[3])*(int(items[2])-int(items[1]))
	infile.close()
	sumOfRead=abs(sumOfRead)
	infile=open(inbdg,'r')
	outfile=open(outbdg,'w')
	rawReadLength=int(Readlength)
	for line in infile:
		items=line.rstrip('\n').split()
		value=float(items[3])/(sumOfRead/rawReadLength)*int(Totalcount)
		outfile.write('\t'.join([items[0],items[1],items[2],str(value)])+'\n')
	infile.close()
	outfile.close()
	return

