# -*- coding: UTF-8 -*-
from __future__ import division
import os
import subprocess
import re
import time
from chunkypipes.components import Software, Parameter, Redirect, Pipe, BasePipeline
import gzip
import time
import math
import numpy as np
import scipy as sp
#import matplotlib.pyplot as plt
import collections
import pandas as pd
from atac_utils import to_log,ref_size,within_genome_boundary,shift_reads_bedpe,extend_bed,PerbasebedGraph,Normbedgraph,footprint,Fragdistribution
from atac_load_cfg import load_cfg
#matplotlib.use('Agg')


MINUS_STRAND_SHIFT = -5
PLUS_STRAND_SHIFT = 4


class Pipeline(BasePipeline):
	def description(self):
		return ['Pipeline for analyzing ATAC-seq data']
	def configure(self):
		return{
			'fastp':{
				'path': 'Full path to fastp executable'
			},
			'bowtie2':{
				'path': 'Full path to bowtie2 executable',
				'index-dir': 'Directory of the bowtie2 reference index [Ex. /path/to/index/genome basename]'
			},
			'fastqc':{
				'path': 'Full path to FastQC'
			},
			'samtools':{
				'path': 'Full path to samtools'
			},
			'picard':{
				'path': 'Full path to Picard [Ex. java -jar /path/to/picard.jar]'
			},
			'bedtools':{
				'path': 'Full path to bedtools >=2.25.0',
				'blacklist-bed': 'Full path to the BED of blacklisted genomic regions,or None',
				'genome-sizes': 'Full path to a genome sizes file'
			},
			'macs2':{
				'path': 'Full path to MACS2'
			},
			'run_spp_nodups':{
				'path':'Full path to run_spp_nodups'
			},
			'chipqc':{
				'path':'/Full/path/to/Rscript /Full/path/to/chipqc.R',
				'sample':'/Full/path/to/sample',
				'annotation':'annotation,eg:hg38'
			},
			'bdg2bw':{
				'path':'/Full/path/to/bedGraphToBigWig'
			},
			'bed_Clip':{
				'path':'/Full/path/to/bed_Clip'
			},
			'idr':{
				'path':'Full/path/to/idr'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('--reads-cfg', required=True, help='configure file,tab separate,1st column for lib name,2nd column for full path of read1,3rd for read2,4th for group,Required', action='store')
		parser.add_argument('--output', required=True,help='output directory,Required')
		parser.add_argument('--step', default=0,type=str,help='1 for fq clean,2 for align,3 for align filter,4 for shift,5 for peak calling,6 for cross-corelation,7 for visualization,8 for idr,9 for difference,10 for ReadLenDistribution,11 for pygenometracks,0 for quality control')
		parser.add_argument('--ref-genome',required=True,help='/path/to/reference genome basename for bowtie2,Required')
		parser.add_argument('--threads',default=1,type=int)
		parser.add_argument('--rm-tmp',default='n',choices=['y','n'],help='need to rm temporary files,y|n?')
		parser.add_argument('--effective-genome-size',default='hs',type=str,help='effective genome size for macs2 peak calling,default human')
		parser.add_argument('--broad',action='store_true',help='call broad peaks?')
		parser.add_argument('--visual-method',required=True,help='visualization method,llr or pileup',type=str,choices=['llr','pileup'])
		parser.add_argument('--extend',required=True,help='extend for both side',type=str,default='250')
		parser.add_argument('--idr-thresh',default='0.05',type=str,help='Report statistics for peaks with a global idr below this value but return all peaks with an idr below --idr.Default: 0.05')
		parser.add_argument('--make-tracks-files',action='append',type=str,nargs='*',help='files to visualize using pygenometracks')
		parser.add_argument('--tracks',type=str,default='tracks.ini',help='pygenometracks configure file')
		parser.add_argument('--region',type=str,nargs='*',help='pygenometracks region,eg: chr2:10,000,000-11,000,000',action='append')
		parser.add_argument('--outFileName',type=str,nargs='*',help='outfilename for pygenometracks',action='append')
		return parser





	def run_pipeline(self, pipeline_args, pipeline_config):
		
		libs,read1s,read2s,groups=load_cfg(pipeline_args['reads_cfg'])
		
		# Instantiate variable from argparse
		if pipeline_args['ref_genome']:
			ref=pipeline_args['ref_genome']
		else:
			ref=pipeline_config['bowtie2']['index-dir']
		threads=str(pipeline_args['threads'])
		effective_genome_size=pipeline_args['effective_genome_size']
		chipqc_sample=pipeline_config['chipqc']['sample']
		chipqc_annotation=pipeline_config['chipqc']['annotation']
		chipqc_blacklist=pipeline_config['bedtools']['blacklist-bed']
		visual_method=pipeline_args['visual_method']
		extend=pipeline_args['extend']
		idr_thresh=pipeline_args['idr_thresh']
		tracks_cfg_file=pipeline_args['tracks']
		tracks_regions=pipeline_args['region']
		tracks_outfiles=pipeline_args['outFileName']
		make_tracks_files=pipeline_args['make_tracks_files']
		diffbind_sample=pipeline_config['diffbind']['sample']
		
		#create dirs
		output_dir = os.path.abspath(pipeline_args['output'])
		logs_dir = os.path.join(output_dir, 'logs')
		clean_dir=os.path.join(output_dir,'clean')
		align_dir=os.path.join(output_dir,'align')
		align_clean_dir=os.path.join(output_dir,'align_clean')
		align_clean_tmp_dir=os.path.join(align_clean_dir,'tmp')
		bed_dir=os.path.join(output_dir,'bed')
		qc_dir=os.path.join(output_dir,'qc')
		peaks_dir=os.path.join(output_dir,'peaks')
		bdgs_dir=os.path.join(output_dir,'bdgs')
		bws_dir=os.path.join(output_dir,'bws')
		diff_dir=os.path.join(output_dir,'diff')
		atac_logfile=os.path.join(logs_dir,'atac.log')
		
		
		step =pipeline_args['step'].split(',')
		


		# Create directories

		subprocess.call(['mkdir','-p',output_dir,logs_dir,clean_dir,align_dir,align_clean_dir,bed_dir,qc_dir,align_clean_tmp_dir,peaks_dir,bdgs_dir,bws_dir])

		#Keep list of items to delete
		staging_delete = []

		#Instatiate software instances
		fastp = Software('fastp', pipeline_config['fastp']['path'])
		fastqc = Software('fastqc', pipeline_config['fastqc']['path'])
		bowtie2 = Software('bowtie2', pipeline_config['bowtie2']['path'])
		samtools_view = Software('samtools view', pipeline_config['samtools']['path'] + ' view')
		samtools_flagstat = Software('samtools flagstat', pipeline_config['samtools']['path'] + ' flagstat')
		samtools_index = Software('samtools index', pipeline_config['samtools']['path'] + ' index')
		samtools_sort = Software('samtools sort', pipeline_config['samtools']['path'] + ' sort')
		picard_mark_dup = Software('Picard MarkDuplicates', pipeline_config['picard']['path'] + ' MarkDuplicates')
		bedtools_bamtobed = Software('bedtools bamtobed', pipeline_config['bedtools']['path'] + ' bamtobed')
		bedtools_sort = Software('bedtools sort', pipeline_config['bedtools']['path'] + 'sort')
		bedtools_merge = Software('bedtools merge', pipeline_config['bedtools']['path'] + ' merge')
		bedtools_intersect = Software('bedtools intersect', pipeline_config['bedtools']['path'] + ' intersect')
		bedtools_genomecov=Software('bedtools genomecov',pipeline_config['bedtools']['path']+' genomecov')
		macs2_callpeak = Software('macs2 callpeak', pipeline_config['macs2']['path'] + ' callpeak')
		macs2_bdgcmp=Software('macs2 bdgcmp',pipeline_config['macs2']['path']+' bdgcmp')
		macs2_bdgdiff=Software('macs2 bdgdiff',pipeline_config['macs2']['path']+' bdgdiff')
		run_spp_nodups=Software('run_spp_nodups',pipeline_config['run_spp_nodups']['path'])
		chipqc=Software('chipqc',pipeline_config['chipqc']['path'])
		bdg2bw=Software('bdg2bw',pipeline_config['bdg2bw']['path'])
		bed_Clip=Software('bed_Clip',pipeline_config['bed_Clip']['path'])
		idr=Software('idr',pipeline_config['idr']['path'])
		pygenometracks=Software('pygenometracks',pipeline_config['pygenometracks']['path'])
		make_tracks=Software('make_tracks',pipeline_config['make_tracks']['path'])
		diffbind=Software('diffbind',pipeline_config['diffbind']['path'])

		if '110120' in step:
			print('test passed')
		if  '1' in step: #trim by fastp and fastqc
			trimmed_read1s=[]
			trimmed_read2s=[]
			for i, lib in enumerate(libs):
				read1=read1s[i]
				read2=read2s[i]
				trimmed_read1_filename = os.path.join(output_dir,'clean',lib,'{}_1.fq'.format(lib))
				trimmed_read2_filename = os.path.join(output_dir,'clean',lib,'{}_2.fq'.format(lib))
				subprocess.call(['mkdir','-p',os.path.join(output_dir,'clean',lib)])
				trimmed_read1s.append(trimmed_read1_filename)
				trimmed_read2s.append(trimmed_read2_filename)
				
				fastp.run(
					Parameter('-i',read1),
					Parameter('-I',read2),
					Parameter('-o',trimmed_read1_filename),
					Parameter('-O',trimmed_read2_filename),
					Parameter('-l','40'),
					Parameter('-q','30'),
					Parameter('--detect_adapter_for_pe'),
					Parameter('--adapter_sequence=CTGTCTCTTATACACATCT --adapter_sequence_r2=AGATGTGTATAAGAGACAG'),
					Parameter('-c'),
					Parameter('-3'),
					Parameter('-w',threads),
					Parameter('-j',os.path.join(qc_dir,lib+'.json')),
					Parameter('-h',os.path.join(qc_dir,lib+'.html')),
					Redirect(stream=Redirect.STDOUT_APPEND, dest=os.path.join(logs_dir, 'fastp.std.log')),
					Redirect(stream=Redirect.STDERR_APPEND, dest=os.path.join(logs_dir, 'fastp.stderr.log'))
					
				)
				fastqc.run(
					Parameter('-t',threads),
					Parameter('-o',qc_dir),
					Parameter(trimmed_read1_filename),
					Parameter(trimmed_read2_filename)
				)
				to_log(read1+' '+read2+' has been trimmed&fastqc '+str(time.time()),atac_logfile)


		if '2' in step:#align
			if 'trimmed_read1s' not in vars():
				trimmed_read1s=[]
				trimmed_read2s=[]
				for i,lib in enumerate(libs):
					trimmed_read1_filename = os.path.join(output_dir,'clean',lib,'{}_1.fq'.format(lib))
					trimmed_read2_filename = os.path.join(output_dir,'clean',lib,'{}_2.fq'.format(lib))
					trimmed_read1s.append(trimmed_read1_filename)
					trimmed_read2s.append(trimmed_read2_filename)
			for i, lib in enumerate(libs):
				read1=trimmed_read1s[i]
				read2=trimmed_read2s[i]
				sam_out=os.path.join(align_dir,lib+'.sam')
				bowtie2.run(
					Parameter('-p',threads),
					Parameter('-q'),
					Parameter('-t'),
					Parameter('-N','1'),
					Parameter('-L','25'),
					Parameter('-X','2000'),
					Parameter('--no-mixed'),
					Parameter('--no-discordant'),
					Parameter('-x',ref),
					Parameter('-1',read1),
					Parameter('-2',read2),
					Parameter('-S',sam_out),
					Redirect(stream=Redirect.STDOUT_APPEND, dest=os.path.join(logs_dir, 'alignment.std.log')),
					Redirect(stream=Redirect.STDERR_APPEND, dest=os.path.join(logs_dir, 'alignment.stderr.log'))
				)
				samtools_view.run(
					Parameter('-bh'),
					Parameter('-o',os.path.join(align_dir,lib+'.bam')),
					Parameter('-@',threads),
					Parameter(sam_out)
				)
				staging_delete.append(read1)
				staging_delete.append(read2)
				staging_delete.append(sam_out)
				to_log(read1+' '+read2+' has been aligned '+str(time.time()),atac_logfile)
		if '3' in step:#alignment filter
			align_quality_threshold='10'
			bams=[os.path.join(align_dir,lib+'.bam') for lib in libs]
			outbams=[]
			for i,bam in enumerate(bams):
				DechrMbam=os.path.join(align_clean_tmp_dir,libs[i]+'.pe.q'+align_quality_threshold+'.sort.bam')
				os.system("samtools view -h -f 0x2 -q %s -@ %s %s|grep -v \"chrM\"|samtools sort -O bam -o %s -@ %s"%(align_quality_threshold,threads,bam,DechrMbam,threads))
				Dedupbam=os.path.join(align_clean_tmp_dir,libs[i]+'.rmdup.bam')
				os.system("picard MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true > %s" %(DechrMbam,Dedupbam,os.path.join(logs_dir,libs[i]+'.Picard_Metrics_unfiltered_bam.txt'),os.path.join(logs_dir,libs[i]+'.Picard.log')))
				sortout=os.path.join(align_clean_tmp_dir,libs[i]+'.rmchrm.dedup.sorted.bam')
				os.system("samtools sort -o %s -@ %s %s"%(sortout,threads,Dedupbam))
				outbams.append(sortout)
				staging_delete.append(DechrMbam)
				staging_delete.append(Dedupbam)
				to_log(bam+'has been filtered '+str(time.time()),atac_logfile)
		if '4' in step:#shift bed
			if 'outbams' not in vars():
				outbams=[]
				for i,lib in enumerate(libs):
					outbams.append(os.path.join(align_clean_tmp_dir,lib+'.rmchrm.dedup.sorted.bam'))
			processed_bedpe_to_beds=[]
			for i,bam in enumerate(outbams):
				# Generate filename for final processed BAM and BED
				rm_black_list_bam = os.path.join(align_clean_dir, '{}.bam'.format(libs[i]))
				sorted_for_PE_bam = os.path.join(align_clean_tmp_dir, '{}.sorted_for_PE.bam'.format(libs[i]))
				unshifted_bedpe = os.path.join(align_clean_tmp_dir,'{}.unshifted_bedpe'.format(libs[i]))
				processed_bedpe_to_bed=os.path.join(bed_dir,'{}.bed'.format(libs[i]))
				processed_bedpe_to_beds.append(processed_bedpe_to_bed)
				
				staging_delete.append(bam)
				staging_delete.append(sorted_for_PE_bam)
				staging_delete.append(unshifted_bedpe)
				staging_delete.append(rm_black_list_bam+'.bai')
				
				# Remove blacklisted genomic regions
				if pipeline_config['bedtools']['blacklist-bed'] != 'None':
					bedtools_intersect.run(
						Parameter('-v'),
						Parameter('-abam', bam),
						Parameter('-b', pipeline_config['bedtools']['blacklist-bed']),
						Parameter('-f', '0.5'),
						Redirect(stream=Redirect.STDOUT, dest=rm_black_list_bam)
					)
				else:
					os.system('cp %s %s'%(bam,rm_black_list_bam))
					rm_black_list_bam = bam
				# Generate index for processed BAM
				samtools_index.run(
					Parameter(rm_black_list_bam)
				)
				# Convert BAM to BEDPE, with specific quality and only properly paired reads, sorted by name
				samtools_view.run(
					Parameter('-uf', '0x2'),
					Parameter('-F', '1548'),
					Parameter('-q', '30'),
					Parameter('-@',threads),
					Parameter(rm_black_list_bam),
					Pipe(
						samtools_sort.pipe(
							Parameter('-n'),
							Parameter('-o',sorted_for_PE_bam),
							Parameter('-@',threads),
							Parameter('-')
						)
					)
				)

				# convert bam to BEDPE
				bedtools_bamtobed.run(
					Parameter('-i', str(sorted_for_PE_bam)),
					Parameter('-bedpe'),
					Redirect(stream=Redirect.STDOUT, dest=unshifted_bedpe)
				)
				unshifted_bedpe_to_bed = open(align_clean_tmp_dir+'/'+'{}.unshifted_bedpe_to_bed'.format(libs[i]), 'w')
				with open(unshifted_bedpe) as convertToBed:
					for line in convertToBed:
						chrpos1, start1, end1, chrpos2, start2, end2, name, score, strand1, strand2=line.split('\t')
						bedformat=[chrpos1, start1, end2, name, score, strand1, strand2.rstrip('\n')]
						unshifted_bedpe_to_bed.write('\t'.join(bedformat)+'\n')
				unshifted_bedpe_to_bed.close()

				staging_delete.append(align_clean_tmp_dir+'/'+'{}.unshifted_bedpe_to_bed'.format(libs[i]))

				#shift reads for 9bp,+4-5
				shift_reads_bedpe(
					input_bed_filepath=align_clean_tmp_dir+'/'+'{}.unshifted_bedpe_to_bed'.format(libs[i]),
					output_bed_filepath=processed_bedpe_to_bed,
					log_filepath=os.path.join(logs_dir,'shift_reads.logs'),
					genome_sizes_filepath=pipeline_config['bedtools']['genome-sizes'],
					minus_strand_shift=MINUS_STRAND_SHIFT,
					plus_strand_shift=PLUS_STRAND_SHIFT
				)


		if '0' in step:#chipqc
			for i,lib in enumerate(libs):
				samtools_index.run(
					Parameter('-@',threads),
					Parameter(os.path.join(align_clean_dir, '{}.bam'.format(libs[i]))),
					Parameter(os.path.join(align_clean_dir,'%s.bam.bai'%libs[i]))
				)
			if chipqc_blacklist != 'None':
				chipqc.run(
					Parameter(chipqc_sample),
					Parameter(chipqc_annotation),
					Parameter(os.path.join(logs_dir,'chipqc.Rdata')),
					Parameter(chipqc_blacklist),
					Redirect(stream=Redirect.STDOUT_APPEND,dest=os.path.join(logs_dir,'chipqc_std.log')),
					Redirect(stream=Redirect.STDERR_APPEND,dest=os.path.join(logs_dir,'chipqc_stderr.log'))
				)
			else:
				chipqc.run(
					Parameter(chipqc_sample),
					Parameter(chipqc_annotation),
					Parameter(os.path.join(logs_dir,'chipqc.Rdata')),
					Redirect(stream=Redirect.STDOUT_APPEND,dest=os.path.join(logs_dir,'chipqc_std.log')),
					Redirect(stream=Redirect.STDERR_APPEND,dest=os.path.join(logs_dir,'chipqc_stderr.log'))
				)



		if '5' in step:#call peaks
			if 'processed_bedpe_to_beds' not in vars():
				processed_bedpe_to_beds=[]
				for lib in libs:
					processed_bedpe_to_beds.append(os.path.join(bed_dir,'{}.bed'.format(lib)))
			for i,bed in enumerate(processed_bedpe_to_beds):
				macs2_callpeak.run(
					Parameter('-t', bed),
					Parameter('-f', 'BEDPE'),
					Parameter('-g', effective_genome_size),
					Parameter('-n', libs[i] + '_regular'),
					Parameter('--outdir',peaks_dir),
					Parameter('--nomodel'),
					Parameter('--extsize', '200'),
					Parameter('--shift', '-100'),
					Parameter('-B', '--SPMR'), # Generates pileup tracks, bedgraph, fragment pileup per million reads
					Parameter('--call-summits'),
					Parameter('--keep-dup', 'all'),
					Redirect(stream=Redirect.STDERR_APPEND,dest=os.path.join(logs_dir,'macs2_narrow.log'))
				)


				if pipeline_args['broad']:
					macs2_callpeak.run(
						Parameter('-t', bed),
						Parameter('-f', 'BEDPE'),
						Parameter('-g', effective_genome_size),
						Parameter('-n', libs[i]+ '_broad'),
						Parameter('--outdir',peaks_dir),
						Parameter('-q', '0.05'),
						Parameter('--nomodel'),
						Parameter('--extsize', '200'),
						Parameter('--shift', '-100'),
						Parameter('--broad'),
						Parameter('--keep-dup', 'all'),
						Redirect(stream=Redirect.STDERR_APPEND,dest=os.path.join(logs_dir,'macs2_broad.log'))
					)
		if '6' in step: #phantompeakqualtools for cross-corelation
			if 'outbams' not in vars():
				outbams=[]
				for i,lib in enumerate(libs):
					outbams.append(os.path.join(align_clean_tmp_dir, '{}.sorted_for_PE.bam'.format(libs[i])))
			for i,bam in enumerate(outbams):
				run_spp_nodups.run(
					Parameter('-p=%s'%threads),
					Parameter('-c=%s'%bam),
					Parameter('-rf'),
					Parameter('-savp=%s'%os.path.join(qc_dir,libs[i]+'_cc.pdf')),
					Parameter('-odir=%s'%qc_dir),
					Parameter('-out=%s'%os.path.join(qc_dir,libs[i]+'.cc.quality')),
					Redirect(stream=Redirect.STDOUT_APPEND,dest=os.path.join(logs_dir,libs[i]+'.quality.log'))
				)
		if '7' in step:#visualization
			for i,lib in enumerate(libs):
				if visual_method=='llr':
					treat=os.path.join(peaks_dir,libs[i] + '_regular_treat_pileup.bdg')
					ctl=os.path.join(peaks_dir,libs[i] + '_regular_control_lambda.bdg')
					macs2_bdgcmp.run(
						Parameter('-t',treat),
						Parameter('-c',ctl),
						Parameter('-m','FE'),
						Parameter('-o',os.path.join(bdgs_dir,libs[i]+'_FE.bdg'))
					)
					macs2_bdgcmp.run(
						Parameter('-t',treat),
						Parameter('-c',ctl),
						Parameter('-m','logLR'),
						Parameter('-p','0.00001'),
						Parameter('-o',os.path.join(bdgs_dir,libs[i]+'_logLR.bdg'))
					)
					#Remove lines from bed file that refer to off-chromosome locations.
					bed_Clip.run(
						Parameter(os.path.join(bdgs_dir,libs[i]+'_FE.bdg')),
						Parameter(pipeline_config['bedtools']['genome-sizes']),
						Parameter(os.path.join(bdgs_dir,libs[i]+'_FE_clip.bdg'))
					)
					bed_Clip.run(
						Parameter(os.path.join(bdgs_dir,libs[i]+'_logLR.bdg')),
						Parameter(pipeline_config['bedtools']['genome-sizes']),
						Parameter(os.path.join(bdgs_dir,libs[i]+'_logLR_clip.bdg'))
					)
					os.system("sort -k 1,1 -k 2,2n %s >%s"%(os.path.join(bdgs_dir,libs[i]+'_FE_clip.bdg'),os.path.join(bdgs_dir,libs[i]+'_FE_clip_sort.bdg')))
					os.system("sort -k 1,1 -k 2,2n %s >%s"%(os.path.join(bdgs_dir,libs[i]+'_logLR_clip.bdg'),os.path.join(bdgs_dir,libs[i]+'_logLR_clip_sort.bdg')))
					bdg2bw.run(
						Parameter(os.path.join(bdgs_dir,libs[i]+'_FE_clip_sort.bdg')),
						Parameter(pipeline_config['bedtools']['genome-sizes']),
						Parameter(os.path.join(bws_dir,libs[i]+'_FE.bw'))
					)
					bdg2bw.run(
						Parameter(os.path.join(bdgs_dir,libs[i]+'_logLR_clip_sort.bdg')),
						Parameter(pipeline_config['bedtools']['genome-sizes']),
						Parameter(os.path.join(bws_dir,libs[i]+'_logLR.bw'))
					)
					staging_delete.append(os.path.join(bdgs_dir,libs[i]+'_FE.bdg'))
					staging_delete.append(os.path.join(bdgs_dir,libs[i]+'_logLR.bdg'))
					staging_delete.append(os.path.join(bdgs_dir,libs[i]+'_FE_clip.bdg'))
					staging_delete.append(os.path.join(bdgs_dir,libs[i]+'_logLR_clip.bdg'))
					staging_delete.append(os.path.join(bdgs_dir,libs[i]+'_FE_clip_sort.bdg'))
					staging_delete.append(os.path.join(bdgs_dir,libs[i]+'_logLR_clip_sort.bdg'))
				elif visual_method=='pileup':
					inbed=os.path.join(bed_dir,'{}.bed'.format(lib))
					outbed=os.path.join(bed_dir,'{}_extend_'.format(lib)+extend+'.bed')
					staging_delete.append(outbed)
					log_filepath=os.path.join(logs_dir,'extend_'+lib+'_'+extend+'.log')
					extend_bed(inbed,outbed,pipeline_config['bedtools']['genome-sizes'],extend,log_filepath)
					extend_sort_bed=os.path.join(bed_dir,lib+'_extend_'+extend+'_sort.bed')
					staging_delete.append(extend_sort_bed)
					os.system('sort -k 1,1 %s >%s'%(outbed,extend_sort_bed))
					outbedgraph=os.path.join(bdgs_dir,lib+'_'+extend+'.bdg')
					staging_delete.append(outbedgraph)
					bedtools_genomecov.run(
						Parameter('-bg'),
						Parameter('-split'),
						Parameter('-i',extend_sort_bed),
						Parameter('-g',pipeline_config['bedtools']['genome-sizes']),
						Redirect(stream=Redirect.STDOUT,dest=outbedgraph)
					)
					#perbase bedgraph
					PerbasebedGraph(outbedgraph)
					#nomalization
					normed_bdg=os.path.join(bdgs_dir,lib+'_'+extend+'.norm.bdg')
					staging_delete.append(normed_bdg)
					Normbedgraph(outbedgraph,normed_bdg)
					normed_sort_bdg=os.path.join(bdgs_dir,lib+'_'+extend+'.norm.sort.bdg')
					staging_delete.append(normed_sort_bdg)
					os.system("sort -k 1,1 -k 2,2n %s >%s"%(normed_bdg,normed_sort_bdg))
					bdg2bw.run(
						Parameter(os.path.join(bdgs_dir,normed_sort_bdg)),
						Parameter(pipeline_config['bedtools']['genome-sizes']),
						Parameter(os.path.join(bws_dir,lib+'.bw'))
					)
				else:
					to_log('wrong visualization method,choose from llr or pileup',atac_logfile)

			
			
		if '10' in step:#fragment length distribution
			for i,bam in enumerate(libs):
				Fragdistribution(os.path.join(align_clean_dir,bam+'.bam'),os.path.join(qc_dir,libs[i]+'_fragLen.pdf'),threads)	

		with open(os.path.join(output_dir,'rm_tmp'),'w+')as f:
			f.write('rm '+' '.join(staging_delete))
		if pipeline_args['rm_tmp']=='y':
			subprocess.call(['rm',staging_delete])