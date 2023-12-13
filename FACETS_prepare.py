import os
import sys

project_dir = os.getcwd()
for line in open(sys.argv[1]):
	line =line.strip().split('\t')
	tumorid = line[0]
	normalid = line[1]
	os.makedirs ('prepare/facets/%s'%tumorid)
    #generate normal control vcf file
	cmd1 = '/bin/java -Djava.io.tmpdir=/jxf/JavaTemp -jar /jxf/softwares/GATK-3.6/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /jxf/databases/hs37d5/bwa_index/hs37d5.fa -I %s/prepare/bam/%s.sorted.rmdup.realigned.recal.bam -glm BOTH --min_indel_count_for_genotyping 3 --min_indel_fraction_per_sample 0.05 -stand_emit_conf 1 -stand_call_conf 1 -L %s/raw.bed -o %s/prepare/facets/%s/%s.ug.vcf -ploidy 2 -nt 8 --disable_auto_index_creation_and_locking_when_reading_rods'%(project_dir,normalid,project_dir,project_dir,tumorid,normalid)
    #obtain snp file containing ref count and alter count
	cmd2 = '/jxf/softwares/facets/facets-master/inst/extcode/snp-pileup -q 10 -Q 15 -P 20 -r 10,0 %s/prepare/facets/%s/%s.ug.vcf %s/prepare/facets/%s/%s.snp %s/prepare/bam/%s.sorted.rmdup.realigned.recal.bam %s/prepare/bam/%s.sorted.rmdup.realigned.recal.bam'%(project_dir,tumorid,normalid,project_dir,tumorid,tumorid,project_dir,normalid,project_dir,tumorid)
    # run facets command
	cmd3 = 'Rscript /jxf/scripts/facets.R %s/prepare/facets/%s/ %s'%(project_dir,tumorid,tumorid)
	cmd = cmd1 + '&&' + cmd2 + '&&' + cmd3
	open('prepare/scripts/'+tumorid+'.sh','w').write(cmd)
	
    


