output_dir = config['output_dir']
bam_file = config['bam_file']
bai_file = config['bai_file']

# extract hla reads with samtools view and sort
rule get_hla_reads:
    input:
        bam=bam_file, # tumor bam
        bam_index=bai_file,
    output:
        hla_bam='{outdir}/{subdir}/get_hla_reads/hla.bam'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
        hla_bai='{outdir}/{subdir}/get_hla_reads/hla.bam.bai'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
    log:
        '{outdir}/{subdir}/get_hla_reads/sample.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log'],
        ),
    params:
        gatk_hla_bai='{outdir}/{subdir}/get_hla_reads/hla.bai'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
    shell:
        '{gatk} PrintReads '.format(gatk=config['gatk']['run_file']) +
        '-I {input.bam} ' +
        '-L {interval} '.format(interval=config['gatk']['interval']) + # chrO:0-9 format
        '-O {output.hla_bam} > {log} 2>&1 && ' +
        'mv {params.gatk_hla_bai} {output.hla_bai}' # fix bai name

# optitype samtools conversion step
rule optitype_convert:
    input:
        hla_bam='{outdir}/{subdir}/get_hla_reads/hla.bam'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
        hla_bai='{outdir}/{subdir}/get_hla_reads/hla.bam.bai'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
    output:
        chr6=temp(
            '{outdir}/{subdir}/optitype_convert_all/sample_chr6.fastq'.format(
                outdir=output_dir,
                subdir=config['outputs']['out']
            ),
        ),
        unmapped=temp(
            '{outdir}/{subdir}/optitype_convert_all/sample_unmapped.fastq'.format(
                outdir=output_dir,
                subdir=config['outputs']['out']
            ),
        ),
        chr6_unmapped=temp(
            '{outdir}/{subdir}/optitype_convert_all/sample_chr6_unmapped.fastq'.format(
                outdir=output_dir,
                subdir=config['outputs']['out']
                ),
        ),
        fastq_R1='{outdir}/{subdir}/optitype_convert_all/sample_R1.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        fastq_R2='{outdir}/{subdir}/optitype_convert_all/sample_R2.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    log:
        '{outdir}/{subdir}/optitype_convert_all/sample.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log']
        ),
    benchmark:
        '{outdir}/{subdir}/optitype_convert_all/sample.txt'.format(
            outdir=output_dir,
            subdir=config['outputs']['bench']
        ),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        # Subset reads mapped to chr 6
        'samtools view -b -f3 {input.hla_bam} \"6\" > {output.chr6} 2> {log}; ' +
        # Subset unmapped reads
        'samtools view -b -f13 {input.hla_bam} > {output.unmapped} 2>> {log}; ' +
        # Merge reads mapped to chr 6 and unmapped reads
        '(samtools merge -f {output.chr6_unmapped} {output.chr6} {output.unmapped}) 2>> {log}; ' +
        # Splitting BAM file into paired-ends and create FASTQs
        'samtools view -bf 0x40 {output.chr6_unmapped} | samtools bam2fq - > {output.fastq_R1} 2>> {log}; ' +
        'samtools view -bf 0x80 {output.chr6_unmapped} | samtools bam2fq - > {output.fastq_R2} 2>> {log}'

# optitype razers
rule optitype_extract:
    input:
        fastq_R1='{outdir}/{subdir}/optitype_convert_all/sample_R1.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        fastq_R2='{outdir}/{subdir}/optitype_convert_all/sample_R2.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    output:
        hla_bam_R1='{outdir}/{subdir}/optitype_razers/sample_R1.bam'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        hla_bam_R2='{outdir}/{subdir}/optitype_razers/sample_R2.bam'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    params:
        name='optitype-razers',
        perc_indent=config['hla']['optitype']['perc_indent'],
        max_hits=config['hla']['optitype']['max_hits'],
        dist_range=config['hla']['optitype']['dist_range'],
        hla_ref=config['hla']['optitype']['hla_ref'],
    log:
        '{outdir}/{subdir}/optitype_razers/sample.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log']
        ),
    benchmark:
        '{outdir}/{subdir}/optitype_razers/sample.txt'.format(
            outdir=output_dir,
            subdir=config['outputs']['bench']
        ),
    singularity:
        "docker://fred2/optitype"
    shell:
        # create the command line using the parameters contained in config['razers']
        # add the output, reference and input paths that Snakemake will automatically replace
        # dump all standard and error output to the log file. Redirection is made outside of the container
        '/usr/local/bin/razers3 -i {params.perc_indent} -m {params.max_hits} -dr {params.dist_range} -tc {threads} '
        '-o {output.hla_bam_R1} {params.hla_ref} {input.fastq_R1} &> {log}; '
        '/usr/local/bin/razers3 -i {params.perc_indent} -m {params.max_hits} -dr {params.dist_range} -tc {threads} '
        '-o {output.hla_bam_R2} {params.hla_ref} {input.fastq_R2} &>> {log}'

# optitype samtools conversion step
rule optitype_convert_extracted:
    input:
        hla_bam_R1='{outdir}/{subdir}/optitype_razers/sample_R1.bam'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        hla_bam_R2='{outdir}/{subdir}/optitype_razers/sample_R2.bam'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    output:
        hla_fastq_R1='{outdir}/{subdir}/optitype_convert_extracted/sample_R1.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        hla_fastq_R2='{outdir}/{subdir}/optitype_convert_extracted/sample_R2.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    log:
        '{outdir}/{subdir}/optitype_convert_extracted/sample.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log']
        ),
    benchmark:
        '{outdir}/{subdir}/optitype_convert_extracted/sample.txt'.format(
            outdir=output_dir,
            subdir=config['outputs']['bench']
        ),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        # create the command line using the input
        'samtools bam2fq {input.hla_bam_R1} > {output.hla_fastq_R1} 2> {log}; '
        'samtools bam2fq {input.hla_bam_R2} > {output.hla_fastq_R2} 2>> {log}'

# optitype step
rule optitype:
    input:
        hla_fastq_R1='{outdir}/{subdir}/optitype_convert_extracted/sample_R1.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        hla_fastq_R2='{outdir}/{subdir}/optitype_convert_extracted/sample_R2.fastq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    output:
        coverage='{outdir}/{subdir}/optitype/sample/sample_coverage_plot.pdf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        hla='{outdir}/{subdir}/optitype/sample/sample_result.tsv'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    params:
        name='optitype-normal',
        outdir='{outdir}/{subdir}/optitype/sample'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    log:
        '{outdir}/{subdir}/optitype/sample.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log']
        ),
    benchmark:
        '{outdir}/{subdir}/optitype/sample.txt'.format(
            outdir=output_dir,
            subdir=config['outputs']['bench']
        ),
    singularity:
        "docker://fred2/optitype"
    shell:
        'OPTITYPE_DIR=/usr/local/bin/OptiType; ' +
        # create the command line using the input paths
        'python2.7 $OPTITYPE_DIR/OptiTypePipeline.py -i {input.hla_fastq_R1} {input.hla_fastq_R2} ' +
        # add the output path and some parameters
        '--dna --config $OPTITYPE_DIR/config.ini --prefix sample -v -o {params.outdir} ' +
        # dump all error output to the log file. Redirection is made outside of the container
        '&> {log}'

rule maf_to_vcf:
    input:
        maf=config['maf_file'],
        ref=config['reference_fasta'],
    output:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    log:
        '{outdir}/{subdir}/pvacseq_input/maf2vcf.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log']
        ),
    params:
        outdir='{outdir}/{subdir}/pvacseq_input'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    singularity:
        '/juno/work/shah/mondrian/singularity/variant_v0.0.26.sif' # contains maf2vcf.pl
    shell:
        'maf2vcf.pl --input-maf {input.maf} --ref-fasta {input.ref} ' +
        '--output-dir {params.outdir} ' +
        '--output-vcf {output} &> {log}' 

rule vep:
    input:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
        ref=config['vep']['fasta'],
    output:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    log:
        '{outdir}/{subdir}/pvacseq_input/vep.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log']
        ),
    shell:
        'module load perl && ' +
        '{vep} '.format(vep=config['vep']['run_file']) +
        '--input_file {input.vcf} --format vcf ' +
        '--output_file {output.vcf} --vcf ' +
        '--symbol --terms SO --tsl --hgvs ' +
        '--fasta {input.ref} ' +
        '--offline --cache ' +
        '--plugin Frameshift --plugin Wildtype ' +
        '--dir_plugins {plugins_dir} '.format(plugins_dir=config['vep']['plugins_dir']) + 
        '--dir_cache {cache_dir} '.format(cache_dir=config['vep']['cache_dir']) +
        '--force_overwrite &> {log}'

rule vt_decompose:
    input:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    output:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.decomposed.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    log:
        vcf='{outdir}/{subdir}/pvacseq_input/vt_decompose.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log']
        ),
    shell:
        '{vt} decompose '.format(vt=config['vt']['run_file']) +
        '-s {input.vcf} ' +
        '-o {output.vcf} &> {log}' 

rule bam_readcount:
    input:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.decomposed.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
        bam=config['bam_file'],
        ref=config['reference_fasta'],
    output:
        snv_tsv='{outdir}/{subdir}/pvacseq_input/{sample_id}_bam_readcount_snv.tsv'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
            sample_id=config['sample_id'],
        ),
        indel_tsv='{outdir}/{subdir}/pvacseq_input/{sample_id}_bam_readcount_indel.tsv'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
            sample_id=config['sample_id'],
        ),
    log:
        '{outdir}/{subdir}/pvacseq_input/bam_readcount.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log'],
        ),
    params:
        outdir='{outdir}/{subdir}/pvacseq_input'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        )
    singularity:
        '/juno/work/shah/users/chois7/apollo/neoantigen/bam-readcount_helper_v0.0.3.sif'
    shell:
        'python /usr/bin/bam_readcount_helper.py ' +
        '{input.vcf} ' +
        '{sample_id} '.format(sample_id=config['sample_id']) +
        '{input.ref} ' +
        '{input.bam} ' +
        '{params.outdir} &> {log}'
    
rule vcf_readcount_annotator_snv:
    input:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.decomposed.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
        snv_tsv='{outdir}/{subdir}/pvacseq_input/{sample_id}_bam_readcount_snv.tsv'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
            sample_id=config['sample_id'],
        ),
    output:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.decomposed.snv_annotated.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
    log:
        '{outdir}/{subdir}/pvacseq_input/vcf_readcount_annotator_snv.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log'],
        ),
    shell:
        '{vra} '.format(vra=config['vcf_readcount_annotator']['run_file']) +
        '{input.vcf} {input.snv_tsv} ' +
        'DNA -s {sample_id} '.format(sample_id=config['sample_id']) +
        '-t snv -o {output.vcf} &> {log}'

rule vcf_readcount_annotator_indel:
    input:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.decomposed.snv_annotated.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
        indel_tsv='{outdir}/{subdir}/pvacseq_input/{sample_id}_bam_readcount_indel.tsv'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
            sample_id=config['sample_id'],
        ),
    output:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.decomposed.both_annotated.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
    log:
        '{outdir}/{subdir}/pvacseq_input/vcf_readcount_annotator_indel.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log'],
        ),
    shell:
        '{vra} '.format(vra=config['vcf_readcount_annotator']['run_file']) +
        '{input.vcf} {input.indel_tsv} ' +
        'DNA -s {sample_id} '.format(sample_id=config['sample_id']) +
        '-t indel -o {output.vcf} > {log} 2>&1'

rule pvacseq:
    input:
        vcf='{outdir}/{subdir}/pvacseq_input/consensus.vep.decomposed.both_annotated.vcf'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
        hla='{outdir}/{subdir}/optitype/sample/sample_result.tsv'.format(
            outdir=output_dir,
            subdir=config['outputs']['out']
        ),
    output:
        tsv ='{outdir}/{subdir}/pvacseq/MHC_Class_I/{sample_id}.filtered.tsv'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
            sample_id=config['sample_id'],
        ),
    params:
        outdir='{outdir}/{subdir}/pvacseq'.format(
            outdir=output_dir,
            subdir=config['outputs']['out'],
        ),
        # TODO: reduce modes; e.g. use NetMHCcons for NetMHC, (NetMHCpan), and PickPocket
        # TODO: use threading option -t N_THREADS
        modes = 'MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign',
    singularity: '/juno/work/shah/users/chois7/apollo/pvactools_latest.sif'
    log:
        '{outdir}/{subdir}/pvacseq/sample.log'.format(
            outdir=output_dir,
            subdir=config['outputs']['log'],
        ),
    resources:
        mem_mb=12000,
        disk_mb=4000,
    shell:
        ## input.hla content:
        #         A1      A2      B1      B2      C1      C2      Reads   Objective
        # 0       A*24:02 A*32:01 B*14:01 B*15:01 C*03:03 C*08:02 1300.0  1241.4799999999998
        'hlas=$(cat {input.hla} | tail -n 1 | cut -f2,3,4,5,6,7) && '
        'hlacs="" && '
        'for hla in ${{hlas}}; do echo $hla; hlacs+="HLA-${{hla}},"; done && '
        'pvacseq run ' +
        '--iedb-install-directory /opt/iedb ' + # local IEDB directory
        '--blastp-path /opt/ncbi-blast-2.12.0+/bin/blastp ' + # local BLASTP binary #'--keep-tmp-files ' + # keep temp
        '-t 8 ' + # threading
        '{input.vcf} ' +
        '{sample_id} '.format(sample_id=config['sample_id']) +
        '${{hlacs::-1}} ' + # HLA-A*24:02,...,HLA-C*08:02, -> rm comma at the end
        '{params.modes} ' + 
        '{params.outdir} &> {log}'

