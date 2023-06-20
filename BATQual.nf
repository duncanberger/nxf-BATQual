#!/usr/bin/env nextflow

// Write help message
def helpMessage() {
	log.info"""
Usage:
	nextflow run BATQual.nf --input read_locations.csv --mode fastq
	nextflow run BATQual.nf --input fasta_locations.csv --mode fasta

Mandatory arguments:
	--input			Path to csv file with read locations
	--mode			Running mode either: fastq or fasta depending on input

Velvet options:
	--min_k			Minimum k-mer length for assembly [90]
	--max_k			Maximum k-mer length for assembly [141]
	--step			Number of hash values to skip each increment from min_k to max_k [2] 

BUSCO options:
	--lineage		Specify the BUSCO lineage to be used ["lactobacillales_odb10"]

Run options:
	--outdir		The output directory where the results will be saved ["results"]

Skip metrics:
	--run_GPSC	Run GPSC estimation [false]
	--pneumo	Run tools specific to Streptococcus pnuemoniae [false]
	--no_CheckM	Don't run CheckM [false]
	--no_Kraken2	Don't run Kraken2 [false]
	--no_HET_call	Don't call heterozygous variants [false]
	--no_Prokka	Don't annotate genomes [false]
    """.stripIndent()
}

// Set default params for mandatory inputs
params.help = null
params.input = null
params.mode = null

// Set up error message for mandatory params
if (params.help){
    helpMessage()
    exit 0
}
if (!params.input){
    log.error "Error: '--input' parameter is missing"
    helpMessage()
    exit 0
}
if (!params.mode){
    log.error "Error: '--mode' parameter is missing"
    helpMessage()
    exit 0
}

// Set default params for optional inputs
params.no_GPSC = true
params.pneumo = false
params.run_GPSC = false
params.index = file(params.input)
params.no_HET_call = false
params.no_CheckM = false
params.no_Prokka = false
params.no_Kraken2 = false

// FASTQ workflow with conditional execution of processes specific to Streptococcus pneumoniae
workflow MAIN_A {
    Boolean flag = false
// Read file and split input csv into columns of sample_id, forward and reverse reads
    Channel
     .fromPath(params.index)
     .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
     .splitCsv(header:true, by: 1)
     .map { row-> tuple(row.sample_id, file(row.read1), file(row.read2)) }
     .set {reads}
    main:
    FASTP(reads)
    FASTQC(reads)
    VELVET(FASTP.out)
    if (params.no_Kraken2 == false ) {
    KRAKEN(reads)
    }
    else {}
    BUSCO(VELVET.out)
    STATS(VELVET.out)
    if (params.no_Prokka == false ) {
    PROKKA(VELVET.out)
    }
    else {}
    if (params.no_CheckM == false ) {
    CHECKM(VELVET.out)
    }
    else {}
    if (params.no_HET_call == false ) {
    SHET(FASTP.out)
    }
    else {}
    MASH(VELVET.out)
    MLST(VELVET.out)

// Set conditions under which serotyping and GPSC assignment will occur. 
    if (params.pneumo == true & params.run_GPSC == true ) {
        GPSC(VELVET.out)
		PK(VELVET.out)
		SBA(reads)
		}
	else if (params.pneumo == true & params.run_GPSC == false ) {
		PK(VELVET.out)
		SBA(reads)
		}
	else {}
    QUAST(VELVET.out)
    }

// FASTA workflow with conditional execution of processes specific to Streptococcus pneumoniae
workflow MAIN_B {
    Boolean flag = false
// Read file and split input csv into columns of sample_id, FASTA file
    Channel
     .fromPath(params.index)
     .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
     .splitCsv(header:true, by: 1)
     .map { row-> tuple(row.sample_id, file(row.fasta)) }
     .set {fastas}
    main:
    BUSCO(fastas)
    STATS(fastas)
    if (params.no_Prokka == false ) {
    PROKKA(fastas)
    }
    else {}
    if (params.no_CheckM == false ) {
    CHECKM(fastas)
    }
    else {}
    MASH(fastas)
// Set conditions under which serotyping and GPSC assignment will occur. 
	if (params.pneumo == true & params.run_GPSC == true ) {
        GPSC(fastas)
		PK(fastas)
		}
	else if (params.pneumo == true & params.run_GPSC == false ) {
		PK(fastas)
		}
	else {}
    QUAST(fastas)
    MLST(fastas)
}

// Run either workflow depending on input mode
workflow {
    if (params.mode == "fastq") {
        MAIN_A()}
    else if (params.mode == "fasta") {
        MAIN_B() }
    else {
        helpMessage()
        exit 1
    }
}

// Run fastp
process FASTP {
    errorStrategy 'ignore'
    cpus = 2
    tag "Running fastp on $sample_id"
 
    input:
    tuple val(sample_id), path(read1), path(read2)
 
    output:
    tuple val(sample_id), path("${sample_id}.R1.fq.gz"), path("${sample_id}.R2.fq.gz")
 
    script:
    """
    fastp -w 2 -i ${read1} -I ${read2} -o ${sample_id}.R1.fq.gz -O ${sample_id}.R2.fq.gz
    """
}

// Run FASTQC, FASTQC doesn't let you specify output file names for each run, so some slightly inelegant renaming is required
process FASTQC {
    errorStrategy 'ignore'
    cpus = 1
    tag "Running FastQC on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple path("${sample_id}.R1_fastqc.zip"), path("${sample_id}.R2_fastqc.zip")

    script:
    """
    fastqc ${read1}
    mv *.zip ${sample_id}.R1_fastqc.tx1
    fastqc ${read2}
    mv *.zip ${sample_id}.R2_fastqc.zip
    mv *.tx1 ${sample_id}.R1_fastqc.zip
    """
}

// Run Velvet for the specific kmer range - default will be 66-90% of read length
process VELVET {
    errorStrategy = 'ignore'
    maxRetries = 1
    cpus 2
    tag "Running Velvet on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.R1.fq.gz"), path("${sample_id}.R2.fq.gz")
    output:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    script:
    """
    VelvetOptimiser.pl -s "$params.min_k" -e "$params.max_k" --p ${sample_id} -x "$params.step" -t 2 --o '-very_clean yes' -f '-shortPaired -separate -fastq.gz ${sample_id}.R1.fq.gz ${sample_id}.R2.fq.gz'
    mv ${sample_id}_data_*/contigs.fa ${sample_id}.velvet_contigs.fa
    rm ${sample_id}_data_*/*
    """
}

// Run kraken2, writing both genus and species assignments to file
process KRAKEN {
    errorStrategy 'ignore'
    cpus 2
    tag "Running Kraken2 on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.R1.fq.gz"), path("${sample_id}.R2.fq.gz")

    output:
    tuple path("${sample_id}.kraken2_results.txt"), path("${sample_id}.kraken2.out")

    script:
    """
    kraken2 --db $baseDir/DB/ --output txc --use-names --report ${sample_id}.kraken2.out --paired ${sample_id}.R1.fq.gz ${sample_id}.R2.fq.gz
    awk '\$4=="G"' ${sample_id}.kraken2.out | sort -k4,4 -gr | awk '{print "$sample_id","kraken2_genus",\$6":"\$1","}' OFS=","| head -1 > ${sample_id}.kraken2_results.txt
    awk '\$4=="S"' ${sample_id}.kraken2.out | sort -k4,4 -gr | awk '{print "$sample_id","kraken2_species",\$6" "\$7":"\$1","}' OFS=","| head -1 >> ${sample_id}.kraken2_results.txt
    """
}

// Run pneumoKITy
process PK {
    errorStrategy 'ignore'
    tag "Running pneumoKITy on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")
    
    output:
    path "${sample_id}.pneumoKITy_results.txt"

    script:
    """
    python $baseDir/PneumoKITy/pneumokity.py pure -a ${sample_id}.velvet_contigs.fa -t 1 -o ./
    grep -A1 Predicted pneumo_capsular_typing/*_serotyping_results.txt | sed 's/ //g' | cut -f2 | paste - - -d, | awk '{print "${sample_id}","pneumoKITy_serotype",\$1,""}' OFS=',' | sed 's/GREEN/PASS/g' | sed 's/AMBER/WARNING/g' | sed 's/RED/FAIL/g' > ${sample_id}.pneumoKITy_results.txt
    """
}

// Run CheckM
process CHECKM {
    errorStrategy 'ignore'
    tag "Running CHECKM on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    path "${sample_id}.checkM_results.txt"

    script:
    """
    checkm taxonomy_wf -x fa --tab_table -f ${sample_id}_check.tbl domain "Bacteria" ./ ${sample_id}_test
    cat ${sample_id}_check.tbl | cut -f12,13,14 | sed 's/ /_/g' | awk -F'\t' '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS \$i: \$i) } END{ for (i in a) print a[i] }'| awk '{print "${sample_id},CHECKM_"\$1,\$2,""}' OFS="," > ${sample_id}.checkM_results.txt
    """
}

// Calculate variant site heterozygosity using minimap2 and BCFtools. Deleting BAMs at the end to save space. 
process SHET {
    errorStrategy 'ignore'
    tag "Calculating heterozygosity of $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.R1.fq.gz"), path("${sample_id}.R2.fq.gz")

    output:
    path "${sample_id}.hetperc_results.txt"

    script:
    """
    minimap2 -R '@RG\\tSM:DRAK\\tID:DRAK' -ax sr $baseDir/DB/$params.ref_assembly ${sample_id}.R1.fq.gz ${sample_id}.R2.fq.gz | samtools sort -O BAM -o ${sample_id}.bam - 
    freebayes -iXu -b ${sample_id}.bam --min-coverage 5 --limit-coverage 30 -f $baseDir/DB/$params.ref_assembly | bcftools view -i 'QUAL>=20' | bcftools query -f '[%GT]\n' - | awk '{if(\$0=="0/1" || \$0=="1/2"){nmw+=1}}END{print ((nmw*100)/NR)}' | awk '{print "${sample_id}","perc_het_vars",\$1,""}' OFS=',' > ${sample_id}.hetperc_results.txt
    """
}

// Run seroBA
process SBA {
    errorStrategy 'ignore'
    tag "Running seroBA on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.seroBA_results.txt")

    script:
    """
    seroba runSerotyping --databases $baseDir/DB/database/ ${read1} ${read2} ${sample_id}_serotype
    awk '{print "${sample_id}","seroBA_serotype", \$2,""}' OFS=',' ${sample_id}_serotype/pred.tsv > ${sample_id}.seroBA_results.txt
    """
}

// Use a custom python script to calculate assembly statistics
process STATS {
    errorStrategy 'ignore'
    tag "Calculating assembly stats of $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    path "${sample_id}.assembly_results.txt"

    script:
    """
    python $baseDir/scripts/asm_stats.py --fasta "${sample_id}".velvet_contigs.fa --isolate "${sample_id}" --gap 1 --output ${sample_id}.assembly_results
    """
}

// Run BUSCO using the specified database (config file)
process BUSCO {
    errorStrategy 'ignore'
    tag "Running BUSCO on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    path "${sample_id}.BUSCO_results.txt"

    script:
    """
    busco -i ${sample_id}.velvet_contigs.fa -m geno -o ${sample_id}_BUSCO --offline -c 1 -l $baseDir/DB/"$params.lineage"
    grep "C:" ${sample_id}_BUSCO/short_summary.specific.*.txt | cut -f2 | tr ':' '\n' | cut -f1 -d "%" | head -6 | tail -4 | paste <(printf "BUSCO_complete_single_copy\nBUSCO_complete_duplicated\nBUSCO_fragmented\nBUSCO_missing\n") - | awk '{print "${sample_id}",\$1,\$2,""}' OFS="," > ${sample_id}.BUSCO_results.txt
    """
}

// Run Prokka
process PROKKA {
    errorStrategy 'ignore'
    cpus 4
    tag "Running Prokka on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    tuple path("${sample_id}.prokka.gff"), path("${sample_id}.prokka.fna"), path("${sample_id}.prokka_results.txt")

    script:
    """
    prokka --cpus 4 --outdir prokka_${sample_id}/ --locustag B_SPN ${sample_id}.velvet_contigs.fa --compliant --centre BDI --prefix ${sample_id}.prokka --strain ${sample_id}    
    mv prokka_${sample_id}/* ./
    awk '{if (\$1=="gene:") {print "${sample_id}","n_genes",\$2,""} }' OFS=','  ${sample_id}.prokka.txt > ${sample_id}.prokka_results.txt
    """
}

// Run Mash
process MASH {
    errorStrategy 'ignore'
    cpus 1
    tag "Running Mash on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    path "${sample_id}.mash_results.txt"

    script:
    """
    mash sketch ${sample_id}.velvet_contigs.fa
    mash dist $baseDir/DB/refseq.genomes.k21s1000.msh ${sample_id}.velvet_contigs.fa.msh | awk '\$4<0.05' | sort -k3,3 -g | head -5  | cut -f1,2 -d"_" | cut -f1 -d"." | while read file; do grep \$file $baseDir/DB/assembly_summary_genbank.txt; done | cut -f8 | cut -f1,2 -d" " | uniq | paste -sd';' | awk -F '\t' '{print "${sample_id}","MASH_hit",\$1,""}' OFS=',' > ${sample_id}.mash_results.txt
    """
}

// Assign GPS clusters
process GPSC {
    errorStrategy 'ignore'
    cpus 1
    tag "Assigning GPSCs to $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    path "${sample_id}.GPSC_results.txt"

    script:
    """
    echo ${sample_id} ${sample_id}.velvet_contigs.fa | tr ' ' '\t' > GPSC_input.txt    
    poppunk_assign --db $baseDir/DB/GPS_v6 --distances $baseDir/DB/GPS_v6/GPS_v6.dists --query GPSC_input.txt --output ${sample_id}_GPSC --external-clustering $baseDir/DB/GPS_v6_external_clusters.csv --threads 1
    cat ${sample_id}_GPSC/${sample_id}_GPSC_external_clusters.csv | tail -1 | awk -F, '{print \$1,"GPSC",\$2,""}' OFS="," > ${sample_id}.GPSC_results.txt 
    """
}

// Run QUAST
process QUAST {
    errorStrategy 'ignore'
    cpus 1
    tag "Running QUAST on $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    path "${sample_id}.quast"

    script:
    """
    quast -o ${sample_id}.quast ${sample_id}.velvet_contigs.fa
    """
}

// Assign MLST alleles
process MLST {
    errorStrategy 'ignore'
    cpus 1
    tag "Assigned MLST alleles to $sample_id"
    publishDir "$params.outdir/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}.velvet_contigs.fa")

    output:
    path "${sample_id}.mlst_results.txt"

    script:
    """
    mlst ${sample_id}.velvet_contigs.fa > ${sample_id}.temp
    cat ${sample_id}.temp | sed 's/,/;/g' | tr '\t' '\n' | grep "("  | sed 's/^/${sample_id},mlst_allele_/g' | sed 's/\$/,/g' | sed 's/(/,/g' | sed 's/)/,/g' >> ${sample_id}.mlst_results.txt
    cat ${sample_id}.temp | head -2 | tail -1 | awk '{print "$sample_id","mlst_species",\$2,","}' OFS="," | sed 's/)//g' | sed 's/(/,/g' >> ${sample_id}.mlst_results.txt
    """
}
