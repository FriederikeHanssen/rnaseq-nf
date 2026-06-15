include { INDEX } from './index'
include { QUANT } from './quant'
include { FASTQC } from './fastqc'
include { FASTP } from './nf-core/fastp'

workflow RNASEQ {
    take:
    read_pairs_ch
    transcriptome

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(read_pairs_ch)

    // Adapt channel shape for nf-core fastp module:
    // Pipeline format: [id, fastq_1, fastq_2]
    // nf-core format:  [meta_map, [fastq_1, fastq_2], adapter_fasta]
    fastp_input_ch = read_pairs_ch.map { id, fastq_1, fastq_2 ->
        [[id: id, single_end: false], [fastq_1, fastq_2], []]
    }

    FASTP(fastp_input_ch, false, false, false)

    // Convert fastp output back to pipeline format: [id, fastq_1, fastq_2]
    trimmed_reads_ch = FASTP.out.reads.map { meta, reads ->
        [meta.id, reads[0], reads[1]]
    }

    quant_ch = QUANT(trimmed_reads_ch, index)
    samples_ch = fastqc_ch.join(quant_ch)

    emit:
    samples = samples_ch
}
