nextflow.enable.dsl = 2

reference_file = file(params.reference_file, type: "file")
if (!reference_file.exists()) {
    println("--reference_file: File doesn't exist, check path.")
    exit 1
}

reference_file_index = file("${reference_file}.fai")
if (!reference_file_index.exists()) {
    println("Missing Fasta Index")
    exit 1
}



alignment_list = file(params.alignment_list, type:"file")
if (!alignment_list.exists()){
   println ("--alignment_list: File doesn't exist, check alignment file list.")
   exit 1
}

data=  Channel
       .fromPath(alignment_list)
       .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
       .map { row -> tuple(row.sample_id, row.alignment_file) }

sv_bed = file(params.sv_bed, type: "file", checkIfExists: true)
snv_vcf = file(params.snv_vcf, type: "file", checkIfExists: true)
pedigree_file = file(params.pedigree_file, type: "file", checkIfExists: true)
clf_folder = file(params.clf_folder, type: "file", checkIfExists: true)
exclude_regions_bed = file(params.exclude_regions_bed, type: "file", checkIfExists: true)

// Run the SV2 genotyper:
process SV2 {
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("sv2_output/$sample_id*.vcf.gz")
    script:
    """
    python /work/emanuela.iovino/sv2_2/run_sv2.py --alignment_file $alignment_file --reference_fasta $reference_file --snv_vcf_file $snv_vcf --sv_bed_file $sv_bed --sample_name $sample_id --ped_file $pedigree_file --exclude_regions_bed $exclude_regions_bed --output_folder sv2_output --clf_folder $clf_folder -M
    """
    stub:
    """
    touch $sample_id*.vcf.gz
    """
}

workflow {
    SV2(data)
}
