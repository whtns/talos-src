process AnnotateWithEchtvar {
    container params.container

    input:
        tuple val(cohort), path(vcf), path(tbi)
        path gnomad_zip
        path am_zip
        path spliceai_zip

    output:
        tuple val(cohort), path("${vcf.simpleName}_echtvar.vcf.bgz"), path("${vcf.simpleName}_echtvar.vcf.bgz.tbi")

    script:
        // optional third source; empty when no SpliceAI zip is configured
        def spliceai_arg = spliceai_zip ? "-e ${spliceai_zip}" : ''
        """
        set -euo pipefail

        echtvar anno \
            -e ${gnomad_zip} \
            -e ${am_zip} \
            ${spliceai_arg} \
            -i "gnomad_AF_joint < 0.05" \
            ${vcf} \
            "${vcf.simpleName}_echtvar.vcf.bgz"
        tabix "${vcf.simpleName}_echtvar.vcf.bgz"
        """
}
