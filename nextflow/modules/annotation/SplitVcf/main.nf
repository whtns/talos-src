process SplitVcf {
    container params.container

    input:
        tuple val(cohort), path(vcf)

    output:
        tuple val(cohort), path("split_*.bgz")

    script:
    """
    set -euo pipefail

    # save the header
    bcftools view -h ${vcf} > header.txt

    # Split the body into fragments of about `vcf_split_n` rows each, but never break a
    # locus across two fragments: a new fragment opens only once the quota is reached AND
    # CHROM:POS has changed. See patches/0002 for why the line-count split was not enough.
    bcftools view -H ${vcf} | awk -v n=${params.vcf_split_n} -v hdr=header.txt '
    {
        key = \$1 ":" \$2
        if (started && count >= n && key != prev) { close(cmd); started = 0; idx++ }
        if (!started) {
            cmd = sprintf("bgzip -c > split_x%02d.vcf.bgz", idx)
            while ((getline line < hdr) > 0) print line | cmd
            close(hdr)
            started = 1
            count = 0
        }
        print | cmd
        count++
        prev = key
    }
    END { if (started) close(cmd) }
    '
    """
}
