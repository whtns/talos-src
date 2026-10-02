#!/usr/bin/env python3

"""
Takes the annotated VCF, samples and all, and reads it as a MatrixTable.
This rearranges all the annotations into the format expected downstream.
"""

import json
from argparse import ArgumentParser

from loguru import logger

import hail as hl

from cpg_utils.hail_batch import init_batch

from talos.models import DownloadedPanelApp
from talos.utils import get_symbol_to_ensg_mapping, read_json_from_path

MISSING_FLOAT = hl.float64(0)
MISSING_INT = hl.int32(0)
MISSING_STRING = hl.str('')


def extract_and_split_csq_string(vcf_path: str) -> list[str]:
    """
    Extract the BCSQ header from the VCF and split it into a list of strings

    Args:
        vcf_path (str): path to the local VCF

    Returns:
        list of strings
    """

    # get the headers from the VCF
    all_headers = hl.get_vcf_metadata(vcf_path)

    # get the '|'-delimited String of all header names
    csq_whole_string = all_headers['info']['BCSQ']['Description'].split('Format: ')[-1]

    # split it all on pipes, return the list
    return csq_whole_string.lower().split('|')


def csq_strings_into_hail_structs(csq_strings: list[str], mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Take the list of BCSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        mt (hl.MatrixTable): the Mt to annotate

    Returns:
        original MatrixTable with the BCSQ annotations re-arranged
    """

    # get the BCSQ contents as a list of lists of strings, per variant
    split_csqs = mt.info.BCSQ.map(lambda csq_entry: csq_entry.split('\|'))  # noqa: W605

    # this looks pretty hideous, bear with me
    # if BCFtools csq doesn't have a consequence annotation, it will truncate the pipe-delimited string
    # this is fine sometimes, but not when we're building a schema here
    # when we find truncated BCSQ strings, we need to add dummy values to the end of the array
    split_csqs = split_csqs.map(
        lambda x: hl.if_else(
            # if there were only 4 values, add 3 missing Strings
            hl.len(x) == 4,
            x.extend([MISSING_STRING, MISSING_STRING, MISSING_STRING]),
            hl.if_else(
                # 5 values... add 2 missing Strings
                hl.len(x) == 5,
                x.extend([MISSING_STRING, MISSING_STRING]),
                hl.if_else(
                    hl.len(x) == 6,
                    x.extend([MISSING_STRING]),
                    x,
                ),
            ),
        ),
    )

    # transform the CSQ string arrays into structs using the header names
    # Consequence | gene | transcript | biotype | strand | amino_acid_change | dna_change
    mt = mt.annotate_rows(
        transcript_consequences=split_csqs.map(
            lambda x: hl.struct(
                **{csq_strings[n]: x[n] for n in range(len(csq_strings)) if csq_strings[n] != 'strand'},
            ),
        ),
    )

    return mt.annotate_rows(
        # amino_acid_change can be absent, or in the form of "123P" or "123P-124F"
        # we use this number when matching to the codons of missense variants, to find codon of the reference pos.
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                codon=hl.if_else(
                    x.amino_acid_change == MISSING_STRING,
                    hl.missing(hl.tint32),
                    hl.if_else(
                        x.amino_acid_change.matches('^([0-9]+).*$'),
                        hl.int32(x.amino_acid_change.replace('^([0-9]+).+', '$1')),
                        hl.missing(hl.tint32),
                    ),
                ),
            ),
            mt.transcript_consequences,
        ),
    )


def get_mane_annotations(mane_path: str) -> hl.DictExpression:
    """
    Parse the MANE file, and get a dict expression of annotations.
    Args:
        mane_path (str): path to the MANE file

    Returns:
        the dict expression of MANE annotations
    """
    # read in the mane table
    with open(mane_path) as handle:
        mane_dict = json.load(handle)

    # convert the dict into a Hail Dict
    return hl.dict(mane_dict)


def annotate_all_transcript_consequences(
    mt: hl.MatrixTable,
    mane: hl.DictExpression,
    ensgs: hl.DictExpression,
) -> hl.MatrixTable:
    """
    In a single loop, update the AM annotations, MANE annotations, and ENSG gene IDs
    Args:
        mt (MatrixTable): the MatrixTable to annotate
        mane (hl.DictExpression): the MANE annotations
        ensgs (hl.DictExpression): the ENSG annotations

    Returns:
        Original MatrixTable, with reformatted + extended annotations
    """
    key_set = mane.key_set()

    return mt.annotate_rows(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                am_class=hl.if_else(
                    x.transcript == mt.info.am_transcript,
                    mt.info.am_class,
                    MISSING_STRING,
                ),
                am_pathogenicity=hl.if_else(
                    x.transcript == mt.info.am_transcript,
                    mt.info.am_score,
                    MISSING_FLOAT,
                ),
                mane_status=hl.if_else(
                    key_set.contains(x.transcript),
                    mane[x.transcript]['mane_status'],
                    MISSING_STRING,
                ),
                ensp=hl.if_else(
                    key_set.contains(x.transcript),
                    mane[x.transcript]['ensp'],
                    MISSING_STRING,
                ),
                mane_id=hl.if_else(
                    key_set.contains(x.transcript),
                    mane[x.transcript]['mane_id'],
                    MISSING_STRING,
                ),
                gene_id=ensgs.get(x.gene, x.gene),
            ),
            mt.transcript_consequences,
        ),
    )


def clear_alphamissense_sentinel(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Replaces echtvar's miss sentinel in `info.am_score` with a genuine missing value.

    Same defect as the one `nest_spliceai_in_struct` handles, in a field nothing guarded.
    echtvar writes `missing_value` into a numeric INFO field for every variant absent from
    its source instead of leaving it missing, so `am_score` carries a large negative number
    -- observed as -2147480064.0, which is the int32 sentinel -2147483648 after a round
    trip through float32, so the test is a range check rather than an equality check
    against any one value.

    The paired string field `am_class` is genuinely missing on those rows, because echtvar
    does leave strings missing. That disagreement is the tell: measured on the 1,392-sample
    run, 609 of 874 reports (70%) carried a sentinel `am_score` while having no `am_class`
    at all.

    `am_score` is display-only -- `categorybooleanalphamissense` is driven by
    `transcript_consequences.am_pathogenicity` (run_hail_filtering.py:387), which is gated
    on `x.transcript == mt.info.am_transcript` and so was already missing on these rows --
    but it reaches the report JSON and the HTML, where -2147480064.0 reads as a real score
    of extreme benignity rather than as no data.

    Applied before `annotate_all_transcript_consequences`, which is the only consumer of
    `info.am_score`, so the sentinel cannot propagate into `am_pathogenicity` even if
    `am_transcript` were populated on such a row.
    """
    return mt.annotate_rows(
        info=mt.info.annotate(
            am_score=hl.or_missing(
                hl.is_defined(mt.info.am_score) & (mt.info.am_score >= 0),
                mt.info.am_score,
            ),
        ),
    )


def nest_gnomad_in_struct(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Tucks all gnomAD annotations into a hl.Struct"""
    return mt.annotate_rows(
        gnomad=hl.struct(
            gnomad_AC=mt.info.gnomad_AC_joint,
            gnomad_AF=mt.info.gnomad_AF_joint,
            gnomad_AC_XY=mt.info.gnomad_AC_joint_XY,
            gnomad_HomAlt=mt.info.gnomad_HomAlt_joint,
        ),
    )


def nest_spliceai_in_struct(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Tucks the SpliceAI annotations into the struct RunHailFiltering expects.

    `annotate_category_spliceai` (run_hail_filtering.py) reads `mt.splice_ai.delta_score`
    and `mt.splice_ai.splice_consequence`, and returns early when the `splice_ai` row field
    is absent -- which it always was here, because upstream removed the annotation step
    that used to create it. This rebuilds the field from the INFO keys written by the
    SpliceAI echtvar source (scripts/make_spliceai_source.py).

    A row with no SpliceAI entry is normalised to 0.0 / '' rather than left carrying
    echtvar's miss sentinel. echtvar writes `missing_value` (-2147483648) into the INFO
    field for every variant absent from the source rather than leaving it missing --
    measured on a 207,482-row slice of chr1, where exactly the 17,242 rows present in the
    source carried a real score and every other row carried the sentinel. The threshold
    comparison would be safe either way, since -2.1e9 never clears 0.5, but
    `splice_ai_delta` is surfaced in the report JSON and a reported variant with no
    SpliceAI entry must not display -2147483648.
    """
    return mt.annotate_rows(
        splice_ai=hl.struct(
            delta_score=hl.if_else(
                hl.is_defined(mt.info.spliceai_ds) & (mt.info.spliceai_ds >= 0),
                hl.float64(mt.info.spliceai_ds),
                hl.float64(0.0),
            ),
            splice_consequence=hl.if_else(
                hl.is_defined(mt.info.spliceai_csq) & (mt.info.spliceai_csq != '.'),
                hl.str(mt.info.spliceai_csq),
                hl.str(''),
            ),
        ),
    )


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a BCSQ annotated VCF and makes it a HT')
    parser.add_argument('--input', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".ht" extension', required=True)
    parser.add_argument('--panelapp', help='PanelApp download')
    parser.add_argument('--mane', help='Hail Table containing MANE annotations', default=None)
    parser.add_argument(
        '--checkpoint',
        help='Whether to use a remote checkpoint. This is an implicit trigger for the batch backend',
        default=None,
    )
    args = parser.parse_args()

    main(
        vcf_path=args.input,
        output_path=args.output,
        panelapp_path=args.panelapp,
        mane=args.mane,
        checkpoint=args.checkpoint,
    )


def main(
    vcf_path: str,
    output_path: str,
    panelapp_path: str,
    mane: str,
    checkpoint: str | None = None,
):
    """
    Takes a BCFtools-annotated VCF, reorganises into a Talos-compatible MatrixTable
    Will annotate at runtime with AlphaMissense annotations

    Args:
        vcf_path (str): path to the annotated sites-only VCF
        output_path (str): path to write the resulting Hail Table to, must
        panelapp_path (str): PanelApp contents
        mane (str): path to a MANE JSON file for enhanced annotation
        checkpoint (str): which hail backend to use. Defaults to
    """

    if checkpoint:
        logger.info(f'Using Hail Batch backend, checkpointing to {checkpoint}')
        init_batch(
            driver_memory='highmem',
            driver_cores=2,
            worker_memory='highmem',
            worker_cores=2,
        )
    else:
        logger.info('Using Hail Local backend, will use a local checkpoint.')
        hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # pull and split the CSQ header line
    csq_fields = extract_and_split_csq_string(vcf_path=vcf_path)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(vcf_path, array_elements_required=False, force_bgz=True, block_size=20)

    # checkpoint locally to make everything downstream faster
    mt = mt.checkpoint(checkpoint or 'checkpoint.ht', overwrite=True)

    # re-shuffle the BCSQ elements
    mt = csq_strings_into_hail_structs(csq_fields, mt)

    # read the PanelApp data from JSON
    panelapp = read_json_from_path(panelapp_path, return_model=DownloadedPanelApp)

    # Convert JSON data sources into a hl.Dict object
    ensg_dict = get_symbol_to_ensg_mapping(panelapp, as_hail=True)
    mane_dict = get_mane_annotations(mane_path=mane)

    # AlphaMissense, only when the echtvar source was applied. Clear echtvar's miss
    # sentinel before anything reads am_score, so neither the report JSON nor
    # am_pathogenicity can carry -2147480064.0 as though it were a score.
    if 'am_score' in mt.info:
        mt = clear_alphamissense_sentinel(mt)

    # in a single loop, update alphamissense annotations, ENSG gene IDs, and MANE status/matched transcripts
    mt = annotate_all_transcript_consequences(mt, mane_dict, ensg_dict)

    # get a hold of the geneIds - use some aggregation
    mt = mt.annotate_rows(gene_ids=hl.set(mt.transcript_consequences.map(lambda c: c.gene_id)))

    # take note of all named gnomad_* fields
    individual_gnomad_fields = [f for f in mt.info if f.startswith('gnomad_')]

    # gather gnomAD annotations into a separate struct
    mt = nest_gnomad_in_struct(mt)

    # SpliceAI, only when the echtvar source was applied - the fields are absent otherwise
    # and every earlier MatrixTable in this project was built without them
    spliceai_fields = [f for f in mt.info if f.startswith('spliceai_')]
    if spliceai_fields:
        mt = nest_spliceai_in_struct(mt)

    # drop the BCSQ field, all individual gnomAD annotations, and the flat SpliceAI keys
    mt = mt.annotate_rows(info=mt.info.drop('BCSQ', *individual_gnomad_fields, *spliceai_fields))

    mt.describe()

    mt.write(output_path, overwrite=True)


if __name__ == '__main__':
    cli_main()
