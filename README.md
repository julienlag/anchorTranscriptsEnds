# NAME

anchorTranscriptEnds

# SYNOPSIS

A utility to prepare a GTF file of transcriptome reads for anchored transcript merging.

**Usage**:

`anchorTranscriptEnds.pl <gff> <5p_supported_reads> <3p_supported_reads> <5p_clusters_bed> <3p_clusters_bed>`

## ARGUMENTS/INPUT

- **arg1**: (&lt;gff>) Path to GTF file containing the aligned reads (or '-' for STDIN).

    &lt;gff> **must** contain exon records (all other features will be skipped), grouped by transcript\_id.

- **arg2**: (<5p\_supported\_reads>) Path to the file containing the transcript\_id's of the reads present in &lt;gff> that are to be 5'-anchored.

    <5p\_supported\_reads> must contain one transcript\_id per line.

- **arg3**: (<3p\_supported\_reads>) Path to the file containing the transcript\_id's of the reads present in &lt;gff> that are to be 3'-anchored.

    <3p\_supported\_reads> must contain one transcript\_id per line.

- **arg4**: (<5p\_clusters\_bed>) Path to the BED6 file containing the coordinates of the TSS clusters, with transcript\_id's present in a comma-separated list in the 4th field (as obtained through _e.g._ `bedtools merge -c 4 -o collapse -d 5 -i <raw_TSSs.bed>`). These will be used to adjust the coordinates of supported TSSs, meaning that all transcript\_id's present in <5p\_supported\_reads> should be present in <5p\_clusters\_bed> (others will be ignored). The TSSs of all transcript\_id's present in <5p\_supported\_reads> will be adjusted according to the corresponding TSS cluster in <5p\_clusters\_bed> (_i.e._, they will be extended or shortened to the cluster's 5' end).
- **arg5**: (<3p\_clusters\_bed>) Path to the BED6 file containing the coordinates of the TTS clusters, with transcript\_id's present in a comma-separated list in the 4th field (as obtained through _e.g._ `bedtools merge -c 4 -o collapse -d 5 -i <raw_TTSs.bed>`). These will be used to adjust the coordinates of supported TTSs, meaning that all transcript\_id's present in <3p\_supported\_reads> should be present in <3p\_clusters\_bed> (others will be ignored). The TTSs of all transcript\_id's present in <3p\_supported\_reads> will be adjusted according to the corresponding TTS cluster in <3p\_clusters\_bed> (_i.e._, they will be extended or shortened to the cluster's 3' end).

## OUTPUT

# DESCRIPTION

# DEPENDENCIES

CPAN: FindBin, Storable qw(dclone)

Custom modules (in current github): gffToHash, hashToGff

# AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com
