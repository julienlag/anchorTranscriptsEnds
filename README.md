# NAME

anchorTranscriptsEnds

# SYNOPSIS

A utility to prepare a GTF file of transcriptome reads for anchored transcript merging.

**Usage**:

`anchorTranscriptsEnds.pl <gff> <5p_supported_reads> <3p_supported_reads> <5p_clusters_bed> <3p_clusters_bed>`

## ARGUMENTS/INPUT

- **arg1** <**gff**>: Path to GTF file containing the aligned reads (or '-' for STDIN).

    &lt;gff> **must** contain exon records (all other features will be skipped), grouped by transcript\_id.

- **arg2** <**5p\_supported\_reads**>: Path to the file containing the transcript\_id's of the reads present in &lt;gff> that are to be 5'-anchored.

    <5p\_supported\_reads> must contain one transcript\_id per line.

- **arg3** <**3p\_supported\_reads**>: Path to the file containing the transcript\_id's of the reads present in &lt;gff> that are to be 3'-anchored.

    <3p\_supported\_reads> must contain one transcript\_id per line.

- **arg4** <**5p\_clusters\_bed**>: Path to the BED6 file containing the coordinates of the TSS clusters, with transcript\_id's present in a comma-separated list in the 4th field (as obtained through _e.g._ `bedtools merge -c 4 -o collapse -s -d 5 -i <raw_TSSs.bed>`).

    These will be used to adjust the coordinates of supported TSSs, meaning that all transcript\_id's present in <5p\_supported\_reads> should be present in <5p\_clusters\_bed> (others will be ignored).

    The TSSs of all transcript\_id's present in <5p\_supported\_reads> will be adjusted according to the corresponding TSS cluster in <5p\_clusters\_bed> (_i.e._, they will be extended to the cluster's 5' end).

- **arg5** <**3p\_clusters\_bed**>: Path to the BED6 file containing the coordinates of the TTS clusters, with transcript\_id's present in a comma-separated list in the 4th field (as obtained through _e.g._ `bedtools merge -c 4 -o collapse -s -d 5 -i <raw_TTSs.bed>`).

    These will be used to adjust the coordinates of supported TTSs, meaning that all transcript\_id's present in <3p\_supported\_reads> should be present in <3p\_clusters\_bed> (others will be ignored).

    The TTSs of all transcript\_id's present in <3p\_supported\_reads> will be adjusted according to the corresponding TTS cluster in <3p\_clusters\_bed> (_i.e._, they will be extended to the cluster's 3' end).

# DESCRIPTION

This program prepares a GTF file containing trancriptome read mappings for "anchored merging" (see Lagarde _et al._, 2017,  https://doi.org/10.1101/105064 for details about this procedure). Its output can be fed to a standard transcript merging program (_e.g._ [Compmerge](https://github.com/sdjebali/Compmerge), [Cuffmerge](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/)) to obtain anchor-merged transcript models.

Any transcript present in the input (&lt;gff>) will be edited according to the following rules:

- If its transcript\_id is **absent from** <**5p\_supported\_reads**>, its **TSS will be left untouched**.
- If its transcript\_id is **absent from** <**3p\_supported\_reads**>, its **TTS will be left untouched**.
- If its transcript\_id is **present in** <**5p\_supported\_reads**>, its 5' end structure will be affected the following way:
    - **First**, its 5' end coordinates will be **adjusted** according to the information contained in <5p\_clusters\_bed>, _i.e._, they will be extended to the 5'end of the TSS cluster the transcript belongs to. In doing so, we ensure that within a cluster, all TSSs align at the exact same position.
    - **Second**, a chain of artificial, "biologically impossible" exons will be added upstream of the adjusted TSS (_i.e._, four 1 nucleotide-long exons, separated by 3 nucleotide-long introns) to the transcript model. These **false exons** (identified with `**fakeExon "yes";**` in the 9th field of the GFF output) serve as **anchors** to supported TSSs during the merging step. In other words, they prevent the transcript they belong to from being merged into a longer transcript "container".
- If its transcript\_id is **present in** <**3p\_supported\_reads**>, its 3' end structure will be affected the following way:
    - **First**, its 3' end coordinates will be **adjusted** according to the information contained in <3p\_clusters\_bed>, _i.e._, they will be extended to the 3'end of the TTS cluster the transcript belongs to. In doing so, we ensure that within a cluster, all TTSs align at the exact same position.
    - **Second**, a chain of artificial, "biologically impossible" exons will be added downstream of the adjusted TTS (_i.e._, four 1 nucleotide-long exons, separated by 3 nucleotide-long introns) to the transcript model. These **false exons** (identified with `**fakeExon "yes";**` in the 9th field of the GFF output) serve as **anchors** to supported TTSs during the merging step. In other words, they prevent the transcript they belong to from being merged into a longer transcript "container".

**IMPORTANT WARNING: Any false exon (identified with `fakeExon "yes";` in anchorTranscriptsEnds' output) generated in the process and fed to the merging program NEEDS TO BE REMOVED AFTER MERGING.**

# OUTPUT

# DEPENDENCIES

CPAN: FindBin, Storable qw(dclone)

Custom modules (in current github): gffToHash, hashToGff

# AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com
