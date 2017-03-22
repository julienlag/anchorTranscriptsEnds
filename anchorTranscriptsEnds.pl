#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;    # find present script
use lib "$FindBin::Bin";
use gffToHash;
use hashToGff;
use Pod::Usage;
use Storable qw(dclone);


my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to
my $sections = "NAME|SYNOPSIS|DESCRIPTION";


=head1 NAME

anchorTranscriptsEnds

=head1 SYNOPSIS

A utility to prepare a GTF file of transcriptome reads for anchored transcript merging.

B<Usage>:

C<< anchorTranscriptsEnds.pl <gff> <5p_supported_reads> <3p_supported_reads> <5p_clusters_bed> <3p_clusters_bed> >>

=head2 ARGUMENTS/INPUT

=over

=item B<arg1> <B<gff>>: Path to GTF file containing the aligned reads (or '-' for STDIN).

<gff> B<must> contain exon records (all other features will be skipped), grouped by transcript_id.

=item B<arg2> <B<5p_supported_reads>>: Path to the file containing the transcript_id's of the reads present in <gff> that are to be 5'-anchored.

<5p_supported_reads> must contain one transcript_id per line.

=item B<arg3> <B<3p_supported_reads>>: Path to the file containing the transcript_id's of the reads present in <gff> that are to be 3'-anchored.

<3p_supported_reads> must contain one transcript_id per line.

=item B<arg4> <B<5p_clusters_bed>>: Path to the BED6 file containing the coordinates of the TSS clusters, with transcript_id's present in a comma-separated list in the 4th field (as obtained through I<e.g.> C<< bedtools merge -c 4 -o collapse -s -d 5 -i <raw_TSSs.bed> >>).

These will be used to adjust the coordinates of supported TSSs, meaning that all transcript_id's present in <5p_supported_reads> should be present in <5p_clusters_bed> (others will be ignored).

The TSSs of all transcript_id's present in <5p_supported_reads> will be adjusted according to the corresponding TSS cluster in <5p_clusters_bed> (I<i.e.>, they will be extended or shortened to the cluster's 5' end).


=item B<arg5> <B<3p_clusters_bed>>: Path to the BED6 file containing the coordinates of the TTS clusters, with transcript_id's present in a comma-separated list in the 4th field (as obtained through I<e.g.> C<< bedtools merge -c 4 -o collapse -s -d 5 -i <raw_TTSs.bed> >>).

These will be used to adjust the coordinates of supported TTSs, meaning that all transcript_id's present in <3p_supported_reads> should be present in <3p_clusters_bed> (others will be ignored).

The TTSs of all transcript_id's present in <3p_supported_reads> will be adjusted according to the corresponding TTS cluster in <3p_clusters_bed> (I<i.e.>, they will be extended or shortened to the cluster's 3' end).


=back


=head1 DESCRIPTION

This program prepares a GTF file containing trancriptome read mappings for "anchored merging" (see Lagarde I<et al.>, 2017,  https://doi.org/10.1101/105064 for details about this procedure). Its output can be fed to standard transcript merging software (I<e.g.> L<https://github.com/sdjebali/Compmerge|Compmerge>, L<http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/|Cuffmerge>)

at is, we merged close and overlapping sites using the bedtools merge utility, with a maximum clustering distance of 5 bases ("-d 5"), and forcing strandedness ("-s"). Each individual TSS/TTS belonging to a cluster was assigned its start/end coordinate, respectively -- meaning that terminal exons were sometimes extended by a few nucleotides when necessary. In doing so, we ensured that within a cluster, all sites aligned at the exact same position. We subsequently added an "anchor" to all high-confidence, adjusted sites. This step consisted in attaching an artificial, biologically impossible chain of exons (i.e., four 1 nucleotide-long exons, separated by 3 nucleotide-long introns) to each transcript model, upstream or downstream of its high-confidence TSS or TTS, respectively. These fake exons served as anchors to supported start and termination sites during the merging step, and were discarded immediately afterwards.


=head1 OUTPUT

=head1 DEPENDENCIES

CPAN: FindBin, Storable qw(dclone)

Custom modules (in current github): gffToHash, hashToGff

=head1 AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com


=cut

unless ($#ARGV == 4){
	$message_text="ERROR: Invalid number or arguments (must be 4)\n";
	pod2usage( { -message => $message_text,
       		     -exitval => $exit_status  ,
           		-verbose => $verbose_level,
               -output  => $filehandle } );

}

my $transcripts_gff=$ARGV[0];
my $fivep_supported_transcripts_list=$ARGV[1];
my $threep_supported_transcripts_list=$ARGV[2];
my $bedOfClusteredFivepEnds=$ARGV[3]; #will only affect 5' supported transcripts
my $bedOfClusteredThreepEnds=$ARGV[4]; # will only affect 3' supported transcripts

my %fivepsupported=();
my %threepsupported=();
my %tssClustersFivepEnds=();
my %ttsClustersThreepEnds=();

print STDERR "Reading list of 5' supported transcripts from $fivep_supported_transcripts_list ...\n";
open F, "$fivep_supported_transcripts_list" or die $!;
while (<F>){
	chomp;
	$fivepsupported{$_}=1;
}
close F;
print STDERR "Done\n";
print STDERR "Reading list of 3' supported transcripts from $threep_supported_transcripts_list ...\n";

open F, "$threep_supported_transcripts_list" or die $!;
while (<F>){
	chomp;
	$threepsupported{$_}=1;
}
close F;
print STDERR "Done\n";

print STDERR "Reading 5' BED coordinates of TSS clusters $bedOfClusteredFivepEnds ...\n";
open F, "$bedOfClusteredFivepEnds" or die $!;
while(<F>){
	chomp;
	my @line=split "\t";
	my @ids=split(",", $line[3]);
	foreach my $id (@ids){
		push(@{$tssClustersFivepEnds{$id}}, $line[0], $line[2], $line[5])
	}
}
close F;
print STDERR "Done\n";
print STDERR "Reading 3' BED coordinates of TTS clusters $bedOfClusteredThreepEnds ...\n";
open F, "$bedOfClusteredThreepEnds" or die $!;
while(<F>){
	chomp;
	my @line=split "\t";
	my @ids=split(",", $line[3]);
	foreach my $id (@ids){
		push(@{$ttsClustersThreepEnds{$id}}, $line[0], $line[2], $line[5])
	}
}
close F;
print STDERR "Done\n";


print STDERR "Parsing input GFF $transcripts_gff\n";
my %gffHash=gffToHash($transcripts_gff, 'transcript_id', 1, 'exon');
print STDERR "Done\n";

foreach my $tx (keys %gffHash){
	if(defined $fivepsupported{$tx}){
		unless(defined $tssClustersFivepEnds{$tx}){
			die "Couldn't find $tx in $bedOfClusteredFivepEnds.\n";
		}
		adjustTranscriptEnd(5, $tx);
	}
	if(defined $threepsupported{$tx}){
		unless(defined $ttsClustersThreepEnds{$tx}){
			die "Couldn't find $tx in $bedOfClusteredThreepEnds.\n";
		}
		adjustTranscriptEnd(3, $tx);
	}
	unless (defined $fivepsupported{$tx} || defined $threepsupported{$tx})	{ # just output corresponding records without editing anything
	}
}

#print only EXON records

my @outGff=hashToGff(\%gffHash);
print join("", @outGff);

sub adjustTranscriptEnd{
	my $endType=$_[0];
	my $transcript=$_[1];
	my $firstRec=${$gffHash{$transcript}}[0]; #extract first record of transcript id just to deduce strand
	my $strand=$$firstRec[6];
	if( ($strand eq '+' && $endType == 5) || ($strand eq '-' && $endType == 3)){ # sort by increasing start
		foreach my $exonRecord (sort { $a->[3] <=> $b->[3] } @{$gffHash{$transcript}}){ # adjust start coordinates of first exon
			my $adjustedStart;
				if($endType == 5){
					$adjustedStart=${$tssClustersFivepEnds{$transcript}}[1]
				}
				elsif($endType == 3){
					$adjustedStart=${$ttsClustersThreepEnds{$transcript}}[1]
				}
				else{
					die;
				}
				$$exonRecord[3]=$adjustedStart;
			last; #exit after first exon
		}
		#add impossible exons
		foreach my $exonRecord (sort { $a->[3] <=> $b->[3] } @{$gffHash{$transcript}}){
			push(@{$gffHash{$transcript}}, addImpossibleExons($exonRecord, 0)); # 0 =left
			last;
		}

	}
	elsif(($strand eq '+' && $endType == 3) || ($strand eq '-' && $endType == 5)){ # sort by decreasing start
		foreach my $exonRecord (sort { $b->[3] <=> $a->[3] } @{$gffHash{$transcript}}){ # adjust end coordinates of first exon
			my $adjustedEnd;
				if($endType == 5){
					$adjustedEnd=${$tssClustersFivepEnds{$transcript}}[1]
				}
				elsif($endType == 3){
					$adjustedEnd=${$ttsClustersThreepEnds{$transcript}}[1]
				}
				else{
					die;
				}
				$$exonRecord[4]=$adjustedEnd;
			last; #exit after first exon
		}
		foreach my $exonRecord (sort { $b->[3] <=> $a->[3] } @{$gffHash{$transcript}}){
			push(@{$gffHash{$transcript}},addImpossibleExons($exonRecord,1)); # 1 = right
			last;
		}
	}
	else{
		die;
	}
}



sub addImpossibleExons{
	my @eRecord=@{$_[0]};
	my $direction=$_[1]; # 0=left, 1=right
	my $exonSize=1;
	my $intronSize=3;
	my $maxExonNumber=4;
	my @currentExonChain;
	push (@currentExonChain, \@eRecord);
	for (my $i=0; $i<$maxExonNumber; $i++){
		my @currentExon=@{$currentExonChain[$#currentExonChain]};
		$currentExon[8]=dclone(${$currentExonChain[$#currentExonChain]}[8]);
		$exonSize++; #increment size each time so that fake exons of a TSS and a TTS from 2 different transcripts are never compatible (this shit happens for real!)
		my $currentExonStart=$currentExon[3];
		my $currentExonEnd=$currentExon[4];
		my $newStart;
		my $newEnd;
		if($direction == 0){
			$newEnd=$currentExonStart-$intronSize;
			$newStart=$currentExonStart-$intronSize-$exonSize;
		}
		elsif($direction ==1){
			$newStart=$currentExonEnd+$intronSize;
			$newEnd=$currentExonEnd+$intronSize+$exonSize;
		}
		else{
			die;
		}
		$currentExon[3]=$newStart;
		$currentExon[4]=$newEnd;
		${$currentExon[8]}{'fakeExon'}='yes';
		push(@currentExonChain, \@currentExon);
	}
	shift(@currentExonChain); # remove first element of array so that last "real" exon is not duplicated
	return @currentExonChain;
}