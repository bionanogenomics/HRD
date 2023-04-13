#!/usr/bin/perl -w 

use strict;
use warnings;
use Getopt::Long;

my $initialCommand = "$0 @ARGV";
print STDOUT "\nInfo: Running the command: $initialCommand\n";

my ($smapFile, $cnvFile, $outFile) = ("", "", "");
my $cnvMaskFile = "";
my $centromereFile = "";
my $chromLengthFile = "";
my $aneuploidyFile = "";

my $minLohSize = 15000000;	# default 15 Mbp for losses and AOHs
my $minTaiSize = 3000000;	# default 3 Mbp for Gains, Losses and AOH from telomere
my $minLstSize = 10000000;	# default 10 Mbp for Gains, Losses, insertions, inversions, AOH
my $maxLstDistance = 3000000;	# default to 3 Mbp from any neighbouring bkpt
my $numChromothripsisFusion = 15;	# number of fusion events to define a chromosome as chromothripsis chromosome

my $telomericWindow = 500000;	# default allows +/- 500 kbp offset from the telomere and still be consider telomeric (also use this window to search for known CNV mask regions)
my $centromericWindow = 500000;	# default allows +/- 500 kbp offset from the centromere and still be consider centromeric (also use this window to search for known CNV mask regions)

my $aneuploidyScore = 0.95;	# default 0.95, the recommended aneuploidy confidence score

my $maxVarFrequency = 0.00;	# default 0, that is, not seen in Bionano control sample SV database
my $indelConfScore = 0.00;	# default 0 (solve3.7)
my $inversionConfScore = 0.7;	# default 0.7 (solve3.7)
my $duplicationConfScore = -1;	# default -1 (solve3.7)
my $translocationConfScore = 0.05;	# default 0.05 (solve3.7)
my $selfMoleculeCheckFlag = "yes";	# default pass self molecule check

my $cnvScore = 0.99;	# default 0.99 (solve3.7)
my $cnvStitchDistance = 500000;	# default based on postnatal study, replication analysis of kh24 (a 29 mbp duplication fragmented)
 

GetOptions	(
	"smapFile=s"			=>	\$smapFile,
	"cnvFile=s"			=>	\$cnvFile,
	"outFile=s"			=>	\$outFile,
	"centromereFile=s"		=>	\$centromereFile,
	"chromLengthFile=s"		=>	\$chromLengthFile,
	"aneuploidyFile=s"		=>	\$aneuploidyFile,
	"cnvMaskFile=s"			=>	\$cnvMaskFile,
	"minLohSize=i"			=> 	\$minLohSize,
	"minTaiSize=i"			=>	\$minTaiSize,
	"telomericWindow=i"		=>	\$telomericWindow,
	"centromericWindow=i"		=>	\$centromericWindow,
	"minLstSize=i"			=>	\$minLstSize,
	"maxLstDistance=i"		=>	\$maxLstDistance,
	"maxVarFrequency=f"		=>	\$maxVarFrequency,
	"numChromothripsisFusion=i"	=>	\$numChromothripsisFusion,
	"aneuploidyScore=f"		=>	\$aneuploidyScore,
	"indelConfScore=f"		=>	\$indelConfScore,
	"inversionConfScore=f"		=>	\$inversionConfScore,
	"duplicationConfScore=f"	=>	\$duplicationConfScore,
	"translocationConfScore=f"	=>	\$translocationConfScore,
	"selfMoleculeCheckFlag=s"	=>	\$selfMoleculeCheckFlag,
	"cnvScore=s"			=>	\$cnvScore,
	"cnvStitchDistance=i"		=>	\$cnvStitchDistance
) or die "ERROR: $0: error in command line arguments.\n";

die "Please enter in input smap file\n" if ($smapFile eq "");
die "Please enter in input cnv file\n" if ($cnvFile eq "");
die "Please enter in input aneuploidy file\n" if ($aneuploidyFile eq "");
die "Please enter in the output file\n" if ($outFile eq "");
die "Please enter in the centromere cytoband file\n" if ($centromereFile eq "");
die "Please enter in the chromosome length file\n" if ($chromLengthFile eq "");

my $chromLengthsRef = getChromLengths($chromLengthFile);
# read in a list of noisy regions (aka Mask regions) of the genome
my $noisyRegionsRef = getCnvMasks($cnvMaskFile, $cnvStitchDistance, $chromLengthsRef);

# extend telomeres including the neighbouring noisy regions
my $extendedTelomeresRef = extendTelomeres($chromLengthsRef, $telomericWindow, $noisyRegionsRef);
# extend centromeres to include the neighbouring noisy regions
my %centromeres = ();
my $centromeresRef = getCentromeres($centromereFile, \%centromeres);
my $extendedCentromeresRef = extendCentromeres($centromeresRef, $centromericWindow, $noisyRegionsRef);

# read in aneuploidy
my $aneuploidyRef = getAneuploidy($aneuploidyFile, $aneuploidyScore);
# read in the SV file first to scan for chromosomes with a significant number of intrachromosomal fusion events
my $chromothripsisChromosomesRef = findChromothripsisChromosomes($smapFile, $numChromothripsisFusion, $extendedCentromeresRef, $maxVarFrequency, $translocationConfScore, $selfMoleculeCheckFlag);

# read in SV and CNV
my $svRef = getSv($smapFile, $aneuploidyRef, $extendedCentromeresRef, $extendedTelomeresRef, $minLohSize, $minTaiSize, $minLstSize, 
	$maxVarFrequency, $indelConfScore, $inversionConfScore, $duplicationConfScore, $translocationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef);
$svRef = getCnv($cnvFile, $aneuploidyRef, $extendedCentromeresRef, $extendedTelomeresRef, $minLohSize, $minTaiSize, $minLstSize, $maxVarFrequency, $cnvScore, $cnvStitchDistance, $chromothripsisChromosomesRef, $svRef);

### compute scores ###
# for LST, we need to record the breakpoints (but distinguish whether they are loh, tai or lst - only lst bkpt can contribute to lst score in the end)
my $bkptsRef = findBkpts($svRef);
$bkptsRef = sortByCoord($bkptsRef, "position", "ascend");

my $lstScoresRef = calculateLstScores($bkptsRef, $maxLstDistance, $chromLengthsRef);
my $lohScoresRef = calculateLohScores($svRef, $chromLengthsRef);
my $taiScoresRef = calculateTaiScores($svRef, $chromLengthsRef);

## output 
printOutput($outFile, $chromLengthsRef, $lstScoresRef, $lohScoresRef, $taiScoresRef, $chromothripsisChromosomesRef, $numChromothripsisFusion);

sub printOutput	{
	my ($file, $chromLengthsRef, $lstScoresRef, $lohScoresRef, $taiScoresRef, $chromothripsisChromosomesRef, $numChromothripsisFusion) = @_;

	open(OUT, ">$file") or die "ERROR: printOutput: cannot write to $file: $!\n";
	my $count = 0;
	print STDOUT "printOutput: writing to $file\n";

	print OUT "chrom\tLST\tLOH\tTAI\tnumIntrachromFusion_chromothripsis\n";

	my @chromInOrder = (1..24);
	foreach my $chrom (@chromInOrder)	{

		my $outLine = join("\t", $chrom, $lstScoresRef->{$chrom}, $lohScoresRef->{$chrom}, $taiScoresRef->{$chrom})."\t";
		$outLine .= ($chromothripsisChromosomesRef->{$chrom} > $numChromothripsisFusion) ? ("$chromothripsisChromosomesRef->{$chrom}*") : ("$chromothripsisChromosomesRef->{$chrom}");
		print OUT "$outLine\n";
		$count += 1;

	} # foreach chrom

	print STDOUT "\twritten $count records\n";
	close OUT;
} # printOutput

sub calculateTaiScores	{
	my ($svRef, $chromLengthsRef) = @_;

	my %taiScores = ();
	# initialize
	foreach my $chrom (keys %$chromLengthsRef)	{
		$taiScores{$chrom} = 0;
	} # foreach chrom

	foreach my $chrom (keys %{$svRef->{tai}})	{
		$taiScores{$chrom} = scalar( @{$svRef->{tai}{$chrom}} );
	} # foreach chrom

	return \%taiScores
} # calculateTaiScores

sub calculateLohScores	{
	my ($svRef, $chromLengthsRef) = @_;

	my %lohScores = ();

	# initialize
	foreach my $chrom (keys %$chromLengthsRef)	{
		$lohScores{$chrom} = 0;
	} # foreach chrom

	foreach my $chrom (keys %{$svRef->{loh}})	{
		for (my $i = 0; $i < scalar(@{$svRef->{loh}{$chrom}}); $i += 1)	{
			# each SV gets one point
			$lohScores{$chrom} += 1;
		} # for i
	} # foreach chrom

	return \%lohScores;
} # calculateLohScores

sub calculateLstScores	{
	my ($bkptsRef, $maxLstDistance, $chromLengthsRef) = @_;
	my %lstScores = ();

	# initialize lstScores
	foreach my $chrom (keys %$chromLengthsRef)	{
		$lstScores{$chrom} = 0;
	} # foreach chrom

	foreach my $chrom (keys %$bkptsRef)	{

		# now along the chromosome
		# each bkpt encountered add 1 to the score, unless it has a noScoreFlag
		# then look at the next bkpt to see if that is within distance from the current bkpt, if yes, mark the next bkpt noScoreFlag

		for (my $i = 0; $i < scalar(@{$bkptsRef->{$chrom}}); $i += 1)	{
			my $curBkptRef = $bkptsRef->{$chrom}[$i];

			if ($curBkptRef->{noScoreFlag} == 0)	{
				$lstScores{$chrom} += 1;
			} # if curBkptRef			
		
			if ($i < $#{$bkptsRef->{$chrom}})	{
				my $nextBkptRef = $bkptsRef->{$chrom}[$i + 1];
				if ( abs($nextBkptRef->{position} - $curBkptRef->{position}) <= $maxLstDistance )	{
					# if next bkpt is too close
					$nextBkptRef->{noScoreFlag} = 1;
				} # if nextBkptRef
			} # if i	
		} # for i
	} # foreach chrom	
	return \%lstScores;
} # calculateLstScores

sub sortByCoord	{
	my ($bkptsRef, $coord, $direction) = @_;

	foreach my $chrom (keys %$bkptsRef)	{
		if ($direction eq "ascend")	{
			@{$bkptsRef->{$chrom}}	= sort	{
				$a->{$coord}	<=>	$b->{$coord}
			} @{$bkptsRef->{$chrom}}
		} else	{
			@{$bkptsRef->{$chrom}} = sort	{
				$b->{$coord}	<=>	$a->{$coord}
			} @{$bkptsRef->{$chrom}}
		} # if direction
	} # foreach chrom

	return $bkptsRef;
} # sortByCoord

sub findBkpts	{
	my ($svRef) = @_;

	my %bkpts = ();

	foreach my $hrdType (keys %$svRef)	{
		foreach my $chrom (keys %{$svRef->{$hrdType}})	{

			for (my $i = 0; $i < scalar(@{$svRef->{$hrdType}{$chrom}}); $i += 1)	{
				my $sRef = $svRef->{$hrdType}{$chrom}[$i];
				if ($hrdType eq "lst")	{
					# only lst type variant
					# first breakpoint
					push(@{$bkpts{$chrom}}, {position => $sRef->{start}, hrdType => $hrdType, id => $sRef->{id}, noScoreFlag => 0});
					# second breakpoint
					push(@{$bkpts{$sRef->{chrom2}}}, {position => $sRef->{end}, hrdType => $hrdType, id => $sRef->{id}, noScoreFlag => 0});
				} # hrdType
			} # for i

		} # foreach chrom
	} # foreach hrdType

	return \%bkpts;
} # findBkpts


##########################################################################################################
### CNV mask ###

sub getCnvMasks	{
	my ($file, $cnvStitchDistance, $chromLengthsRef) = @_;
	my %noisyRegions = ();
	my %tempNoisyRegions = ();
	
	if ($file eq "")	{
		# if no cnv mask (aka noisy cnv regions) file has been provided, then initialize and return
		foreach my $chrom (keys %$chromLengthsRef)	{
			push(@{$noisyRegions{$chrom}{noisyRegion}}, {id => -999999, start => -999999, end => -999999, size => 0, confidence => 1, line => ""});
		} # foreach chrom
		return \%noisyRegions;
	} # if file
	
	open(IN, $file) or die "ERROR: getCnvMasks: cannot read in $file: $!\n";
	my $count = 0;
	print STDOUT "getCnvMasks: reading in $file\n";
	my %headerContent = ();
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		
		# empty line
		next if ($line =~ /^$/);
		
		# header line
		if ($line =~ /^chr/i)	{
			my @content = split(/\t/, $line);
			
			for (my $i = 0; $i < scalar(@content); $i += 1)	{
				$headerContent{$content[$i]} = $i;
			} # for i
			next;
		} # if line
		
		
		# data line
		my @content = split(/\t/, $line);
		my ($chrom, $start, $end) = ($content[$headerContent{Chr}], int($content[$headerContent{StartPos}]), int( $content[$headerContent{EndPos}]) );
		($start, $end) = ($end, $start) if ($end < $start);
		my $size = $end - $start + 1;
		
		$count += 1;
		push(@{$tempNoisyRegions{$chrom}{noisyRegion}}, {id => $count, start => $start, end => $end, size => $size, svType => "noisyRegion", confidence => 1, line => "$line", nearByFlag => 0});

	} # while line
	print STDOUT "\tread in $count records\n";
	close IN;
	
	# sort by coordinates
	my $tempNoisyRegionsRef = sortByCnvCoord(\%tempNoisyRegions, "start");	
	
	my $noisyRegionsRef = stitchCnvs($tempNoisyRegionsRef, $cnvStitchDistance);	

	return $noisyRegionsRef;
} # getCnvMasks


############## CNV #########################

sub getCnv	{
	my ($file, $aneuploidyRef, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $minLstSize, $maxVarFrequency, $cnvScore, $cnvStitchDistance, $chromothripsisChromosomesRef, $svRef) = @_;
	
	open(IN, $file) or die "ERROR: cannot open $file: $!\n";
	my $count = 0;
	print STDOUT "getCnv: reading in $file\n";
	my %headerContent = ();
	my %tempCnvs = ();		# temporary stores CNV until stitching is done
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		
		# empty line
		next if ($line =~ /^$/);
		
		if ($line =~ /^#/)	{
			if ($line =~ /^#\s*Id/i)	{
				$line =~ s/^#\s*//;
				
				my @content = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@content); $i += 1)	{
					$headerContent{$content[$i]} = $i;
				} # for i
			} # if line
			next;
		} # if line
		
		# data lines
		my @content = split(/\t/, $line);
		my ($id, $chrom, $start, $end, $size, $svType, $confidence) = ($content[$headerContent{Id}], $content[$headerContent{Chromosome}], $content[$headerContent{Start}], $content[$headerContent{End}], $content[$headerContent{Width}], $content[$headerContent{Type}], $content[$headerContent{Confidence}]);
		
		# score
		next if ($confidence < $cnvScore);
		
		($start, $end) = (int($start), int($end));
		($start, $end) = ($end, $start) if ($end < $start);
		
		my $type = ($svType =~ /gain/i) ? ("dup") : ( ($svType =~ /loss/i) ? ("del") : ("") );
		next if ($type eq "");

		push(@{$tempCnvs{$chrom}{$type}}, {id => $id, start => $start, end => $end, size => $size, svType => $svType, confidence => $confidence, line => $line, nearByFlag => 0});
		$count += 1;
	} # while line

	# sort by coordinates
	my $tempCnvsRef = sortByCnvCoord(\%tempCnvs, "start");

	# stitch (assuming that neighbouring CNVs do not overlap)
	# note that this does not take into consideration of masking 
	# (lone masked CNV, would not be stitched, and would not be included in HRD)
	# if post stitched CNV consisted of ONLY masked call, also, will not go into HRD
	# if mask is mixed with non-masked, then the whole region would not have been masked, thus justifying the stitch, will continue to HRD
	my $stitchedCnvsRef = stitchCnvs($tempCnvsRef, $cnvStitchDistance);

	# filter
	foreach my $chrom (keys %$stitchedCnvsRef)	{
		# skip if this chromosome has chromothripsis signature
		# next if (exists $chromothripsisChromosomesRef->{$chrom});
		
		foreach my $type (keys %{$stitchedCnvsRef->{$chrom}})	{
			# check aneuploidy, disregard CNV calls that constitute to an aneuploidy call
			next if ( (exists $aneuploidyRef->{$chrom} && $aneuploidyRef->{$chrom}{type} =~ /del/ && $type =~ /del/) ||
				(exists $aneuploidyRef->{$chrom} && $aneuploidyRef->{$chrom}{type} =~ /dup/ && $type =~ /dup/));

			# parse
			for (my $i = 0; $i < scalar(@{$stitchedCnvsRef->{$chrom}{$type}}); $i += 1)	{
				my $iRef = $stitchedCnvsRef->{$chrom}{$type}[$i];
				if ($type =~ /del/)	{
					$svRef = parseDeletion($iRef->{line}, $chrom, $iRef->{start}, $iRef->{end}, $iRef->{id}, $type, $iRef->{size}, $iRef->{confidence}, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $svRef);
				} else	{
					$svRef = parseDuplication($iRef->{line}, $chrom, $iRef->{start}, $iRef->{end}, $iRef->{id}, $type, $iRef->{confidence}, $iRef->{size}, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $svRef);
				} # if type
			} # for i

		} # foreach type
	} # foreach chrom
	
	return $svRef;
} # getCnv

sub stitchCnvs	{
	my ($tempCnvsRef, $cnvStitchDistance) = @_;
	my %stitchedCnvs = ();

	# stitch near by neighbouring CNVs of the same type
	# note that this stitch does not consider whether a call has been masked - because masking leads to fragmentation and not premature termination

	foreach my $chrom (keys %$tempCnvsRef)	{
		foreach my $type (keys %{$tempCnvsRef->{$chrom}})	{
			my $theIndex = 0;
			my ($theStart, $theEnd) = ($tempCnvsRef->{$chrom}{$type}[$theIndex]{start}, $tempCnvsRef->{$chrom}{$type}[$theIndex]{end});
			my $i = 1;
			for (; $i < scalar(@{$tempCnvsRef->{$chrom}{$type}}); $i += 1)	{
				my $iRef = $tempCnvsRef->{$chrom}{$type}[$i];

				if ($iRef->{start} <= $theEnd + $cnvStitchDistance)	{
					$iRef->{nearByFlag} = 1;
					$theEnd = $iRef->{end};
					next;
				} # if close by

				# not close by to "theIndex" call
				# make a new record of variants starting from "theIndex"
				# reset "theIndex" call
				my ($numMaskedCalls, $line) = (0, "");
				for (my $j = $theIndex; $j <= $i - 1; $j += 1)	{
					my $jRef = $tempCnvsRef->{$chrom}{$type}[$j];
					$numMaskedCalls += 1 if ($jRef->{svType} =~ /mask/i);
					$line .= "$jRef->{line}\n";
				} # for j
				$line =~ s/\n$//;

				# now see if all the calls to be stitched had been masked, if so, do not stitch, otherwise, stitch
				if ($numMaskedCalls < ($i - 1) - $theIndex + 1)	{
					# stitch
					my $firstCallRef = $tempCnvsRef->{$chrom}{$type}[$theIndex]; # first call
					my $lastCallRef = $tempCnvsRef->{$chrom}{$type}[$i - 1]; # first call
					push(@{$stitchedCnvs{$chrom}{$type}}, {id => $firstCallRef->{id}, start => $firstCallRef->{start}, end => $lastCallRef->{end}, 
						size => $lastCallRef->{end} - $firstCallRef->{start} + 1, confidence => $firstCallRef->{confidence}, line => $line});
					
				} else	{
					# do not stitch, and throw out the masked calls
				} # if numMaskedCalls
				($theIndex, $theStart, $theEnd) = ($i, $iRef->{start}, $iRef->{end});
			
			} # for i
			my ($numMaskedCalls, $line) = (0, "");
			for (my $j = $theIndex; $j <= $i - 1; $j += 1)	{
				my $jRef = $tempCnvsRef->{$chrom}{$type}[$j];
				$numMaskedCalls += 1 if ($jRef->{svType} =~ /mask/i);
				$line .= "$jRef->{line}\n";
			} # for j
			$line =~ s/\n$//;
			
			if ($numMaskedCalls < ($i - 1) - $theIndex + 1)	{
				# stitch, since not all calls had been masked
				my $firstCallRef = $tempCnvsRef->{$chrom}{$type}[$theIndex]; # first call
				my $lastCallRef = $tempCnvsRef->{$chrom}{$type}[$i - 1]; # first call
				push(@{$stitchedCnvs{$chrom}{$type}}, {id => $firstCallRef->{id}, start => $firstCallRef->{start}, end => $lastCallRef->{end},
					size => $lastCallRef->{end} - $firstCallRef->{start} + 1, confidence => $firstCallRef->{confidence}, line => $line});
			} # if numMaskedCalls

		} # foreach type
	} # foreach chrom

	return \%stitchedCnvs;
} # stitchCnvs

sub sortByCnvCoord	{
	my ($tempCnvsRef, $coord) = @_;

	foreach my $chrom (keys %$tempCnvsRef)	{
		foreach my $type (keys %{$tempCnvsRef->{$chrom}})	{
			@{$tempCnvsRef->{$chrom}{$type}} = sort	{
				$a->{$coord}	<=>	$b->{$coord}
			} @{$tempCnvsRef->{$chrom}{$type}}
		} # foreach type
	} # foreach chrom

	return $tempCnvsRef;
} # sortByCoord

##########################################################################################################
############### SV #########################

sub getSv	{
	my ($file, $aneuploidyRef, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $minLstSize, $maxVarFrequency, $indelConfScore, $inversionConfScore, $duplicationConfScore, $translocationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef) = @_;
	my %sv = ();
	my $svRef = \%sv;
	
	$svRef = getIndels($file, $aneuploidyRef, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $minLstSize, $maxVarFrequency, $indelConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef);

	$svRef = getDuplications($file, $aneuploidyRef, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $minLstSize, $maxVarFrequency, $duplicationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef);

	$svRef = getInversions($file, $centromeresRef, $telomeresRef, $minLstSize, $maxVarFrequency, $inversionConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef);
	
	$svRef = getTranslocations($file, $minLstSize, $maxVarFrequency, $translocationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef);
	
	return $svRef;
} # getSv

################ Insertion/Deletion #######################

sub getIndels	{
	my ($file, $aneuploidyRef, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $minLstSize, $maxVarFrequency, $indelConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef) = @_;	
	
	open(IN, $file) or die "ERROR: getIndels: cannot open file: $!\n";
	my $count = 0;
	print STDOUT "getIndels: reading in $file\n";
	my %headerContent = ();
	
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		
		# empty lines
		next if ($line =~ /^$/);
		
		if ($line =~ /^#/)	{
			if ($line =~ /^#h/)	{
				$line =~ s/^#h\s+//;
				my @content = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@content); $i += 1)	{
					$headerContent{$content[$i]} = $i;
				} # for i
			} # if line
			# just skip header line
			next;
		} # if line
		
		my @content = split(/\t/, $line);
		my ($smapId, $chrom, $start, $end, $confidence, $svType, $size) = ($content[$headerContent{SmapEntryID}], $content[$headerContent{RefcontigID1}], $content[$headerContent{RefStartPos}], $content[$headerContent{RefEndPos}], $content[$headerContent{Confidence}], $content[$headerContent{Type}], $content[$headerContent{SVsize}]);
		my ($percentDb, $percentDbSameEnzyme, $selfMoleculeCheck) = ($content[$headerContent{'Present_in_%_of_BNG_control_samples'}], $content[$headerContent{'Present_in_%_of_BNG_control_samples_with_the_same_enzyme'}], $content[$headerContent{'Found_in_self_molecules'}]);
		
		$svRef = filterSvInsertionDeletion($line, $chrom, $start, $end, $smapId, $svType, $size, $confidence, $percentDb, $percentDbSameEnzyme, $selfMoleculeCheck, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $maxVarFrequency, $indelConfScore, $selfMoleculeCheckFlag, $aneuploidyRef, $chromothripsisChromosomesRef, $svRef);

		
		$count += 1;
	} # while line
	print STDOUT "\tread in $count records\n";
	close IN;
	
	return $svRef
} # getIndels

sub filterSvInsertionDeletion	{
	my ($line, $chrom, $start, $end, $id, $type, $size, $confidence, $percentDb, $percentDbSameEnzyme, $selfMoleculeCheck, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $maxVarFrequency, $indelConfScore, $selfMoleculeCheckFlag, $aneuploidyRef, $chromothripsisChromosomesRef, $svRef) = @_;

	# insertion and deletion, no tiny
	return $svRef if (! ($type =~ /insertion|deletion/ && $type !~ /tiny/) );
	# nbase
	return $svRef if ($type =~ /insertion_nbase|deletion_nbase/);


	# filter by control sample SV database, if available
	# if not available, then keep
	return $svRef if ($percentDb > $maxVarFrequency || $percentDbSameEnzyme > $maxVarFrequency);

	# score
	return $svRef if ($confidence < $indelConfScore);

	# self molecule check
	return $svRef if (! ($selfMoleculeCheck =~ /^$selfMoleculeCheckFlag$/i) );
	

	# disregard CNV calls that constitute to an aneuploidy call
	#return $svRef if ( (exists $aneuploidyRef->{$chrom} && $aneuploidyRef->{$chrom}{type} =~ /del/ && $type =~ /del/) ||
	#	(exists $aneuploidyRef->{$chrom} && $aneuploidyRef->{$chrom}{type} =~ /dup/ && $type =~ /dup/) );

	# skip if this chromosome has chromothripsis signature
	# return $svRef if (exists $chromothripsisChromosomesRef->{$chrom});	

	($start, $end) = (int($start), int($end));	
	($start, $end) = ($end, $start) if ($end < $start);		
	
	if ($type =~ /del/)	{
		$type = "del";
		$svRef = parseDeletion($line, $chrom, $start, $end, $id, $type, $size, $confidence, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $svRef);
	} else	{
		$type = "ins";
		$svRef = parseInsertion($line, $chrom, $start, $end, $id, $type, $size, $confidence, $minLstSize, $centromeresRef, $telomeresRef, $svRef);
	} # if type
		
	return $svRef;
} # filterSvInsertionDeletion

sub parseInsertion	{
	my ($line, $chrom, $start, $end, $id, $type, $size, $confidence, $minLstSize, $centromeresRef, $telomeresRef, $svRef) = @_;	

	# insertion >10MB, exclusive of intervals <3MB (not accounting for centromeric breaks)
	if ($size > $minLstSize)	{
		# interstitial
		if ( ($telomeresRef->{$chrom}{p}{end} < $start && $end < $centromeresRef->{$chrom}{start}) ||
			($centromeresRef->{$chrom}{end} < $start && $end < $telomeresRef->{$chrom}{q}{start}) )	{
			# 2 breaks
			# store temporarily, will need to break into breakpoints to count the number of breaks
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST i INS interstitial: $line\n";
			return $svRef;
		} # if interstitial
		
		# whole chromosome arm
		if ( ($start <= $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) || 
			( $centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $telomeresRef->{$chrom}{q}{start}  <= $end ) )	{
			# no op, because technically there is no chrom break outside the centromere
			return $svRef;
		} # if whole chromosome arm

		# centromeric
		if ( $start > $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) {
			# no break, centromeric breaks do not count
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST j INS p centromeric: $line\n";
			return $svRef;
		} # if p centromeric
		if ($centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $end < $telomeresRef->{$chrom}{q}{start}) {
			# 1 break: end
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST k INS q centromeric: $line\n";
			return $svRef;
		} # if q centromeric 

		# telomeric 
		if ( $start <= $telomeresRef->{$chrom}{p}{end} && $end < $telomeresRef->{$chrom}{q}{start} )       {
			# 1 break: end
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST l INS p telomeric: $line\n";
			return $svRef;
		} # if p telomeric
		if ( $telomeresRef->{$chrom}{p}{end} < $start && $telomeresRef->{$chrom}{q}{start} <= $end )    {
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST m INS q telomeric: $line\n";
			return $svRef;
		} # if q telomeric

		# across centromere
		if ($telomeresRef->{$chrom}{p}{end} < $start && $start < $centromeresRef->{$chrom}{start} && $centromeresRef->{$chrom}{end} < $end && $end < $telomeresRef->{$chrom}{q}{start})  {
			# not possible, no op
			return $svRef;
		} # if across centromere
		
		# whole chromosome
		if ($start < $telomeresRef->{$chrom}{p}{end} && $telomeresRef->{$chrom}{q}{start} < $end)  {
			# not possible, no op
			return $svRef;
		} # while chromosome
	} # if size
	
	return $svRef;
} # parseInsertion

sub parseDeletion	{
	my ($line, $chrom, $start, $end, $id, $type, $size, $confidence, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $svRef) = @_;
	
	# deletion
	# interstitial
	if ( ($telomeresRef->{$chrom}{p}{end} < $start && $end < $centromeresRef->{$chrom}{start}) ||
		($centromeresRef->{$chrom}{end} < $start && $end < $telomeresRef->{$chrom}{q}{start}) )	{

		if ($size > $minLstSize)	{
			# 2 breaks
			# store temporarily, will need to break into breakpoints to count the number of breaks
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST a DEL interstitial: $line\n";
		} # if size

		if ($size > $minLohSize)	{
			# LOH: # regions losses & AOH > 15 Mbp (not whole chromosome)
			push(@{$svRef->{loh}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LOH a DEL intersitital: $line\n";
		} # if size
		return $svRef;
	} # if interstitial

	# whole chromosome arm
	if ( ($start <= $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) ||
		( $centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $telomeresRef->{$chrom}{q}{start} <= $end ) )       {
		# no TAI when a breakpoint overlaps the centromere
		#if ($size > $minTaiSize)      {
		#	# get a score for TAI
		#	push(@{$svRef->{tai}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
		#	print STDOUT "TAI a DEL whole chromosome arm: $line\n";
		#} # if size	

		if ($size > $minLohSize)	{
			# LOH
			push(@{$svRef->{loh}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LOH b DEL whole chromosome arm: $line\n";	
		} # if size

		# no break, no score for LST
		return $svRef;
	} # if whole chromosome arm

	# centromeric
	if ( $start > $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) {
		if ($size > $minLstSize)	{
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST b DEL p centromeric: $line\n";
		} # if size

		if ($size > $minLohSize)	{
			# LOH
			push(@{$svRef->{loh}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LOH c DEL p centromeric: $line\n";	
		} # if size
		return $svRef;
	} # if p centromeric
	if ($centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $end < $telomeresRef->{$chrom}{q}{start} ) {
		if ($size > $minLstSize)	{
			# 1 break: end
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
                        print STDOUT "LST c DEL q centromeric: $line\n";
		} # if size

		if ($size > $minLohSize)	{
			# LOH
			push(@{$svRef->{loh}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LOH d DEL q centromeric: $line\n";	
		} # if size
		return $svRef;
	} # if q centromeric

	# telomeric
	if ( $start <= $telomeresRef->{$chrom}{p}{end} && $end < $telomeresRef->{$chrom}{q}{start} )       {
		if ($size > $minTaiSize && ($end < $centromeresRef->{$chrom}{start}))	{
			# regions Gains, Losses and AOH > 3 Mbp from telomere (not crossing the centromere)
			push(@{$svRef->{tai}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "TAI b DEL p telomeric: $line\n";	
		} # if size

		if ($size > $minLstSize)	{
			# 1 break: end
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST d DEL p telomeric: $line\n";
		} # if size

		if ($size > $minLohSize)	{
			# LOH
			push(@{$svRef->{loh}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LOH e DEL p telomeric: $line\n";
		} # if size

		return $svRef;
	} # if p telomeric
	if ( $telomeresRef->{$chrom}{p}{end} < $start && $telomeresRef->{$chrom}{q}{start} <= $end )    {
		if ($size > $minTaiSize && ($centromeresRef->{$chrom}{end} < $start))        {
			# regions Gains, Losses and AOH > 3 Mbp from telomere (not crossing the centromere)
			push(@{$svRef->{tai}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "TAI c DEL q telomeric: $line\n";
		} # if size

		if ($size > $minLstSize)	{
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST e DEL q telomeric: $line\n";
		} # if size

		if ($size > $minLohSize)	{
			# LOH
			push(@{$svRef->{loh}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LOH f DEL q telomeric: $line\n";
		} # if size	
		return $svRef;
	} # if q telomeric

	# across centromere
	if ($telomeresRef->{$chrom}{p}{end} < $start && $start < $centromeresRef->{$chrom}{start} && $centromeresRef->{$chrom}{end} < $end && $end < $telomeresRef->{$chrom}{q}{start} )  {

		if ($size > $minLstSize)	{
			# 2 breaks
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST f DEL across centromere: $line\n";
		} # if size

		if ($size > $minLohSize)	{
			# LOH
			push(@{$svRef->{loh}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LOH g DEL across centromere: $line\n";
		} # if size
		return $svRef;
	} # if across centromere
	
	# whole chromosome
	if ($start < $telomeresRef->{$chrom}{p}{end} && $telomeresRef->{$chrom}{q}{start} < $end)  {
		# no op
		return $svRef;
	} # if whole chromosome

	return $svRef;
} # parseDeletion

########### Duplication #################

sub getDuplications	{
	my ($file, $aneuploidyRef, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $minLstSize, $maxVarFrequency, $duplicationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef) = @_;
	
	open(IN, $file) or die "ERROR: getDuplications: cannot open $file: $!\n";
	my $count = 0;
	print STDOUT "getDuplications: reading in $file\n";
	my %headerContent = ();
	
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		
		# empty lines
		next if ($line =~ /^$/);
		
		# header line 
		if ($line =~ /^#/)	{
			if ($line =~ /^#h/)	{
				$line =~ s/^#h\s+//;
				my @content = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@content); $i += 1)	{
					$headerContent{$content[$i]} = $i;
				} # for i
			} # if line
			# just skip header line
			next;
		} # if line
		
		# data line
		my @content = split(/\t/, $line);
		my ($smapId, $chrom, $start, $end, $confidence, $svType, $size) = ($content[$headerContent{SmapEntryID}], $content[$headerContent{RefcontigID1}], $content[$headerContent{RefStartPos}], $content[$headerContent{RefEndPos}], $content[$headerContent{Confidence}], $content[$headerContent{Type}], $content[$headerContent{SVsize}]);
		my ($percentDb, $percentDbSameEnzyme, $selfMoleculeCheck) = ($content[$headerContent{'Present_in_%_of_BNG_control_samples'}], $content[$headerContent{'Present_in_%_of_BNG_control_samples_with_the_same_enzyme'}], $content[$headerContent{'Found_in_self_molecules'}]);
		
		$svRef = filterSvDuplication($line, $chrom, $start, $end, $smapId, $svType, $confidence, $size, $percentDb, $percentDbSameEnzyme, $selfMoleculeCheck, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $maxVarFrequency, $duplicationConfScore, $selfMoleculeCheckFlag, $aneuploidyRef, $chromothripsisChromosomesRef, $svRef);	
		
		$count += 1;		
	} # while line
	close IN;
	print STDOUT "\tread in $count records\n";
	return $svRef;
} # getDuplications

sub filterSvDuplication	{
	my ($line, $chrom, $start, $end, $id, $type, $confidence, $size, $percentDb, $percentDbSameEnzyme, $selfMoleculeCheck, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $maxVarFrequency, $duplicationConfScore, $selfMoleculeCheckFlag, $aneuploidyRef, $chromothripsisChromosomesRef, $svRef) = @_;
	
	return $svRef if ($type !~ /duplication/);
	$type = "dup";

	# filter by control sample SV database, if available
        # if control sample SV database frequency is not available, then keep
        return $svRef if ($percentDb > $maxVarFrequency || $percentDbSameEnzyme > $maxVarFrequency);
            
        # score
        return $svRef if ( $confidence < $duplicationConfScore );
        
        # molecule check
        return $svRef if (! ($selfMoleculeCheck =~ /^$selfMoleculeCheckFlag$/i) );
        
	# disregard CNV calls that constitute to an aneuploidy call
	#return $svRef if ( (exists $aneuploidyRef->{$chrom} && $aneuploidyRef->{$chrom}{type} =~ /del/ && $type =~ /del/) ||
	#	(exists $aneuploidyRef->{$chrom} && $aneuploidyRef->{$chrom}{type} =~ /dup/ && $type =~ /dup/) );

	# skip if this chromosome has chromothripsis signature
	# return $svRef if (exists $chromothripsisChromosomesRef->{$chrom});	        
        
	($start, $end) = (int($start), int($end));
        ($start, $end) = ($end, $start) if ($end < $start);

	$svRef = parseDuplication($line, $chrom, $start, $end, $id, $type, $confidence, $size, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $svRef);
	
	return $svRef;
} # filterSvDuplications

sub parseDuplication	{
	my ($line, $chrom, $start, $end, $id, $type, $confidence, $size, $minLstSize, $centromeresRef, $telomeresRef, $minLohSize, $minTaiSize, $svRef) = @_;

	# duplication
	# interstitial
	if ( ($telomeresRef->{$chrom}{p}{end} < $start && $end < $centromeresRef->{$chrom}{start} ) ||
		($centromeresRef->{$chrom}{end} < $start && $end < $telomeresRef->{$chrom}{q}{start}) )        {
		# check if pass LST size (10 Mbp)
		if ($size > $minLstSize)	{
			# 2 breaks
			# store temporarily, will need to break into breakpoints to count the number of breaks
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST k DUP interstitial: $line\n";
		} # if size
		return $svRef;
	} # if interstitial

	# whole chromosome arm
	if ( ($start <= $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) ||
		( $centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $telomeresRef->{$chrom}{q}{start} <= $end ) )       {
		# no TAI, if a breakpoint overlaps the centromere
		#if ($size > $minTaiSize)	{
		#	# get a score for TAI
		#	push(@{$svRef->{tai}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
		#	print STDOUT "TAI: $line\n";
		#} # if size > minTaiSize

		# no score for LST, no break
		return $svRef;
	} # if whole chromosome arm

	# centromeric
	if ( $start > $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) {
		if ($size > $minLstSize)	{
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST l DUP p centromeric: $line\n";
		} # if size
		return $svRef;
	} # if p centromeric
	if ($centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $end < $telomeresRef->{$chrom}{q}{start} ) {
		if ($size > $minLstSize)	{
                       	# 1 break: end
                  	push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
                       	print STDOUT "LST m DUP q centromeric: $line\n";
		} # if size
		return $svRef;
  	} # if q centromeric

	# telomeric
	if ( $start <= $telomeresRef->{$chrom}{p}{end} && $end < $telomeresRef->{$chrom}{q}{start} )       {
		if ( $size > $minTaiSize && ($end < $centromeresRef->{$chrom}{start}) )	{
			# regions Gains, Losses and AOH > 3 Mbp from telomere (not crossing the centromere)
			push(@{$svRef->{tai}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "TAI a DUP: $line\n";
		} # if size > minTaiSize
			
		if ($size > $minLstSize)	{
			# 1 break: end
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST n DUP p telomeric: $line\n";
		} # if size
		return $svRef;
	} # if p telomeric
	if ( $telomeresRef->{$chrom}{p}{end} < $start && $telomeresRef->{$chrom}{q}{start} <= $end )    {
		if ($size > $minTaiSize && ($centromeresRef->{$chrom}{end} < $start) )	{
			# regions Gains, Losses and AOH > 3 Mbp from telomere (not crossing the centromere)
			push(@{$svRef->{tai}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "TAI b DUP: $line\n";
		} # if size > $minTaiSize

		if ($size > $minLstSize)	{
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST o DUP q telomeric: $line\n";
		} # if size
		return $svRef;
	} # if q telomeric

	# across centromere
	if ($telomeresRef->{$chrom}{p}{end} < $start && $start < $centromeresRef->{$chrom}{start} && $centromeresRef->{$chrom}{end} < $end && $end < $telomeresRef->{$chrom}{q}{start})  {
		if ($size > $minLstSize)	{
			# 2 breaks
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST p DUP across centromere: $line\n";
		} # if size
		return $svRef;
	} # if across centromere

	# whole chromosome
	if ($start < $telomeresRef->{$chrom}{p}{end} && $telomeresRef->{$chrom}{q}{start} < $end)	{
		# no op, because technically there is no chrom break
		return $svRef;
	} # if whole chromosome

	return $svRef;
} # parseDuplication

############ Inversions ###################

sub getInversions	{
	my ($file, $centromeresRef, $telomeresRef, $minLstSize, $maxVarFrequency, $inversionConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef) = @_;
	
	my %invTemp = ();
	open(IN, $file) or die "ERROR: getInversions: cannot open $file: $!\n";
	my $count = 0;
	print STDOUT "getInversions: reading in $file\n";
	my %headerContent = ();
	
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		
		# empty lines
		next if ($line =~ /^$/);
		
		# header lines
		if ($line =~ /^#/)	{
			if ($line =~ /^#h/)	{
				$line =~ s/^#h\s+//;
				my @content = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@content); $i += 1)	{
					$headerContent{$content[$i]} = $i;
				} # for i
			} # if line
			# just skip header line
			next;
		} # if line
		
		# data line
		my @content = split(/\t/, $line);
		my ($smapId, $chrom, $start, $end, $score, $svType, $linkSmapId) = ($content[$headerContent{SmapEntryID}], $content[$headerContent{RefcontigID1}], $content[$headerContent{RefStartPos}], $content[$headerContent{RefEndPos}], $content[$headerContent{Confidence}], $content[$headerContent{Type}], $content[$headerContent{LinkID}]);
		my ($percentDb, $percentDbSameEnzyme, $selfMoleculeCheck) = ($content[$headerContent{'Present_in_%_of_BNG_control_samples'}], $content[$headerContent{'Present_in_%_of_BNG_control_samples_with_the_same_enzyme'}], $content[$headerContent{'Found_in_self_molecules'}]);
		
		next if (! ($svType =~ /inversion/i) );
		my $excludeFlag = ($svType =~ /overlap/i || $svType =~ /duplicate/i ) ? (1) : (0);	# exclude such inversion types 
		my $type = "inv";
		
		my $size = -1;
		my $id = ($smapId < $linkSmapId) ? ("$smapId\t$linkSmapId") : ("$linkSmapId\t$smapId");
		($start, $end) = (int($start), int($end));
		
		# combines two lines (per inversion) into one (in other words, flattening the inversion)
		if (! exists $invTemp{$id})	{
			# if this id pair has not been encountered, populates positions 1 and 2
			$invTemp{$id} = {chrom => $chrom, type => $type, position1 => $start, position2 => $end, position3 => -2, position4 => -2, size => $size, score => $score, percentDb => $percentDb, percentDbSameEnzyme => $percentDbSameEnzyme, selfMoleculeCheck => $selfMoleculeCheck, line => $line, excludeFlag => $excludeFlag};
		} else	{
			# populates positions 3 and 4
			($invTemp{$id}{position3}, $invTemp{$id}{position4}) = ($start, $end);
			$invTemp{$id}{excludeFlag} += $excludeFlag;	# remember whether to exclude, as soon as the value is > 0
			$invTemp{$id}{line} .= "\n$line";	# concatenate two lines
		} # if exists
		
		$count += 1;
	} # while line
	
	print STDOUT "\tread in $count records\n";
	close IN;
	
	
	# now determine the outer most non-negative co-ordinates, thus flattening
	foreach my $id (keys %invTemp)	{
		my $iRef = $invTemp{$id};

		# skip those flagged as exclude
		next if ($iRef->{excludeFlag} > 0);
		
		my ($id1, $id2) = split(/\t/, $id);	# note that id1 will be recorded
		
		my ($outerStart, $outerEnd) = getOuterCoords($iRef->{position1}, $iRef->{position2}, $iRef->{position3}, $iRef->{position4});

		if ($outerStart < 0 || $outerEnd < 0 || $outerEnd < $outerStart)	{
			die "ERROR: getInversions: cannot determine the outerStart and outerEnd position for $iRef->{position1}, $iRef->{position2}, $iRef->{position3}, $iRef->{position4} where variant id1=$id1, variant id2=$id2\n";
		} # if outerStart
		
		# recalculate the size of the inversion
		$iRef->{size} = $outerEnd - $outerStart + 1;

		$svRef = parseInversion($iRef->{line}, $iRef->{chrom}, $outerStart, $outerEnd, $id, $iRef->{type}, $iRef->{score}, $iRef->{size}, $iRef->{percentDb}, $iRef->{percentDbSameEnzyme}, $iRef->{selfMoleculeCheck}, $minLstSize, $centromeresRef, $telomeresRef, $maxVarFrequency, $inversionConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef);		
		
	} # foreach id
	
	return $svRef;
} # getInversions

sub getOuterCoords	{
	my ($pos1, $pos2, $pos3, $pos4) = @_;
	my @theArray = ($pos1, $pos2, $pos3, $pos4);
	my ($outerStart, $outerEnd) = (-999999, -999999);
	for (my $i = 0; $i < scalar(@theArray); $i += 1)	{

		if ($outerStart == -999999 || ($theArray[$i] > 0 && $theArray[$i] < $outerStart) )	{
			$outerStart = $theArray[$i];
		} # if 

		if ($outerEnd == -999999 || ($theArray[$i] > 0 && $outerEnd < $theArray[$i]))	{
			$outerEnd = $theArray[$i];
		} # if 
		
	} # for i
	return ($outerStart, $outerEnd);
} # getOuterCoords

sub parseInversion	{
	my ($line, $chrom, $start, $end, $id, $type, $score, $size, $percentDb, $percentDbSameEnzyme, $selfMoleculeCheck, $minLstSize, $centromeresRef, $telomeresRef, $maxVarFrequency, $inversionConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef) = @_;
	
	# score
	return $svRef if ($score < $inversionConfScore);
	
	# percent db
	return $svRef if ($percentDb > $maxVarFrequency || $percentDbSameEnzyme > $maxVarFrequency);
	
	# self molecule check
	return $svRef if (! ($selfMoleculeCheck =~ /^$selfMoleculeCheckFlag$/i) );
	
	# skip if this chromosome has chromothripsis signature
	# return $svRef if (exists $chromothripsisChromosomesRef->{$chrom});	

	($start, $end) = ($end, $start) if ($end < $start);

	# inverions >10MB, exclusive of intervals <3MB (not accounting for centromeric breaks)
	if ($size > $minLstSize)	{
		# interstitial 
		if ( ($telomeresRef->{$chrom}{p}{end} < $start && $end < $centromeresRef->{$chrom}{start}) || 
			($centromeresRef->{$chrom}{end} < $start && $end < $telomeresRef->{$chrom}{q}{start}) )	{
			# 2 breaks
			# store temporarily, will need to break into breakpoints to count the number of breaks
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST c inv interstitial: $line\n";
			return $svRef;
		} # if interstitial

		# whole chromosome arm
		if ( ($start <= $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) || 
			( $centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $telomeresRef->{$chrom}{q}{start} <= $end ) )	{
			# no op, because technically there is no chrom break outside the centromere
			return $svRef;
		} # if whole chromosome arm

		# centromeric
		if ( $start > $telomeresRef->{$chrom}{p}{end} && $centromeresRef->{$chrom}{start} <= $end && $end < $centromeresRef->{$chrom}{end} ) {
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST d inv p centromeric: $line\n";
			return $svRef;
		} # if p centromeric
		if ($centromeresRef->{$chrom}{start} < $start && $start <= $centromeresRef->{$chrom}{end} && $end < $telomeresRef->{$chrom}{q}{start})	{
			# 1 break: end
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST e inv q centromeric: $line\n";
			return $svRef;
		} # if q centromeric

		# telomeric 
		if ( $start <= $telomeresRef->{$chrom}{p}{end} && $end < $telomeresRef->{$chrom}{q}{start} )	{
			# 1 break: end
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $end, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST f inv p telomeric: $line\n";
			return $svRef;
		} # if p telomeric
		if ( $telomeresRef->{$chrom}{p}{end} < $start && $telomeresRef->{$chrom}{q}{start} <= $end )	{
			# 1 break: start
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $start, size => $size, id => $id});
			print STDOUT "LST g inv q telomeric: $line\n";
			return $svRef;
		} # if q telomeric

		# across centromere
		if ($telomeresRef->{$chrom}{p}{end} < $start && $start < $centromeresRef->{$chrom}{start} && $centromeresRef->{$chrom}{end} < $end && $end < $telomeresRef->{$chrom}{q}{start})	{
			# 2 breaks
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom, end => $end, size => $size, id => $id});
			print STDOUT "LST h inv pericentric: $line\n";
			return $svRef;
		} # if across centromere

		# whole chromosome
		if ($start < $telomeresRef->{$chrom}{p}{end} && $telomeresRef->{$chrom}{q}{start} < $end)	{
			# no op, because technically there is no chrom break
			return $svRef;
		} # if whole chromosome
	} # if size

	return $svRef;
} # parseInversion

######### Translocations #############

sub getTranslocations	{
	my ($file, $minLstSize, $maxVarFrequency, $translocationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef) = @_;
	
	open(IN, $file) or die "getTranslocations: cannot open $file: $!\n";
	my $count = 0;
	print STDOUT "getTranslocations: reading in $file\n";
	my %headerContent = ();	
	
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;
		
		# empty lines
		next if ($line =~ /^$/);
		
		# header lines
		if ($line =~ /^#/)      {
			if ($line =~ /^#h/)     {
				$line =~ s/^#h\s+//;
				my @content = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@content); $i += 1) {
					$headerContent{$content[$i]} = $i;
				} # for i
			} # if line

			# just skip header line
			next;
		} # if line		
		
		# data lines
		my @content = split(/\t/, $line);
		my ($smapId, $chrNum1, $chrNum2, $start1, $start2, $score, $svType) = ($content[$headerContent{SmapEntryID}], $content[$headerContent{RefcontigID1}], $content[$headerContent{RefcontigID2}], $content[$headerContent{RefStartPos}], $content[$headerContent{RefEndPos}], $content[$headerContent{Confidence}], $content[$headerContent{Type}]);
		my ($percentDb, $percentDbSameEnzyme, $selfMoleculeCheck) = ($content[$headerContent{'Present_in_%_of_BNG_control_samples'}], $content[$headerContent{'Present_in_%_of_BNG_control_samples_with_the_same_enzyme'}], $content[$headerContent{'Found_in_self_molecules'}]);
		
		next if ($svType !~ /trans/i);	# only keep translocation calls
		# filter out those that are not simple translocations
		next if (! ($svType =~ /translocation_interchr|trans_interchr/ || $svType =~ /translocation_intrachr|trans_intrachr/));
		# whether to include translocation common and segdup
		next if ($svType =~ /common|segdup/);
		my $type = "tra";

		( $start1, $start2 ) = ( int($start1), int($start2) );
		
		$svRef = parseTranslocation($line, $chrNum1, $chrNum2, $start1, $start2, $smapId, $type, $score, $percentDb, $percentDbSameEnzyme, $selfMoleculeCheck, $minLstSize, $maxVarFrequency, $translocationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef);
		
		$count += 1;
	} # while line
	
	print STDOUT "\tread in $count records\n";
	close IN;
	
	return $svRef;
} # getTranslocations

sub parseTranslocation	{
	my ($line, $chrom, $chrom2, $start, $end, $id, $type, $score, $percentDb, $percentDbSameEnzyme, $selfMoleculeCheck, $minLstSize, $maxVarFrequency, $translocationConfScore, $selfMoleculeCheckFlag, $chromothripsisChromosomesRef, $svRef) = @_;
	
	die "ERROR: parseTranslocation: cannot parse end position for translocation variant = $line\n" if ($end == -1);

	# score
	return $svRef if ($score < $translocationConfScore);

	# filter by control sample SV database, if available
	# if control sample SV database frequency is not available, then keep
	return $svRef if ($percentDb > $maxVarFrequency || $percentDbSameEnzyme > $maxVarFrequency);
	
	# molecule check
	return $svRef if (! ($selfMoleculeCheck =~ /^$selfMoleculeCheckFlag$/i) );
	
	# skip if this chromosome has chromothripsis signature
	# return $svRef if (exists $chromothripsisChromosomesRef->{$chrom});

	my $size = -1;
	

	# LST: for intrachromosomal translocation/fusion, make sure that the start and end are > 10 Mbp apart
	if ($chrom eq $chrom2)  {
		# intrachromosomal fustion/translocation
		($start, $end) = ($end, $start) if ($end < $start);
		if ($end - $start >= $minLstSize)       {
			push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom2, end => $end, size => $size, id => $id});
			print STDOUT "LST a: $line\n";
		} # if the two positions are far apart
	} else  {
		# inter chromosomal translocation
		
		# skip if chrom2 has chromothripsis signature
		#return $svRef if (exists $chromothripsisChromosomesRef->{$chrom2});

		push(@{$svRef->{lst}{$chrom}}, {type => $type, start => $start, chrom2 => $chrom2, end => $end, size => $size, id => $id});
		print STDOUT "LST b: $line\n";
	} # if chrom

	return $svRef;
} # parseTranslocation

sub findChromothripsisChromosomes	{
	my ($file, $numChromothripsisFusion, $centromeresRef, $maxVarFrequency, $translocationConfScore, $selfMoleculeCheckFlag) = @_;
	my %chromothripsisChromosomesCount = ();

	# initialize
	foreach my $chrom (keys %$centromeresRef)	{
		$chromothripsisChromosomesCount{$chrom} = 0;
	} # foreach chrom

	my %headerContent = ();
	open(IN, $file) or die "ERROR: findChromothripsisChromosomes: cannot open $file: $!\n";
	my $count = 0; 
	print STDOUT "findChromothripsisChromosomes: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;

		# empty lines
		next if ($line =~ /^$/);
		
		# header lines
		if ($line =~ /^#/)      {
			if ($line =~ /^#h/)     {
				$line =~ s/^#h\s+//;
				my @content = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@content); $i += 1) {
					$headerContent{$content[$i]} = $i;
				} # for i
			} # if line

			# just skip header line
			next;
		} # if line	

		# data lines
		my @content = split(/\t/, $line);
		my ($smapId, $chrNum1, $chrNum2, $start1, $start2, $score, $svType) = ($content[$headerContent{SmapEntryID}], $content[$headerContent{RefcontigID1}], $content[$headerContent{RefcontigID2}], $content[$headerContent{RefStartPos}], $content[$headerContent{RefEndPos}], $content[$headerContent{Confidence}], $content[$headerContent{Type}]);
		my ($percentDb, $percentDbSameEnzyme, $selfMoleculeCheck) = ($content[$headerContent{'Present_in_%_of_BNG_control_samples'}], $content[$headerContent{'Present_in_%_of_BNG_control_samples_with_the_same_enzyme'}], $content[$headerContent{'Found_in_self_molecules'}]);
		
		next if ($svType !~ /trans/i);	# only keep translocation calls
		# filter out those that are not simple translocations
		next if (! ($svType =~ /translocation_interchr|trans_interchr/ || $svType =~ /translocation_intrachr|trans_intrachr/));
		# whether to include translocation common and segdup
		next if ($svType =~ /common|segdup/);

		( $start1, $start2 ) = ( int($start1), int($start2) );

		# only high confidence variant
		next if ($score < $translocationConfScore);
		
		# ctrl db
		next if ($percentDb > $maxVarFrequency || $percentDbSameEnzyme > $maxVarFrequency );
		
		# self molecule check
		next if (! ($selfMoleculeCheck =~ /^$selfMoleculeCheckFlag$/i) );
		

		# only look at intrachromosomal fusion
		next if (! ($chrNum1 =~ /^$chrNum2$/) );

		# now record the number of intrachromosomal fusion
		$chromothripsisChromosomesCount{$chrNum1} += 1;
		$count += 1;
	} # while line
	print STDOUT "\tread in $count records\n";
	close IN;

	# now iterate the count and see if that chromosome can be considered as chromothriptic
	my %chromothripsisChromosomes = ();
	foreach my $chrom (keys %chromothripsisChromosomesCount)	{
		$chromothripsisChromosomes{$chrom} = $chromothripsisChromosomesCount{$chrom};
		
		
		if ($chromothripsisChromosomesCount{$chrom} >= $numChromothripsisFusion)	{
			print STDOUT "\t$chrom has chromothripsis with $chromothripsisChromosomesCount{$chrom} fusions\n";
		} # if chromothripsisChromosomes
		
	} # freach chrom

	return \%chromothripsisChromosomes;
} # findChromothripsisChromosomes

sub getAneuploidy	{
	my ($file, $aneuploidyScore) = @_;
	my %aneuploidy = ();

	open(IN, $file) or die "ERROR: getAneuploidy: cannot open $file: $!\n";
	my $count = 0;

	print STDOUT "getAneuploidy: reading in $file\n";
	my %headerContent = ();
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;

		# empty lines
		next if ($line =~ /^$/);

		# header lines
		if ($line =~ /^#/)	{
			if ($line =~ /^#chr/)	{
				$line =~ s/^#//;
				my @content = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@content); $i += 1)	{
					$headerContent{$content[$i]} = $i;
				} # for i
			} # if line
			next;
		} # if line

		# data line
		my @content = split(/\t/, $line);
		my ($chromNum, $type, $fractChrLen, $score, $fractCN) = ($content[$headerContent{chr}], $content[$headerContent{types}], $content[$headerContent{fractChrLen}], $content[$headerContent{score}], ,$content[$headerContent{fractCN}]);

		$type = ($type eq "gain") ? ("DUP") : ("DEL");
		$type = lc ($type);	# lower case

		# skip of the event is low confidence
		next if ($score < $aneuploidyScore);
		
		$aneuploidy{$chromNum} = {type => $type, fractChrLen => $fractChrLen, score => $score, fractCN => $fractCN};
		$count += 1;	
	} # while line

	print STDOUT "\tread in $count records\n";
	close IN;

	return \%aneuploidy;
} # getAneuploidy

################ centromeres #########################
sub extendCentromeres	{
	my ($centromeresRef, $window, $noisyRegionsRef) = @_;
	
	my %extendedCentromeres = ();
	
	foreach my $chrom (keys %$centromeresRef)	{
		my ($theStart, $theEnd) = ($centromeresRef->{$chrom}{start} - $window, $centromeresRef->{$chrom}{end} + $window);
		my $oIndeciesRef = searchForOverlap($theStart, $theEnd, $noisyRegionsRef->{$chrom}{noisyRegion} );
		
		if (scalar(@$oIndeciesRef) > 0)	{
			my ($leftMostPos, $rightMostPos) = findLeftRightMostPos($noisyRegionsRef->{$chrom}{noisyRegion}, $oIndeciesRef);
			$theStart = $leftMostPos if ($leftMostPos < $theStart);
			$theEnd = $rightMostPos if ($theEnd < $rightMostPos);
		} # if scalar
		
		$extendedCentromeres{$chrom} = {start => $theStart, end => $theEnd};
	} # foreach chrom
	
	return \%extendedCentromeres;
} # extendCentromeres

sub getCentromeres      {
	my ($file, $centromeresRef) = @_;

	open(IN, $file) or die "ERROR: getCentromeres: cannot open $file: $!\n";
	my %tempCentromeres = ();
	my $count = 0;
	print STDOUT "getCentromeres: reading in $file\n";
	while (my $line = <IN>) {
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;

		my @content = split(/\t/, $line);
		my ($chrom, $start, $end, $band, $name) = @content;

		# take only the typical chromosomes
		next if ($chrom =~ /ran|chrM|hap/i);
		my $chromNum = $chrom;
		$chromNum =~ s/^chr//;
		$chromNum = ($chromNum =~ /^X$/i) ? (23) : ( ($chromNum =~ /^Y$/i ) ? (24) : ($chromNum) );

		# only record the centromeres
		next if ($name ne "acen");

		# fix co-ordinate if end < start
		($start, $end) = ($end, $start) if ($end < $start);

		push(@{$tempCentromeres{$chromNum}}, {start => $start, end => $end});

		$count += 1;
	} # while line
	print STDOUT "\tread in $count records\n";
	close IN;

	# now combine all bands into one cytoband per chromosome
	foreach my $chrom (keys %tempCentromeres)       {
		@{$tempCentromeres{$chrom}} = sort      {
			$a->{start}     <=>     $b->{start}
		} @{$tempCentromeres{$chrom}};

		my ($theStart, $theEnd) = ($tempCentromeres{$chrom}[0]{start}, $tempCentromeres{$chrom}[$#{$tempCentromeres{$chrom}}]{end});

		$centromeresRef->{$chrom} = {start => $theStart, end => $theEnd};
	} # foreach chrom

	return $centromeresRef;
} # getCentromeres        

########## overlap ##########
sub searchForOverlap	{
	my ($theStart, $theEnd, $dataChromRef) = @_;
	
	my @oIndecies = ();
	
	for (my $i = 0; $i < scalar(@$dataChromRef); $i += 1)	{
		my $iRef = $dataChromRef->[$i];
		
		if ( ($theStart <= $iRef->{start} && $iRef->{start} <= $theEnd) || ($iRef->{start} <= $theStart && $theStart <= $iRef->{end}) )	{
			push(@oIndecies, $i);
		} # if overlap
	} # for i
	
	return \@oIndecies;
} # searchForOverlap

sub findLeftRightMostPos	{
	my ($dataChromRef, $oIndeciesRef) = @_;
	
	my ($leftMostPos, $rightMostPos) = ($dataChromRef->[$oIndeciesRef->[0]]{start}, $dataChromRef->[$oIndeciesRef->[0]]{end});
	
	for (my $o = 1; $o < scalar(@$oIndeciesRef); $o += 1)	{
		my $theDataRef = $dataChromRef->[$oIndeciesRef->[$o]];
		$leftMostPos = $theDataRef->{start} if ($theDataRef->{start} < $leftMostPos);
		$rightMostPos = $theDataRef->{end} if ($rightMostPos < $theDataRef->{end});
	} # for o
	
	return ($leftMostPos, $rightMostPos);
} # findLeftRightMostPos

############## telomeres ###################
sub extendTelomeres	{
	my ($chromLengthsRef, $window, $noisyRegionsRef) = @_;
	
	# find overlaps between telomeres and the noisy regions
	# then find the leftmost or rightmost co-ordinates
	
	my %extendedTelomeres = ();
	
	foreach my $chrom (keys %$chromLengthsRef)	{
		
		# p ter
		my ($pStart, $pEnd) = (1, 1 + $window);
		my $oIndeciesRef = searchForOverlap($pStart, $pEnd, $noisyRegionsRef->{$chrom}{noisyRegion} );
		if (scalar(@$oIndeciesRef) > 0)	{
			my ($leftMostPos, $rightMostPos) = findLeftRightMostPos($noisyRegionsRef->{$chrom}{noisyRegion}, $oIndeciesRef);
			$pEnd = $rightMostPos if ($pEnd < $rightMostPos);
		} # if scalar 
		
		# q ter
		my ($qStart, $qEnd) = ($chromLengthsRef->{$chrom} - $window, $chromLengthsRef->{$chrom});
		my $o2IndeciesRef = searchForOverlap($qStart, $qEnd, $noisyRegionsRef->{$chrom}{noisyRegion} );
		if (scalar(@$o2IndeciesRef) > 0)	{
			my ($leftMostPos, $rightMostPos) = findLeftRightMostPos($noisyRegionsRef->{$chrom}{noisyRegion}, $o2IndeciesRef);
			$qStart = $leftMostPos if ($leftMostPos < $qStart);
		} # if scalar
		
		$extendedTelomeres{$chrom}{p} = {start => $pStart, end => $pEnd};
		$extendedTelomeres{$chrom}{q} = {start => $qStart, end => $qEnd};
	} # foreach chrom
	
	return \%extendedTelomeres;
} # extendTelomeres

sub getChromLengths	{
	my ($file) = @_;
	my %chromLengths = ();

	open(IN, $file) or die "ERROR: getChromLengths: reading in $file: $!\n";
	my $count = 0;
	print STDOUT "getChromLengths: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r+//g;
		$line =~ s/^\s+|\s+$//g;

		# empty lines
		next if ($line =~ /^$/);

		# header lines
		next if ($line =~ /^#/);

		# data lines
		my ($chromNum, $chromLength) = split(/\t/, $line);

		$chromLengths{$chromNum} = $chromLength;

		$count += 1;
	} # while line
	print STDOUT "\tread in $count records\n";
	close IN;
	return \%chromLengths;
} # getChromLengths
