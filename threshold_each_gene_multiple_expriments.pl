#!/usr/bin/perl -w

use List::Util qw(max min);
my $indir = shift || die $! ;
my $outdir = shift || die $!;
my $experiments_count = shift || die $! ;
my @files = glob "$indir/*";
print "## All " ;
print $#files + 1;
print " to be processed\n";
print "## Files from $indir\n";
print "## Reports to $outdir\n";
print "## $experiments_count times random samples for each file\n";
my $fileNth = 0;
while ($fileNth <= @files - 1){
	my $in = $files[$fileNth];
	open IN,$in;

	my ($gName,$gStart,$gEnd,$tagStart,$tagEnd);
	my (@tagStartS,@tagEndS);
	while ( <IN> ) {
		chomp;
		my @info = (split /\t/,$_);
		($gStart,$gEnd,$tagStart,$tagEnd) = ($info[7],$info[8],$info[1],$info[2]);	
		push (@tagStartS, $tagStart);
		push (@tagEndS,	$tagEnd);
		$gName=$info[10];	
	}
	#my $outRep = $outdir . "/" . $gName . ".report";
	#open OUT,">$outRep";
	my $tagsInGeneNumSum=@tagStartS;
	die $! if (@tagStartS != @tagEndS);

	my $experiment_Nth = $experiments_count;
	my @gMaxs;
	do {
		my (@gSigDis); #单次实验,基因内随机reads num/base
		my $i = 0; 
		while ($i <= $gEnd-$gStart){
				$gSigDis[$i]=0;
				$i ++;
		}
		
		my $tagNth = 0; #每一个tag随机位置
		while ($tagNth <= (@tagStartS - 1 ) ) {
			my $tagLength = $tagEndS[$tagNth] - $tagStartS[$tagNth] + 1;
			my $gRandRange = $gEnd - $gStart - $tagLength + 1;
			my ($randTagStart,$randTagEnd) ;
			while (1){
				$randTagStart = int(rand($gRandRange+1));
				$randTagEnd = $randTagStart + $tagLength -1;
				if ($randTagEnd <= $gEnd){
					last;
				}
			}

			my $j = $randTagStart; #每一个tag赋值基因位置
			while ($j <= $randTagEnd){
				$gSigDis[$j] ++;
				$j ++;
			}

			$tagNth ++;
		}

		my $maxRandPeak = max(@gSigDis);
		#print "$experiment_Nth,$maxRandPeak\n";
		push (@gMaxs,$maxRandPeak);
		#print $maxRandPeak . "\n";
		$experiment_Nth --;
	} until ($experiment_Nth == 0 ); # count 100,99,98 ... 1, sum 100 times
#print @gMaxs;
#print "\n";
	die $! if (@gMaxs != $experiments_count);

	my %count;
	my @gThresholds = (grep {++$count{$_} < 2} @gMaxs);
	@gThresholds = sort {$b <=> $a} @gThresholds;
#print $gMaxThreshold . "\n";
	my @gMaxsDis ;#数组序号是threshold 序号对应的则是出现的次数(概率)累加
	my $indicator = 0;
#print "$experiments_count\n";
#print "Threshold\tSumCounts(>=Threshold)\n";
	my @realCumSumProb;
	while ($indicator <= (@gThresholds -1)){
		#$gMaxsDis[$threshold] = 0;
		my $tmp = $gThresholds[$indicator];
		$gMaxsDis[$indicator] = (grep /^$tmp$/,@gMaxs) ;
		if ($indicator != 0){
			$gMaxsDis[$indicator] = ($gMaxsDis[$indicator] +$gMaxsDis[$indicator-1] ) ;
		}
		#print "$gThresholds[$indicator]\t$gMaxsDis[$indicator]\n";
		$realCumSumProb[$indicator] = $gMaxsDis[$indicator] / $experiments_count;
#=head
		if ($realCumSumProb[$indicator] > 0.05){
			my $output = $gThresholds[$indicator] + 1;
			print  "$gName\t$output\t$realCumSumProb[$indicator-1]\t$gThresholds[0]\t$realCumSumProb[0]\t$tagsInGeneNumSum\n";
			print STDERR "$gName\t$output\t$realCumSumProb[$indicator-1]\t$gThresholds[0]\t$realCumSumProb[0]\t$tagsInGeneNumSum\n";
			#print OUT "$gName\t$output\t$realCumSumProb[$indicator-1]\t$gThresholds[0]\t$realCumSumProb[0]\t$tagsInGeneNumSum\n";
			#$gThresholds[$indicator]\t$realCumSumProb\t$tagsInGeneNumSum\n";
			close IN;
			#close OUT;
			last;
		}
#=cut
		$indicator ++;
	}

$fileNth ++;
print "## $fileNth files have been processed\n" if ($fileNth % 1000 == 0);
}

print "## End processing, Thank you\n";

__DATA__
[Usage] 
$ perl threshold_each_gene_multiple_expriments.pl <in_files_dir> <out_report_dir> <random_tims_for_each_file> 
[Note] <in_files_dir> should not store other unrelevant files
[Author] Yun YAN
[Data] Feb/25/2013

PTL&TMG
