#!/usr/bin/perl -w
#输入文件是每一个gene内的信号(bed)格式(由bedgraph转换而来)
#以及每个基因的阈值表#
#升级: 通过样条插值获得最大值的peak的具体点位置
use List::Util qw (max);
use Statistics::R;
my $R = Statistics::R-> new(shared => 1);
my $indir = shift || die $!;
my $ref = shift || die $!; #每一个基因的阈值
#my $m = shift;
#$m = "0" unless ($m) ;
my @in = glob "$indir/*";
open REF, $ref;

my (%ucgThre,%ucgPval);
while ( <REF> ) {
	chomp;
	my ($ucg,$gthreshold,$pval) = (split /\t/, $_)[0,1,2];
	$ucgThre{$ucg} = $gthreshold;
	$ucgPval{$ucg} = $pval;
}
my $fileNth=0;
print "All "; print $#in + 1 ; print " files to be processed\n";

while ($fileNth <= @in -1){
	my $in = $in[$fileNth];
	open IN, $in;
	open TMPOUT, ">/tmp/tmpfile";
	my $threshold = "";
	while (my $line = <IN>){
		chomp ($line);
		my @line = (split /\t/,$line);
		#my ($peakStart,$peakEnd,$peak,$value,$strand,$chr,$ucgStart,$ucgEnd,$ucgSymbol) = @line[1,2,3,4,5,6,7,8,9];
		my $clutValue = $line[4];
		$threshold = $ucgThre{$line[10]};
		
		
		#print $threshold . "\n";
		#if (defined $threshold and $clutValue > $threshold){
		if (defined $threshold and $clutValue > ($threshold /2)){
			print TMPOUT "$line\n";
		}
	}
	#my ($peakStarts,$peakEnds);
	my $k = 0;
	my ($clutStart, $clutEnd);
	my @clutValue = ();
	close IN;
	close TMPOUT;
	#print "$in file larger than threshold signal beds file made\n";
	open IN2, "/tmp/tmpfile";
	my $t = 0;
	while (my $info = <IN2>){
		chomp ($info);$
		$t ++;
		my @infos = (split /\t/,$info);
		my ($peakStart,$peakEnd,$peak,$value,$strand,$chr,$ucgStart,$ucgEnd,$ucgSymbol,$ucgid) = @infos[1,2,3,4,5,6,7,8,9,10];
		my $sigValueCounts = $peakEnd - $peakStart ; 
	
		my ($pkpos,$pkheight,$report) = ("","",""); 

		if ($k ==0){
			do {
				push @clutValue,$value;
				$sigValueCounts --;
			}until ($sigValueCounts == 0);
			$clutStart = $peakStart;
			$clutEnd = $peakEnd;
			$k = 1;
			#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\n" if (eof);
			if (eof){
				if (max(@clutValue) > $threshold ){
					
					push @clutValue, $clutStart;
					push @clutValue, $clutEnd;
					$report = &peak_finder(@clutValue);
					($pkpos,$pkheight) = (split /\t/,$report)[0,1];
					$report = &sub_peak_filter ($report, $threshold);
					pop @clutValue; pop @clutValue;
					print STDERR "$chr\t$clutStart\t$clutEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$ucgStart\t$ucgEnd\t$report\t$ucgid\n" ;
					#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\n" ;
					#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$clutStart\t$clutEnd\t$ucgid\n" if ($m ==1);
				}
			}
			#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$clutStart\t$clutEnd\n" if (eof);



		}
		elsif ($k ==1){
			if ($clutEnd == $peakStart){
				#@clutValue = &repeat_push(@clutValue, $sigValueCounts, $value);
				do {
				push @clutValue,$value;
				$sigValueCounts --;
				}until ($sigValueCounts == 0);
				#push @clutValue,$value;
				$clutEnd =$peakEnd;
				
				#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\n" if (eof);
				#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$clutStart\t$clutEnd\n" if (eof);
				if (eof){
					if (max(@clutValue) > $threshold ){
						push @clutValue, $clutStart;
						push @clutValue, $clutEnd;
					$report = &peak_finder(@clutValue);
					($pkpos,$pkheight) = (split /\t/,$report)[0,1];
					$report = &sub_peak_filter ($report, $threshold);
					
						#($pkpos,$pkheight) = (split /,/,&peak_finder(@clutValue))[0,1];
						pop @clutValue; pop @clutValue;
						print STDERR "$chr\t$clutStart\t$clutEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$ucgStart\t$ucgEnd\t$report\t$ucgid\n" ;
						#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\t$report\t$ucgid\n" ;
						#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$clutStart\t$clutEnd\t$ucgid\n" if ($m ==1);
						#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\n" ;
					}
				}
			}
			elsif ($clutEnd != $peakStart){
				#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue). "\t$strand\t$clutStart\t$clutEnd\n";
				if (max(@clutValue) > $threshold ){
					push @clutValue, $clutStart;
					push @clutValue, $clutEnd;
					$report = &peak_finder(@clutValue);
					($pkpos,$pkheight) = (split /\t/,$report)[0,1];
					$report = &sub_peak_filter ($report, $threshold);
					
					#($pkpos,$pkheight) = (split /,/,&peak_finder(@clutValue))[0,1];
					pop @clutValue; pop @clutValue;
					print STDERR "$chr\t$clutStart\t$clutEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$ucgStart\t$ucgEnd\t$report\t$ucgid\n" ;
					#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\t$report\t$ucgid\n" ;
					#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$clutStart\t$clutEnd\t$ucgid\n"if ($m ==1);
						#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\n" ;
				}	
				#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$clutStart\t$clutEnd\n" ;
				$clutStart = $peakStart;
				$clutEnd = $peakEnd;
				@clutValue =();
				do {
				push @clutValue,$value;
				$sigValueCounts --;
				}until ($sigValueCounts == 0);
				#push @clutValue,$value;
				
				#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue). "\t$strand\t$peakStart\t$peakEnd\n" if (eof);
				if (eof){
					if (max(@clutValue) > $threshold ){
						push @clutValue, $clutStart;
						push @clutValue, $clutEnd;
					$report = &peak_finder(@clutValue);
					($pkpos,$pkheight) = (split /\t/,$report)[0,1];
					$report = &sub_peak_filter ($report, $threshold);
					
					#	($pkpos,$pkheight) = (split /,/,&peak_finder(@clutValue))[0,1];
						pop @clutValue; pop @clutValue;
						#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$clutStart\t$clutEnd\t$ucgid\n"if ($m ==1);
						print STDERR "$chr\t$clutStart\t$clutEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$ucgStart\t$ucgEnd\t$report\t$ucgid\n" ;
						#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\t$report\t$ucgid\n" ;
						#print STDERR "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\t" . max(@clutValue) . "\t$strand\t$clutStart\t$clutEnd\n" ;
					}
				}#print "$chr\t$ucgStart\t$ucgEnd\t$ucgSymbol\tmax(@clutValue)\t$strand\t$peakStart\t$peakEnd\n" if (eof);
			}
		}
	}
	close IN2;
	$fileNth ++;
	
	print "$fileNth have been processed\n" if ($fileNth % 100 ==0);

}

$R->stop();
print "Finish "; print $#in +1; print " files processed, Hooray!!!\n";
sub peak_finder{
#因为perl的子程序里不能把array赋值给array, 所以把起始位置信息之前就已经push在cluster 信号值的末尾两个
#输出为峰值位置,峰值数值高度
my @values = @_;
#my $R = Statistics::R->new();
my $end = pop @values;
my $start = pop @values;
$R->set('start',$start);
$R->set('end',$end);
$R->run(q`x<-seq(start,end-1)`);
$R->set('y',\@values);

my $dfValue=($end - 1 - $start)/10;
if ($dfValue < 3){
	$dfValue =3;
}
$R->set('dfValue',$dfValue);
$R->run(q`f<-smooth.spline(x,y,df=dfValue)`);
$R->run(q`deriv_pred<-predict(f,x,deriv=1)`);
my $clustDerivs=$R->get('deriv_pred$y'); #这是一个array ref, 所以条目还是沿用perl的0,1,2风格

#print $clustDerivs->[78];

my $position =1 ;
my $left = $clustDerivs->[0];
my $peakPos = "";
my $peakValue = "";
my $maxmin = "";
while ($position <= @values - 1 ){
	my $right = $clustDerivs->[$position];
	if ($left * $right <= 0){
		if ($left != $right){
			if (abs($left) > abs($right)){
				$peakValue .= $values[$position] . ",";
				$peakPos .= ($start + $position) . ",";
			}
			else {
				$peakValue .= $values[$position - 1] . ",";
				$peakPos .= ($start + $position -1) . ",";
			}
			if ($left < $right){
				$maxmin .="0".",";
			}
			else {
				$maxmin .= "1".",";
			}
		}
	}
	$left = $right;

	$position ++;

}
my $report = "";
#if ($peakPos && $peakValue){
#	$report = $peakPos . "," . $peakValue;
#}
#else  {
#	$report = "N,N";
#}
#$R->stop();
$report = "$peakPos" . "\t" . "$peakValue" . "\t" . $maxmin;
return $report;
}

sub sub_peak_filter{
	my ($report,$cutoff) = ($_[0],$_[1]);
	my ($realcorr,$realheight) = ("","");
	#峰值坐标, 峰值高度, 峰值是否为极大值
	my ($corr,$height,$maxmin) = (split /\t/,$report)[0,1,2];
	my @corr = split /,/,$corr;
	my @height = split /,/,$height;
	my @maxmin = split /,/,$maxmin;
	for (my $i = 0; $i <= @maxmin -1 ; $i ++){
		if ($maxmin[$i] == 1){
			#说明是极大值
			if ($height[$i] > $cutoff){
				$realcorr .= $corr[$i] . ",";
				$realheight .= $height[$i] . ",";
			}
		}
	}
	return ($realcorr . "\t" . $realheight );
}
=head
sub repeat_push {
	my (@array, $count, $charactor) =@_;
	do {
		push @array, $charactor;
		$count --;
	}until ($count == 0);
	return (@array);
}
=cut

