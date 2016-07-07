use strict;

my %seqs2;
my $r1_file=$ARGV[0];
my $r2_file=$ARGV[1];
my $r1_design=$ARGV[2];
my $seq_batch=$ARGV[3];
my $oligos_fn=$ARGV[4];
my $out_fastq_file=$ARGV[5];
my $out_stats_file=$ARGV[6];
my $scdb_path=$ARGV[7];
my $i5_design;
if ($#ARGV>7){
  $i5_design=$ARGV[8];
}
my $index;
my $amp_batches_fn = $scdb_path."/annotations/amp_batches.txt";
my @plate_barcode_intervals=();
my @plate_barcode_i5_intervals=();
my @mRNA_intervals=();

my %column_num;
my %well_barcode_intervals_hash;
my %RMT_intervals_hash;
my %polyT_intervals_hash;
my %pool_barcodes;

open(AMP_BATCHES_FILE, $amp_batches_fn) || die "cannot open file $amp_batches_fn to read.\n";
my $line=<AMP_BATCHES_FILE>;
chomp $line;
my @column_names=split("\t",$line);
for (my $i=0;$i<=$#column_names;$i++){
  $column_num{$column_names[$i]}=$i;
}
print "pool_barcodes:\n";
while (<AMP_BATCHES_FILE>){
  chomp;
  my @row= split("\t");
  my $cur_seq_batch=$row[$column_num{"Seq_batch_ID"}];
  if ($cur_seq_batch eq $seq_batch){
	my @well_barcode_intervals=();
	my @RMT_intervals=();
	my @polyT_intervals=();

	my $pool_barcode=$row[$column_num{"Pool_barcode"}];
	$pool_barcodes{$pool_barcode}=1;
	print STDOUT $pool_barcode."\n";
	my $cur_r2_design=$row[$column_num{"R2_design"}];
	my @r2_split= split(/\./,$cur_r2_design);
	
	my $cur_loc=0;
	for (my $ri=0;$ri<=$#r2_split;$ri++){
	  my $field_type=substr($r2_split[$ri],length($r2_split[$ri])-1);
	  my $field_size=int substr($r2_split[$ri],0,length($r2_split[$ri])-1);
	  my @interval=($cur_loc,$field_size);
	  if ($field_type eq "W"){
		@well_barcode_intervals=(@well_barcode_intervals,@interval);
	  }
	  if ($field_type eq "R"){
		@RMT_intervals=(@RMT_intervals,@interval);
	  }
	  if ($field_type eq "T"){
		@polyT_intervals=(@polyT_intervals,@interval);
	  }
	  $cur_loc+=$field_size;
	}
	
	$well_barcode_intervals_hash{$pool_barcode}=\@well_barcode_intervals;
	@RMT_intervals_hash{$pool_barcode}=\@RMT_intervals;
	@polyT_intervals_hash{$pool_barcode}=\@polyT_intervals;
  }
}

my @column_names=split("\t",<INDEX_FILE>);
for (my $i=0;$i<=$#column_names;$i++){
  $column_num{$column_names[$i]}=$i;
}

my @r1_split= split(/\./,$r1_design);
my $cur_loc=0;
for (my $ri=0;$ri<=$#r1_split;$ri++){
  my $field_type=substr($r1_split[$ri],length($r1_split[$ri])-1);
  my $field_size=int substr($r1_split[$ri],0,length($r1_split[$ri])-1);
  my @interval=($cur_loc,$field_size);
  if ($field_type eq "P"){
	@plate_barcode_intervals=(@plate_barcode_intervals,@interval);
	print $interval[0]."\t".$interval[1]."\n";
  }
  if ($field_type eq "M"){
	@mRNA_intervals=(@mRNA_intervals,@interval);
  }
  $cur_loc+=$field_size;
}

my @i5_split= split(/\./,$i5_design);
my $cur_loc=0;
for (my $ri=0;$ri<=$#i5_split;$ri++){
  my $field_type=substr($i5_split[$ri],length($i5_split[$ri])-1);
  my $field_size=int substr($i5_split[$ri],0,length($i5_split[$ri])-1);
  my @interval=($cur_loc,$field_size);
  if ($field_type eq "P"){
	@plate_barcode_i5_intervals=(@plate_barcode_i5_intervals,@interval);
	print $interval[0]."\t".$interval[1]."\n";
  }
  $cur_loc+=$field_size;
}
my $k=7;
my @oligos=();
my %kmers_hash={};
my $oi=0;
open(OLIGOS_FILE, $oligos_fn) || die "cannot open file $oligos_fn to read.\n";
while (<OLIGOS_FILE>){
  chomp;
  my @line= split("=");
  my $oligo_name=$line[0];
  $oligos[$oi]=$oligo_name;
  my $oligo_seq=$line[1];
  $oligo_seq =~ s/^\s+|\s+$//g;
  for (my $si=0;$si<length($oligo_seq)-$k+1;$si++){
	my $kmer=substr($oligo_seq,$si,$k);
	if (!exists $kmers_hash{$kmer}){
	  $kmers_hash{$kmer}={};
	}
	$kmers_hash{$kmer}->{$oligo_name}=1;
  }
  $oi++;
}


#my $low_RMT_quality_counter=0;
my $non_polyT_counter=0;
my $unknown_pool_barcode_counter=0;
my $counter=0;

my $r1_open_cmd = ($r1_file =~ /\.gz$/) ? "/bin/gzip -d -c $r1_file |" : $r1_file; 
my $r2_open_cmd = ($r2_file =~ /\.gz$/) ? "/bin/gzip -d -c $r2_file |" : $r2_file;
my $out_open_cmd = ($out_fastq_file =~ /\.gz$/) ? "| /bin/gzip -c > $out_fastq_file" : "> $out_fastq_file";
open(DATA, $r1_open_cmd) || die "cannot open file $r1_file.\n";
open(DATA2,$r2_open_cmd) || die "cannot open file $r2_file.\n";
open(OUT, $out_open_cmd) || die "cannot open file $out_fastq_file to write.\n";


while(<DATA>) {
  chomp;
  my $header1 = $_;
  my $seq1 = <DATA>;
  chomp($seq1);
  my $l31 = <DATA>;
  chomp($l31);
  my $quality1 = <DATA>;
  chomp($quality1);
  my $header2 = <DATA2>;
  chomp($header2);
  my $seq2 = <DATA2>;
  chomp($seq2);
  my $l32 = <DATA2>;
  chomp($l32);
  my $quality2 = <DATA2>;
  chomp($quality2);
  
  my @header1_split=split(":",$header1);
  my $i5=$header1_split[$#header1_split];
  my $low_RMT_quality_flag=0;
  my $low_RMT_quality_bp_counter=0;

  my $low_well_barcode_quality_flag=0;
  my $low_well_barcode_quality_bp_counter=0;
  my $non_polyT_flag=0;

  
  my $mRNA_seq="";
  my $mRNA_quality="";
  my $RMT_seq="";
  my $polyT_seq="";

  my $RMT_quality="";
  my $pool_barcode="";
  my $pool_barcode_quality="";

  my $well_barcode="";
  my $well_barcode_quality="";
  my $outs;
  my $index;
  
  for (my $ci=0; $ci< $#plate_barcode_i5_intervals;$ci+=2){
	$pool_barcode=$pool_barcode.substr($i5,$plate_barcode_i5_intervals[$ci],$plate_barcode_i5_intervals[$ci+1]);
  }

  for (my $ci=0; $ci< $#plate_barcode_intervals;$ci+=2){
	$pool_barcode=$pool_barcode.substr($seq1,$plate_barcode_intervals[$ci],$plate_barcode_intervals[$ci+1]);
	$pool_barcode_quality=$pool_barcode_quality.substr($quality1,$plate_barcode_intervals[$ci],$plate_barcode_intervals[$ci+1]);
  }
  $pool_barcode_quality=~ s/-/^/g;

  if (!( exists $pool_barcodes{$pool_barcode})){
	$unknown_pool_barcode_counter++;
	next;
  }
  my @well_barcode_intervals=@{$well_barcode_intervals_hash{$pool_barcode}};
  my @polyT_intervals=@{$polyT_intervals_hash{$pool_barcode}};
  my @RMT_intervals=@{$RMT_intervals_hash{$pool_barcode}};

  for (my $ci=0; $ci< $#well_barcode_intervals;$ci+=2){
	$well_barcode=$well_barcode.substr($seq2,$well_barcode_intervals[$ci],$well_barcode_intervals[$ci+1]);
	$well_barcode_quality=$well_barcode_quality.substr($quality2,$well_barcode_intervals[$ci],$well_barcode_intervals[$ci+1]);
  }
  $well_barcode_quality=~ s/-/^/g;

  for (my $ci=0; $ci< $#polyT_intervals;$ci+=2){
	$polyT_seq=$polyT_seq.substr($seq2,$polyT_intervals[$ci],$polyT_intervals[$ci+1]);
  }
	    
  for (my $ci=0; $ci< $#RMT_intervals;$ci+=2){
	$RMT_seq=$RMT_seq.substr($seq2,$RMT_intervals[$ci],$RMT_intervals[$ci+1]);
	$RMT_quality=$RMT_quality.substr($quality2,$RMT_intervals[$ci],$RMT_intervals[$ci+1]);
  }
  $RMT_quality=~ s/-/^/g;

  $polyT_seq=~ s/T//;

  if (length($polyT_seq)>0){
	$non_polyT_flag=1;
  }


  for (my $ci=0; $ci< $#mRNA_intervals;$ci+=2){
	$mRNA_seq=$mRNA_seq.substr($seq1,$mRNA_intervals[$ci],$mRNA_intervals[$ci+1]);
	$mRNA_quality=$mRNA_quality.substr($quality1,$mRNA_intervals[$ci],$mRNA_intervals[$ci+1]);
  }

  my %oligo_kmer_counter={};
  for (my $si=0;$si<length($mRNA_seq)-$k+1;$si++){
	my $kmer=substr($mRNA_seq,$si,$k);
	if (exists $kmers_hash{$kmer}){
	  my @cur_oligos=keys %{$kmers_hash{$kmer}};
	  for (my $oi=0;$oi<=$#cur_oligos;$oi++){
		$oligo_kmer_counter{$cur_oligos[$oi]}++;
	  }
	}
  }

  my $similar_oligo="NA";
  for (my $oi=0;$oi<=$#oligos;$oi++){
	if ((0+$oligo_kmer_counter{$oligos[$oi]})/(length($mRNA_seq)-$k+1)>.3){
	  $similar_oligo=$oligos[$oi];
	  last;
	}
  }
   my $new_header = $header1."_barcode=".$similar_oligo."-".$pool_barcode_quality."-".$well_barcode_quality."-".$RMT_quality."-".$pool_barcode."-".$well_barcode."-".$RMT_seq;
  $new_header =~ s/ /_/;
  $outs="$outs$new_header\n$mRNA_seq\n$l31\n$mRNA_quality\n";
  print OUT $outs;

  if ($non_polyT_flag==1){
	$non_polyT_counter++;
  }
  $counter++;
  #  if ($counter==10000){last;}
}
close(DATA);
close(OUT);  

open(OUT2, ">$out_stats_file") || die "cannot open file $out_stats_file to write.\n";
print OUT2 "unkown_pool_barcode\ttotal\n";
print OUT2 $unknown_pool_barcode_counter."\t".$counter."\n";
close(OUT2);
