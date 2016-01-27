use strict;

if ($#ARGV != 1) {
 die("usage: demultiplex.pl scdb_path amp_batch\n");
}

my $scdb_path=$ARGV[0];
my $amp_batch=$ARGV[1];


my $genomic_bin_size=100;
my $hamming_thresh_cell_barcode=1;

my $hamming_thresh_RMT=1;

my $qual_string="#$%&'()*+,^./0123456789:;<=>?\@ABCDEFGHIJ";

my %quality_hash={};
for (my $qi=0;$qi<length($qual_string);$qi++){
  $quality_hash{substr($qual_string,$qi,1)}=$qi+2;
}

my $seq_batch="";
my $wells_cells_fn = $scdb_path."/annotations/wells_cells.txt";
my $seq_batches_fn = $scdb_path."/annotations/seq_batches.txt";
my $amp_batches_fn = $scdb_path."/annotations/amp_batches.txt";
my $config_fn = $scdb_path."/config/config.txt";

my $umitab_fn = $scdb_path."/output/umi.tab/".$amp_batch.".txt";
my $offsetab_fn = $scdb_path."/output/offset.tab/".$amp_batch.".txt";
my $singleofftab_fn = $scdb_path."/output/singleton_offset.tab/".$amp_batch.".txt";

if (! -d $scdb_path."/output/umi.tab/" ) {
  mkdir $scdb_path."/output/umi.tab/";
}
if (! -d $scdb_path."/output/offset.tab/" ) {
  mkdir $scdb_path."/output/offset.tab/";
}
if (! -d $scdb_path."/output/singleton_offset.tab/" ) {
  mkdir $scdb_path."/output/singleton_offset.tab/";
}
if (! -d $scdb_path."/output/QC/read_stats" ) {
  mkdir $scdb_path."/output/QC/read_stats";
}
if (! -d $scdb_path."/output/QC/read_stats_amp_batch" ) {
  mkdir $scdb_path."/output/QC/read_stats_amp_batch";
}
if (! -d $scdb_path."/output/QC/rmt_stats" ) {
  mkdir $scdb_path."/output/QC/rmt_stats";
}
if (! -d $scdb_path."/output/QC/noffsets_per_rmt_distrib" ) {
  mkdir $scdb_path."/output/QC/noffsets_per_rmt_distrib";
}
if (! -d $scdb_path."/output/QC/nreads_per_rmt_distrib" ) {
  mkdir $scdb_path."/output/QC/nreads_per_rmt_distrib";
}
if (! -d $scdb_path."/output/QC/rmt_nuc_per_pos" ) {
  mkdir $scdb_path."/output/QC/rmt_nuc_per_pos";
}



my $debug_offsets_fn = $scdb_path."/_debug/".$amp_batch."/offsets.txt";
#my $debug_offsets2_fn = $scdb_path."/processed_data/".$amp_batch."/offsets2.txt";
my $debug_rmts_fn = $scdb_path."/_debug/".$amp_batch."/RMTs.txt";

my $read_stats_fn = $scdb_path."/output/QC/read_stats/".$amp_batch.".txt";
my $read_stats_amp_batch_fn = $scdb_path."/output/QC/read_stats_amp_batch/".$amp_batch.".txt";
my $rmt_stats_fn = $scdb_path."/output/QC/rmt_stats/".$amp_batch.".txt";
my $rmt_nuc_per_pos_fn = $scdb_path."/output/QC/rmt_nuc_per_pos/".$amp_batch.".txt";
my $noffsets_per_rmt_fn = $scdb_path."/output/QC/noffsets_per_rmt_distrib/".$amp_batch.".txt";
my $nreads_per_rmt_fn = $scdb_path."/output/QC/nreads_per_rmt_distrib/".$amp_batch.".txt";

###########################################################################
#Read config file
my %config_hash;
open(CONFIG_FILE, $config_fn) || die "ERROR: cannot open file $config_fn to read.\n";
while(<CONFIG_FILE>) {
  chomp;
  my @row= split("=");
  $config_hash{$row[0]}=$row[1];
  print $row[0]."\t".$row[1]."\n";
}
my $well_barcode_min_quality_thresh=27;
if (exists $config_hash{"well_barcode_min_quality_thresh"}){
  $well_barcode_min_quality_thresh=$config_hash{"well_barcode_min_quality_thresh"};
}
my $well_barcode_max_num_of_bad_quality_bp=1;
if (exists $config_hash{"well_barcode_max_num_of_bad_quality_bp"}){
  $well_barcode_max_num_of_bad_quality_bp=$config_hash{"well_barcode_max_num_of_bad_quality_bp"};
}
my $pool_barcode_min_quality_thresh=27;
if (exists $config_hash{"pool_barcode_min_quality_thresh"}){
  $pool_barcode_min_quality_thresh=$config_hash{"pool_barcode_min_quality_thresh"};
}
my $pool_barcode_max_num_of_bad_quality_bp=1;
if (exists $config_hash{"pool_barcode_max_num_of_bad_quality_bp"}){
  $pool_barcode_max_num_of_bad_quality_bp=$config_hash{"pool_barcode_max_num_of_bad_quality_bp"};
}

my $fdr=0.25;
if (exists $config_hash{"fdr_offset_err"}){
  $fdr=$config_hash{"fdr_offset_err"};
}

my $gene_intervals_fn=$scdb_path."/".$config_hash{"gene_intervals_file"};
my $spike_seq_fn=$scdb_path."/".$config_hash{"spike_seq_file"};
my $oligos_fn=$scdb_path."/".$config_hash{"oligos_file"};
my @oligos=();
open(OLIGO_FILE, $oligos_fn) || die "ERROR: cannot open file $oligos_fn to read.\n";
while(<OLIGO_FILE>) {
  chomp;
  my @line= split("=");
  push(@oligos,$line[0]);
}

print "fdr=$fdr\n";
print "well_barcode_min_quality_thresh=$well_barcode_min_quality_thresh\n";



###################################################################

sub map_to_gene{
  
  my $chr=$_[0];
  my $coor=$_[1];
  my $strand=$_[2];
  my $gene_hash_ref=$_[3];
  my $bin=int($coor/$genomic_bin_size);
  my $key=$chr."_".$strand."_".$bin;
 
  if (!exists $gene_hash_ref->{$key}){
	return ("");
  }
  else{
	my $gene=$gene_hash_ref->{$key};
	return ($gene);
  }
}

##################################################################
sub median {
    my @a = sort {$a <=> $b} @_;
    my $length = scalar @a;
    return undef unless $length;
    ($length % 2)
        ? $a[$length/2]
        : ($a[$length/2] + $a[$length/2-1]) / 2.0;
}

##################################################################
sub fdr_thresh {

  my @pvalues =@_;

  my @sorted_pvalues = sort { $a <=> $b } @pvalues;
  my $n=($#sorted_pvalues+1);
  for (my $j=($#sorted_pvalues+1);$j>=1;$j--) {
#	print $j."/".$n."\t".$sorted_pvalues[$j-1]."\t".($thresh*$j/$n)."\n";

	if ($sorted_pvalues[$j-1]<=($fdr*$j/$n)){
	  return($sorted_pvalues[$j-1]);
	}
  }
  return(0);
}


sub logfact {
   return gammln(shift(@_) + 1.0);
}

sub hypergeom {
   # There are m "bad" and n "good" balls in an urn.
   # Pick N of them. The probability of i or more successful selections:
   # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
   my ($n, $m, $N, $i) = @_;
#   if ($m==0||$n==0){
#	 return(1);
#   }
  
   my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
   my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
   return exp($loghyp1 - $loghyp2);
}

sub gammln {
  my $xx = shift;
  my @cof = (76.18009172947146, -86.50532032941677,
             24.01409824083091, -1.231739572450155,
             0.12086509738661e-2, -0.5395239384953e-5);
  my $y = my $x = $xx;
  my $tmp = $x + 5.5;
  $tmp -= ($x + .5) * log($tmp);
  my $ser = 1.000000000190015;
  for my $j (0..5) {
     $ser += $cof[$j]/++$y;
  }
  -$tmp + log(2.5066282746310005*$ser/$x);
}


sub hamming {
	my($s1, $s2,$max) = @_;

	my($dist) = 0;
	for(my($i) = 0; $i < length($s1); $i++) {
		if(substr($s1, $i, 1) ne substr($s2, $i, 1)) {
			$dist++;
			if ($dist>=$max){
			  return($dist)
			}
		  }
	  }
	return($dist);
  }


############################################################################################################
#
#  filter_errors
#
#  This function makes filtering decisions for RMT and barcode errors and return p-values for offset errors.
#
#  Parameter:
#  $hash_in -       ref for hash that stores all the input information of a gene
#                   keys: tab delimited (offset, well ID, barcode sequence, RMT sequence)
#                   values: number of reads
#  $hash_out -      ref for hash that  stores all the input+output information of a gene
#                   keys: tab delimited (offset, well ID, barcode sequence, RMT sequence, ...)
#  $cur_gene -      The name of the gene
#  $offsets_p1 -    ref for hash that stores the hypergeometic p-values of the "lonely offset" test.
#                   keys: tab delimited (gene,offset)
#                   values: raw p-value
#  $offsets_p2 -    ref for hash that stores the hypergeometic p-values of the "read poor offset" test.
#                   keys: tab delimited (gene,offset)
#                   values: raw p-value

sub filter_errors {
  my($hash_in, $hash_out,$RMT_hash,$cur_gene,$offsets_p1,$offsets_p2,$hamming_thresh_RMT,$hamming_thresh_cell_barcode) = @_;
  my %filt_RMT_hash=();
  my %filt_barcode_hash=();
  my @data=keys %{$hash_in};
  my(%RMT_off_map);
  my(%RMT_off_list);
  my(%RMT_count);
 

  my(%barcode_off_map);
  my(%barcode_off_list);
  my(%barcode_count);
  my(%barcode_to_wellid);

  my(%cell_offset_to_umis);
  my(%offset_to_mol);
  my(%mol_to_n_offsets);
  my(%offset_stats_hash);

  
  my($n_total_offsets_with_single_read)=0;
  my($n_total_offsets_with_multiple_reads)=0;
  my($tot);
 
  for(my($i) = 0; $i <= $#data; $i++) {
	my($offset,$wellid,$barcode, $RMT) = split("\t", $data[$i]);
	$barcode_to_wellid{$barcode}=$wellid;
	if(!exists($RMT_off_map{$wellid})){
	  $RMT_off_map{$wellid}={};
	}
	$RMT_off_map{$wellid}->{"$RMT\t$offset"} = 1;
	if(!exists($RMT_off_list{$wellid}->{$RMT})) {
	  $RMT_off_list{$wellid}->{$RMT} = ();
	}
	push(@{$RMT_off_list{$wellid}->{$RMT}}, $offset);
	if (!exists $RMT_count{$wellid}){
	  $RMT_count{$wellid}={};
	}
	$RMT_count{$wellid}->{$RMT}++;
 
	if(!exists($barcode_off_map{$RMT})){
	  $barcode_off_map{$RMT}={};
	}
	$barcode_off_map{$RMT}->{"$barcode\t$offset"} = 1;
	if(!exists($barcode_off_list{$RMT}->{$barcode})) {
	  $barcode_off_list{$RMT}->{$barcode} = ();
	}
	push(@{$barcode_off_list{$RMT}->{$barcode}}, $offset);
	if (!exists $barcode_count{$RMT}){
	  $barcode_count{$RMT}={};
	}
	$barcode_count{$RMT}->{$barcode}++;

	my $reads=%$hash_in->{$data[$i]};
	if(!exists($cell_offset_to_umis{$wellid})) {
		$cell_offset_to_umis{$wellid}= {};
	  }
	  if(!exists($cell_offset_to_umis{$wellid}->{$offset})) {
		$cell_offset_to_umis{$wellid}->{$offset}= {};
	  }
	$cell_offset_to_umis{$wellid}->{$offset}->{$RMT}=1;
# Need to decide whether to use stats from empty wells
#	if (exists(%$single_cells->{$wellid})){
	  
	  if(!exists( $offset_to_mol{$offset})) {
		$offset_to_mol{$offset} = {};
	  }
	  $offset_to_mol{$offset} ->{$wellid."\t".$RMT}=$reads;
	  if ($reads==1){
		$n_total_offsets_with_single_read++;
	  }
	  if ($reads>1){
		$n_total_offsets_with_multiple_reads++;
		}
	  $mol_to_n_offsets{$wellid."\t".$RMT}++;
#	}
	$tot++;
  }
  
  my @all_wellids=keys %RMT_count;
  for (my $i=0;$i<=$#all_wellids;$i++){
	my $wellid=$all_wellids[$i];
	my(@all_RMTs) = keys %{$RMT_count{$wellid}};
	#	print STDERR "Total #RMTs=".$#all_RMTs."\n";
	my(@RMTs) = sort { $RMT_count{$wellid}->{$a} <=> $RMT_count{$wellid}->{$b} } @all_RMTs;
	for(my($i) = 0; $i <= $#RMTs; $i++) {
	  my($RMT) = $RMTs[$i];
	  my($i_offs) = $RMT_off_list{$wellid}->{$RMT};
	  my($filt) = 0;
	  for(my($j) = $i + 1; $j <=$#RMTs;$j++) {
	    my($dom_RMT) = $RMTs[$j];
		if(hamming($RMT, $dom_RMT,$hamming_thresh_RMT+1) > $hamming_thresh_RMT) {
	  	  next;
		}
		my($off);
		my($indep) = 0;
		foreach $off (@{$i_offs}) {
		  if(!exists($RMT_off_map{$wellid}->{"$dom_RMT\t$off"})) {
			$indep = 1;
			last;
		  }
		}
		if($indep == 0) {
		  $filt = 1;
		}
	  }
	  $filt_RMT_hash{$wellid."_".$RMT} = $filt;
	}
  }
	
  my @all_RMTs=keys %barcode_count;
  for (my $i=0;$i<=$#all_RMTs;$i++){
	my $RMT=$all_RMTs[$i];
	my(@all_barcodes) = keys %{$barcode_count{$RMT}};
	my(@barcodes) = sort { $barcode_count{$RMT}->{$a} <=> $barcode_count{$RMT}->{$b} } @all_barcodes;
	for(my($i) = 0; $i <= $#barcodes; $i++) {
	  my($barcode) = $barcodes[$i];
	  my($i_offs) = $barcode_off_list{$RMT}->{$barcode};
	  my($filt) = 0;
	  for(my($j) = $i + 1; $j <= $#barcodes; $j++) {
		my($dom_barcode) = $barcodes[$j];
		  if(hamming($barcode, $dom_barcode,$hamming_thresh_cell_barcode+1) >  $hamming_thresh_cell_barcode | $barcode_to_wellid{$barcode} eq $barcode_to_wellid{$dom_barcode}) {
			next;
		  }
		my($off);
		my($indep) = 0;
		foreach $off (@{$i_offs}) {
		  if(!exists($barcode_off_map{$RMT}->{"$dom_barcode\t$off"})) {
			$indep = 1;
			last;
		  }
		}
		if($indep == 0) {
		  $filt = 1;
		}
	  }
	  $filt_barcode_hash{$barcode."_".$RMT} = $filt;
	
	}	 

  }

   
  
  my($n_total_lonely_offsets)=0;
  my($n_total_offsets_with_friends)=0;
  
  foreach my $mol (keys %mol_to_n_offsets) {
	if ($mol_to_n_offsets{$mol}==1){
	  $n_total_lonely_offsets++;
	}
	else{
	  $n_total_offsets_with_friends+=$mol_to_n_offsets{$mol};
	}
  }

  foreach my $offset (keys %offset_to_mol) {
	
	my($n_current_lonely_offsets)=0;
	my($n_current_offsets_with_friends)=0;
	
	my($n_current_offsets_with_single_read)=0;
	my($n_current_offsets_with_multiple_read)=0;
	my (@v_n_reads_per_offset,@v_n_offsets_per_molecules_participating_in_this_offset,$i);
	
	foreach my $mol (keys %{$offset_to_mol{$offset}}) {
	  $v_n_reads_per_offset[$i]=%{$offset_to_mol{$offset}}->{$mol};
	  $v_n_offsets_per_molecules_participating_in_this_offset[$i]=$mol_to_n_offsets{$mol};
	  
	  if ($mol_to_n_offsets{$mol}==1){
		$n_current_lonely_offsets++;
	  }
	  else {
		$n_current_offsets_with_friends++;
	  }

	  if (%{$offset_to_mol{$offset}}->{$mol}==1){
		$n_current_offsets_with_single_read++;
	  }
	  else {
		$n_current_offsets_with_multiple_read++;
	  }
	  
	  
	 $i++;
	}
  
	my($n_molecules_with_this_offset) =scalar keys %{$offset_to_mol{$offset}};
	my $median_n_reads_per_offset=median(@v_n_reads_per_offset);
	my $median_n_offsets=median(@v_n_offsets_per_molecules_participating_in_this_offset);
	my $ratio=$n_molecules_with_this_offset/$median_n_offsets;
	my $hypergeom_p1=hypergeom($n_total_lonely_offsets,$n_total_offsets_with_friends,$n_current_lonely_offsets+$n_current_offsets_with_friends,$n_current_lonely_offsets);	
	my $hypergeom_p2=hypergeom($n_total_offsets_with_single_read,$n_total_offsets_with_multiple_reads,$n_current_offsets_with_single_read+$n_current_offsets_with_multiple_read,$n_current_offsets_with_single_read);

	my $stats="$n_total_lonely_offsets\t$n_total_offsets_with_friends\t$n_current_lonely_offsets\t$n_current_offsets_with_friends\t$n_total_offsets_with_single_read\t$n_total_offsets_with_multiple_reads\t$n_current_offsets_with_single_read\t$n_current_offsets_with_multiple_read";
#	print $cur_gene."\t".$offset."\t".$stats."\n";
	
	$offset_stats_hash{$cur_gene."\t".$offset}=$stats;
	%$offsets_p1->{$cur_gene."\t".$offset}=$hypergeom_p1;
	%$offsets_p2->{$cur_gene."\t".$offset}=$hypergeom_p2;
  }
 

  foreach my $s (keys %{$hash_in}){
	my($offset,$wellid,$barcode, $RMT) = split("\t", $s);

	my $s2=$s."\t".%$hash_in->{$s}."\t".$filt_RMT_hash{$wellid."_".$RMT}."\t".$filt_barcode_hash{$barcode."_".$RMT}."\t".%$offsets_p1->{$cur_gene."\t".$offset}."\t".%$offsets_p2->{$cur_gene."\t".$offset}."\t".$offset_stats_hash{$cur_gene."\t".$offset};
	
	$hash_out->{$s2}=$hash_in->{$s};
  }
}

#############################################################################



my $index;

my %hash1;
my %hash2;

my %number_of_cells;

my %RMT_hash;
my %rmt_stats;
my %nuc_per_pos_stats;
my @well_list=();


# read sample_index
# make barcode hash
# read amplification batch, map reads to genes, store in hashes for filtering
# write debug files
# write cell gene matrix.


my %cell_barcode_to_well_id_hash=();
my %extended_cell_barcode_to_well_id_hash=();
my %extended_cell_barcode_in_amp_batch_hash=();
my %extended_pool_barcode_hash=();
my %column_num;
my %gene_hash;
my %read_counter;
my %oligo_counter;



###########################################################################
#Read gene intervals
my %binned_coordinate_to_gene_hash;
open(GENE_INTERVALS_FILE, $gene_intervals_fn) || die "ERROR: cannot open file $gene_intervals_fn to read.\n";
my $line=<GENE_INTERVALS_FILE>;
chomp $line;
my @column_names=split("\t",$line);
for (my $i=0;$i<=$#column_names;$i++){
  $column_num{$column_names[$i]}=$i;
}
while(<GENE_INTERVALS_FILE>) {
  chomp;
  my @row= split("\t");
  my $chr=$row[$column_num{"chrom"}];
  my $start=$row[$column_num{"start"}];
  my $end=$row[$column_num{"end"}];
  my $strand=$row[$column_num{"strand"}];
  my $gene_name=$row[$column_num{"gene_name"}];
  $gene_hash{$gene_name}=1;
  for (my $bin=int($start/$genomic_bin_size);$bin<=int($end/$genomic_bin_size);$bin++){  
	my $key=$chr."_".$strand."_".$bin;
	$binned_coordinate_to_gene_hash{$key}=$gene_name;
  }
}

my %spike_ins;
open(SPIKE_FILE, $spike_seq_fn) || die "ERROR: cannot open file $spike_seq_fn to read.\n";
<SPIKE_FILE>;
while(<SPIKE_FILE>) {
  my @l=split("\t",$_);
  $spike_ins{$l[0]}=1;
}

my @genes=sort keys %gene_hash;
for (my $j=0; $j<=$#genes;$j++){
  $hash1{$genes[$j]}={};
  $hash2{$genes[$j]}={};
}

###########################################################################
#Read amplification batch details
open(AMP_BATCHES_FILE, $amp_batches_fn) || die "ERROR: cannot open file $amp_batches_fn to read.\n";
$line=<AMP_BATCHES_FILE>;
chomp $line;
@column_names=split("\t",$line);
for (my $i=0;$i<=$#column_names;$i++){
  $column_num{$column_names[$i]}=$i;
}
my $pool_barcode;
my $protocol_version;
my $seq_batch;

my $cur_amp_batch="";
while (<AMP_BATCHES_FILE>){
  chomp;
  my @row= split("\t");
  $cur_amp_batch=$row[$column_num{"Amp_batch_ID"}];
  if ($cur_amp_batch eq $amp_batch){
	$pool_barcode=$row[$column_num{"Pool_barcode"}];
	$protocol_version=$row[$column_num{"Protocol_version_ID"}];
	$seq_batch=$row[$column_num{"Seq_batch_ID"}];
	last;
  }
}
if ($cur_amp_batch eq ""){
  die("ERROR: $amp_batch does not exist!");
}

###########################################################################
#Read barcodes

open(INDEX_FILE, $wells_cells_fn) || die "ERROR: cannot open file $wells_cells_fn to read.\n";

%column_num=();
my @column_names=split("\t",<INDEX_FILE>);
for (my $i=0;$i<=$#column_names;$i++){
  $column_num{$column_names[$i]}=$i;
}

while(<INDEX_FILE>) {
  chomp;
  my @row= split("\t");
  my $cur_amp_batch=$row[$column_num{"Amp_batch_ID"}];
#  print $cur_amp_batch."\t".$amp_batch."\t".$column_num{"Amp_batch_ID"}."\n";
  if ($cur_amp_batch eq $amp_batch){
	my $sample_index=$row[$column_num{"Well_ID"}];
	my $sample_barcode=$row[$column_num{"Cell_barcode"}];
	$number_of_cells{$sample_index}=$row[$column_num{"Number_of_cells"}];
	push(@well_list,$sample_index);

	$extended_cell_barcode_in_amp_batch_hash{$sample_barcode}=1;
	$cell_barcode_to_well_id_hash{$sample_barcode}=$sample_index;
	$extended_cell_barcode_to_well_id_hash{$sample_barcode}=$sample_index;
	$read_counter{$sample_index}={};
	$oligo_counter{$sample_index}={};
  }
}



my @nucs=("A","C","G","T");
foreach my $barcode (keys %cell_barcode_to_well_id_hash){
  for (my $i=0;$i<length($barcode);$i++){
	if ((substr($barcode,$i,1) ne "_") && (substr($barcode,$i,1) ne "N")){
	  for (my $to_nuci=0;$to_nuci<4;$to_nuci++){
		my $barcode2=$barcode;
		substr($barcode2,$i,1)=$nucs[$to_nuci];
		if (exists $cell_barcode_to_well_id_hash{$barcode2}){
		  next;
		}
		if (!exists $extended_cell_barcode_to_well_id_hash{$barcode2}){
		  $extended_cell_barcode_to_well_id_hash{$barcode2}=$cell_barcode_to_well_id_hash{$barcode};
		  $extended_cell_barcode_in_amp_batch_hash{$barcode2}= $extended_cell_barcode_in_amp_batch_hash{$barcode};
		}
		else{
		  $extended_cell_barcode_to_well_id_hash{$barcode2}=-1;
		  $extended_cell_barcode_in_amp_batch_hash{$barcode2}=-1;
		}
#		print $barcode."\t".$barcode2."\t".$cell_barcode_to_well_id_hash{$barcode}."\t".$extended_cell_barcode_to_well_id_hash{$barcode}."\t".$extended_cell_barcode_in_amp_batch_hash{$barcode2}."\n";
	  }
	}
  }
}


my @nucs=("A","C","G","T");
for (my $i=0;$i<length($pool_barcode);$i++){
  if ((substr($pool_barcode,$i,1) ne "_") && (substr($pool_barcode,$i,1) ne "N")){
	for (my $to_nuci=0;$to_nuci<4;$to_nuci++){
	  my $barcode2=$pool_barcode;
	  substr($barcode2,$i,1)=$nucs[$to_nuci];
	  $extended_pool_barcode_hash{$barcode2}=1;
	}
  }
}


#############################################################
# Read Reads
#

my @fn = <$scdb_path/_trimmed_mapped_reads/$seq_batch/*.sam>;

if ($#fn<=0) {
 die("ERROR: $scdb_path/_trimmed_mapped_reads/$seq_batch/ dir doesn't contain sam files\n");
}
print "reading $#fn sam files\n";



my %rmt_status;

my %reads_per_rmt;
my %offsets_per_rmt;



# Reading mapped reads and filtering reads that are not mapped to gene intervals
my @parsed_line;
my @parsed_header;
my $counter_reads_with_unknown_barcodes=0;
my $counter_bad_quality_well_barcodes=0;
my $counter_bad_quality_pool_barcodes=0;
my $counter_bad_quality_RMT=0;
my $counter_amp_batch_reads=0;
my $counter_seq_batch_reads=0;

my $fh=select(STDOUT);
$|=1;
select($fh);

foreach my $fn (@fn) {
  open(DATA, $fn) || die "ERROR: cannot open file $fn to read.\n";
  print "Reading $fn\n";
 # my $ii=0;
  while(<DATA>) {
	chomp;
	my $line=$_;
	if ($line =~ /^@/){
	   next;
	 }
#	$ii++;   
#	if ($ii==1e4){
#	  last;
#	}

	@parsed_line=split("\t",$line);
	my $barcode_info=$parsed_line[0];
	$barcode_info =~ s/^.*barcode=//;
	@parsed_header=split("-",$barcode_info); 
	
	my $cur_pool_barcode=$parsed_header[$#parsed_header-2];

	$counter_seq_batch_reads++;
	if ($extended_pool_barcode_hash{$cur_pool_barcode}==0){
	  next;
	}

	my $cur_pool_quality=$parsed_header[$#parsed_header-5];
	my $n_bad_poolbarcode_quality_bps=0;
	my @cur_poolbarcode_arr=split(//,$cur_pool_quality);

	for (my $qi=0;$qi<$#cur_poolbarcode_arr;$qi++){
	  my $bp_qual=(0+$quality_hash{$cur_poolbarcode_arr[$qi]});
	  if ( ($bp_qual)<$pool_barcode_min_quality_thresh){
		$n_bad_poolbarcode_quality_bps++;
	  }
	}

	if ($n_bad_poolbarcode_quality_bps>$pool_barcode_max_num_of_bad_quality_bp){
	  $counter_bad_quality_pool_barcodes++;
	  next;
	}

	$counter_amp_batch_reads++;

	my $cur_wellbarcode_quality=$parsed_header[$#parsed_header-4];
	my $cur_RMT_quality=$parsed_header[$#parsed_header-3];

    my $cell_barcode=$parsed_header[$#parsed_header-1];
	my $RMT= $parsed_header[$#parsed_header];
    my $flag=$parsed_line[1];
	my $mapq=$parsed_line[4];

	my $strand=1;  
	my $well_id;
#	my $well_barcode_min_bp_quality=1000;
	my $n_bad_quality_bps=0;
	my @cur_wellbarcode_arr=split(//,$cur_wellbarcode_quality);

	for (my $qi=0;$qi<$#cur_wellbarcode_arr;$qi++){
	  my $bp_qual=(0+$quality_hash{$cur_wellbarcode_arr[$qi]});
	  if ( ($bp_qual)<$well_barcode_min_quality_thresh){
		$n_bad_quality_bps++;
	  }
	}

	if ($n_bad_quality_bps>$well_barcode_max_num_of_bad_quality_bp){
	  $counter_bad_quality_well_barcodes++;
	  next;
	}

	if ($RMT =~ m/N/){
 	  $counter_bad_quality_RMT++;
	  next;
	}

    if ($extended_cell_barcode_in_amp_batch_hash{$cell_barcode}<1){
	  $counter_reads_with_unknown_barcodes++;
	  next;
    }

	$well_id= $extended_cell_barcode_to_well_id_hash{$cell_barcode};
	if ($well_id<0){
	  next;
	}

# At this point we know that this read has a valid pool and well barcode
	$read_counter{$well_id}->{"total"}++;
	my $oligo_est=$parsed_header[$#parsed_header-6];
	if ($flag&4){
	  if ($oligo_est eq "NA"){
		$read_counter{$well_id}->{"unmapped"}++;
	  }
	  else {
		$oligo_counter{$well_id}->{$oligo_est}++;
	  }
      next;
    }
## checking multi-mapping
	my $as=-1000;
	my $xs=-1000;
	for (my $j=9;$j<=$#parsed_line;$j++){

	  if ($parsed_line[$j]=~/^AS:i:/){
		($as=$parsed_line[$j]) =~ s/^AS:i://;
	  }
	  if ($parsed_line[$j]=~/^XS:i:/){
		($xs=$parsed_line[$j]) =~ s/^XS:i://;
	  }
	}
	if (($xs>=$as)&($as>-1000)){
	  $read_counter{$well_id}->{"non_unique_mapping"}++;
	  next;
	}
	
	if ($mapq<30){
	  $read_counter{$well_id}->{"lowmapq"}++;
	  next;
	}
  
	if ($flag&16){
      $strand=-1;
    }

    my $chr=$parsed_line[2];
    my $coor=$parsed_line[3];
    my $gene=map_to_gene($chr,$coor,$strand,\%binned_coordinate_to_gene_hash);
	
	if ($gene eq ""){
	  $read_counter{$well_id}->{"mapped_to_nongenic"}++;
	  next;
	}

	my $pos=$chr."_".$strand."_".$coor;
	if (!exists $RMT_hash{$RMT}){
	  $RMT_hash{$RMT}={};
	}
	$RMT_hash{$RMT}->{"$well_id\t$gene\t$pos\t$chr\t$strand\t$coor"}=1;
	# incrementing read counts;
	$hash1{$gene}->{"$pos\t$well_id\t$cell_barcode\t$RMT"}++;
	
  }
}

print "Finished reading sam files\n";


##############################################################################
#open(TMP,">$debug_offsets2_fn" )|| die "could not open debug file\n";
#my @RMTS=keys %RMT_hash;
#for (my $ii=0;$ii<=$#RMTS;$ii++){
#  my @offs= keys %{$RMT_hash{$RMTS[$ii]}};
#  for (my $jj=0;$jj<=$#offs;$jj++){
#	print TMP $RMTS[$ii]."\t".$offs[$jj]."\n";
#  }
#}
#close(TMP);
##############################################################################

# Now we will filter reads with RMT or barcode single bp error.
print "Filtering single bp erros\n";

open(OUT_DEBUG_OFFSETS, ">$debug_offsets_fn") || die "ERROR: cannot open file $debug_offsets_fn to read.\n";
print OUT_DEBUG_OFFSETS "gene_name\toffset\twell_id\tbarcode\tRMT\treads\tfiltered_RMT_err\tfiltered_barcode_err\tp_lonely_offset\tp_readpoor_offset\tn_total_lonely_offsets\tn_total_offsets_with_friends\tn_current_lonely_offsets\tn_current_offsets_with_friends\tn_total_offsets_with_single_read\tn_total_offsets_with_multiple_reads\tn_current_offsets_with_single_read\tn_current_offsets_with_multiple_read\n";

my @genes=sort keys %gene_hash;
my %offsets_p1;
my %offsets_p2;
for (my $j=0; $j<=$#genes;$j++){
  my $nrows=keys %{$hash1{$genes[$j]}};
  #print $genes[$j]."\t$j\t$nrows\n";
  filter_errors(\%{$hash1{$genes[$j]}},\%{$hash2{$genes[$j]}},\%RMT_hash,$genes[$j],\%offsets_p1,\%offsets_p2,$hamming_thresh_RMT,$hamming_thresh_cell_barcode);
  my @data=keys %{$hash2{$genes[$j]}};
  for (my $i=0;$i<=$#data;$i++){
	print OUT_DEBUG_OFFSETS $genes[$j]."\t".$data[$i]."\n";
  }
}

close(OUT_DEBUG_OFFSETS);

my @p1_vec=values %offsets_p1;
my @p2_vec=values %offsets_p2;

my $fdr_thresh1=fdr_thresh(@p1_vec);
my $fdr_thresh2=fdr_thresh(@p2_vec);
print "FDR_THRESHES:\t$fdr_thresh1\t$fdr_thresh2\n";

##############################################################################

# Writing umitab, offsetab, singleofftab
print "Writing output files\n";
open(UMITAB_FILE, ">$umitab_fn") || die "ERROR: cannot open file $umitab_fn to read.\n";
open(OFFSETAB_FILE, ">$offsetab_fn") || die "ERROR: cannot open file $offsetab_fn to read.\n";
open(SINGLEOFFTAB_FILE, ">$singleofftab_fn") || die "ERROR: cannot open file $singleofftab_fn to read.\n";
print UMITAB_FILE $well_list[0];
print OFFSETAB_FILE $well_list[0];
print SINGLEOFFTAB_FILE $well_list[0];
for (my $i=1;$i<=$#well_list;$i++){
  print UMITAB_FILE "\t".$well_list[$i];
  print OFFSETAB_FILE "\t".$well_list[$i];
 print SINGLEOFFTAB_FILE "\t".$well_list[$i];
}
print UMITAB_FILE "\n";
print OFFSETAB_FILE "\n";
print SINGLEOFFTAB_FILE "\n";




print $#genes." genes\n";
for (my $j=0; $j<=$#genes;$j++){
  

  my $gene=$genes[$j];
  my $spike_or_gene="gene";
  if (exists $spike_ins{$gene}){
	$spike_or_gene="spike";
  }

  
  my @data=keys %{$hash2{$gene}};

  my %well_to_RMTs={};
  my %well_to_RMTs_ok={};

  my %wellgene_to_noffsets={};
  my %wellgene_to_nsingleton_offsets={};
  for (my $i=0;$i<=$#data;$i++){
# Iterating over all IVT products of the gene
	my($off,$wellid,$barcode, $RMT,$nreads,$filt1,$filt2,$p1,$p2) = split("\t", $data[$i]);
	my $well_gene=$wellid."\t".$gene;
	my $well_gene_rmt=$wellid."\t".$gene."\t".$RMT;
	if (!exists $well_to_RMTs{$wellid}){
	   $well_to_RMTs{$wellid}={};
	}
	if (!exists $well_to_RMTs_ok{$wellid}){
	   $well_to_RMTs_ok{$wellid}={};
	}
	$well_to_RMTs{$wellid}->{$RMT}=1;

	# Offset filtering: filtering RMT if all its IVT products are filtered	
	if ($p1<$fdr_thresh1){
	  $read_counter{$wellid}->{$spike_or_gene."\tfilt3"}+=$nreads;
	  if (!exists $rmt_status{$well_gene_rmt}){
		$rmt_status{$well_gene_rmt}=-3;
	  }
	}
	elsif ($p2<$fdr_thresh2){
	  $read_counter{$wellid}->{$spike_or_gene."\tfilt4"}+=$nreads;
	  if (!exists $rmt_status{$well_gene_rmt}){
		$rmt_status{$well_gene_rmt}=-4;
	  }
	}
	# RMT filtering
	elsif ($filt1==1){
	  $read_counter{$wellid}->{$spike_or_gene."\tfilt1"}+=$nreads;
	  if ($rmt_status{$well_gene_rmt}>=0){
		$rmt_status{$well_gene_rmt}=-1;
	  }
	}
	# Well barcode filtering
	elsif ($filt2==1){
	  $read_counter{$wellid}->{$spike_or_gene."\tfilt2"}+=$nreads;
	  if ($rmt_status{$well_gene_rmt}>=0){
		$rmt_status{$well_gene_rmt}=-2;
	  }
	}
	else{
	  ##This (gene,well,RMT,OFFSET) is not going to be filtered
	  if ($spike_or_gene eq "spike"){
		$read_counter{$wellid}->{"spike_mapped"}+=$nreads;
	  }
	  else{
		$read_counter{$wellid}->{"gene_mapped"}+=$nreads;
	  }
	  
	  $well_to_RMTs_ok{$wellid}->{$RMT}=1;
	  $wellgene_to_noffsets{$well_gene}++;
	  $rmt_status{$well_gene_rmt}=1;
	  my @rmtarr=split("",$RMT);
	  for (my $ni=0;$ni<=$#rmtarr;$ni++){
		$nuc_per_pos_stats{$ni."\t".$rmtarr[$ni]}++;
	  }
	  
	  $reads_per_rmt{$well_gene_rmt}+=$nreads;
	  $offsets_per_rmt{$well_gene_rmt}++;
	  
	  if ($offsets_per_rmt{$well_gene_rmt}==1) {
		$wellgene_to_nsingleton_offsets{$well_gene}++;
	  }elsif ($offsets_per_rmt{$well_gene_rmt}==2){
		$wellgene_to_nsingleton_offsets{$well_gene}--;
	  }
	}
  }

  print UMITAB_FILE $gene;
  print OFFSETAB_FILE $gene;
  print SINGLEOFFTAB_FILE $gene;
  for (my $i=0;$i<=$#well_list;$i++){
	my $wellid=$well_list[$i];
	my @RMTs=keys %{$well_to_RMTs{$wellid}};
	my @RMTs_ok=keys %{$well_to_RMTs_ok{$wellid}};
	
	foreach my $RMT (@RMTs){
	  $rmt_stats{$wellid."\t".$spike_or_gene."\t".$rmt_status{$wellid."\t".$gene."\t".$RMT}}++;
	}
	
	my $number_of_RMTs=1+$#RMTs_ok;
	my $number_of_offsets=0+$wellgene_to_noffsets{$well_list[$i]."\t".$gene};
	my $number_of_singleton_offsets=0+$wellgene_to_nsingleton_offsets{$well_list[$i]."\t".$gene};
	
	print UMITAB_FILE "\t".$number_of_RMTs;
	print OFFSETAB_FILE "\t".$number_of_offsets;
	print SINGLEOFFTAB_FILE "\t".$number_of_singleton_offsets;
  }
  print UMITAB_FILE "\n"; 
  print OFFSETAB_FILE "\n";
  print SINGLEOFFTAB_FILE "\n";
}
close(UMITAB_FILE);
close(OFFSETAB_FILE);
close(SINGLEOFFTAB_FILE);

# Writing stats


open(OUT_READ_STATS_AMP_BATCH, ">$read_stats_amp_batch_fn") || die "ERROR: cannot open file $read_stats_amp_batch_fn to read.\n";
print OUT_READ_STATS_AMP_BATCH "nreads_with_unknown_cell_barcode\tnreads_bad_quality_well_barcode\tnreads_bad_quality_pool_barcode\tnreads_bad_quality_RMT\tnreads_amp_batch\tnreads_seq_batch\n";

print OUT_READ_STATS_AMP_BATCH "$counter_reads_with_unknown_barcodes\t$counter_bad_quality_pool_barcodes\t$counter_bad_quality_well_barcodes\t$counter_bad_quality_RMT\t$counter_amp_batch_reads\t$counter_seq_batch_reads\n";
close(OUT_READ_STATS_AMP_BATCH);

my $oligos_header="";
for (my $oi=0;$oi<=$#oligos;$oi++){
  $oligos_header=$oligos_header.$oligos[$oi]."\t";
}
 
open(OUT_READ_STATS_PER_CELL, ">$read_stats_fn") || die "ERROR: cannot open file $read_stats_fn to read.\n";
print OUT_READ_STATS_PER_CELL "well_id\tlow_quality\tlow_mapq\tunmapped\tnon_unique_mapping\tmapped_to_nongenic\tspike_RMT_err\tspike_barcode_err\tspike_lonely_offset_err\tspike_readpoor_offset_err\tspike_mapped\tgene_RMT_err\tgene_barcode_err\tgene_lonely_offset_err\tgene_readpoor_offset_err\tgene_mapped\t".$oligos_header."total\n";

for (my $j=0; $j<=$#well_list;$j++){
  my $well_id=$well_list[$j];
  
  print OUT_READ_STATS_PER_CELL "$well_id\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"low_quality"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"lowmapq"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"unmapped"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"non_unique_mapping"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"mapped_to_nongenic"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"spike\tfilt1"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"spike\tfilt2"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"spike\tfilt3"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"spike\tfilt4"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"spike_mapped"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"gene\tfilt1"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"gene\tfilt2"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"gene\tfilt3"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"gene\tfilt4"})."\t";
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"gene_mapped"})."\t";
  for (my $oi=0;$oi<=$#oligos;$oi++){
	 print OUT_READ_STATS_PER_CELL (0+$oligo_counter{$well_id}->{$oligos[$oi]})."\t";
  }
  print OUT_READ_STATS_PER_CELL (0+$read_counter{$well_id}->{"total"})."\n";
}

close(OUT_READ_STATS_PER_CELL);


open(OUT_RMT_STATS_PER_CELL, ">$rmt_stats_fn") || die "ERROR: cannot open file $rmt_stats_fn to read.\n";
print OUT_RMT_STATS_PER_CELL "well_id\tfilt1\tfilt2\tfilt3\tfilt4\tok\tspike_filt1\tspike_filt2\tspike_filt3\tspike_filt4\tspike_ok\n";
for (my $j=0; $j<=$#well_list;$j++){
  my $well_id=$well_list[$j];
  print OUT_RMT_STATS_PER_CELL "$well_id\t".(0+$rmt_stats{"$well_id\tgene\t-1"})."\t".(0+$rmt_stats{"$well_id\tgene\t-2"})."\t".(0+$rmt_stats{"$well_id\tgene\t-3"})."\t".(0+$rmt_stats{"$well_id\tgene\t-4"})."\t".(0+$rmt_stats{"$well_id\tgene\t1"})."\t".(0+$rmt_stats{"$well_id\tspike\t-1"})."\t".(0+$rmt_stats{"$well_id\tspike\t-2"})."\t".(0+$rmt_stats{"$well_id\tspike\t-3"})."\t".(0+$rmt_stats{"$well_id\tspike\t-4"})."\t".(0+$rmt_stats{"$well_id\tspike\t1"})."\n";
}

close(OUT_RMT_STATS_PER_CELL);

open(OUT_NUC_PER_POS_STATS, ">$rmt_nuc_per_pos_fn") || die "ERROR: cannot open file $rmt_nuc_per_pos_fn to read.\n";
my $RMT_size=(keys %nuc_per_pos_stats)/4;
print OUT_NUC_PER_POS_STATS "position\tA\tC\tG\tT\n";
for (my $ni=0;$ni<$RMT_size;$ni++){
  print OUT_NUC_PER_POS_STATS $ni."\t".$nuc_per_pos_stats{$ni."\tA"}."\t".$nuc_per_pos_stats{$ni."\tC"}."\t".$nuc_per_pos_stats{$ni."\tG"}."\t".$nuc_per_pos_stats{$ni."\tT"}."\n";
}
close(OUT_NUC_PER_POS_STATS);

my @rmts=keys %rmt_status;
my %noffsets_per_rmt_hash={};
my %nreads_per_rmt_hash={};
open(OUT_DEBUG_RMTS, ">$debug_rmts_fn") || die "ERROR: annot open file $debug_rmts_fn to read.\n";
print OUT_DEBUG_RMTS "well_id\tgene\trmt\tis_spike\tstatus\tnreads\tn_offsets\n";
for (my $i=0;$i<=$#rmts;$i++){
  my $is_spike=0;
  my @arr=split("\t",$rmts[$i]);
  if (exists $spike_ins{$arr[1]}){
	$is_spike=1;
  }
  my $noffsets_per_cur_rmt=$offsets_per_rmt{$rmts[$i]};
  my $nreads_per_cur_rmt=$reads_per_rmt{$rmts[$i]};
 
  my $well_id=$arr[0];

  if ($rmt_status{$rmts[$i]}==1){
	if ($noffsets_per_cur_rmt>20){
	  $noffsets_per_cur_rmt=20
	}
	$noffsets_per_rmt_hash{$is_spike."\t".$number_of_cells{$well_id}."\t".$noffsets_per_cur_rmt}++;
	
	if ($nreads_per_cur_rmt>100){
	  $nreads_per_cur_rmt=100
	}
	$nreads_per_rmt_hash{$is_spike."\t".$number_of_cells{$well_id}."\t".$nreads_per_cur_rmt}++;
  }
  print OUT_DEBUG_RMTS $rmts[$i]."\t".$is_spike."\t".$rmt_status{$rmts[$i]}."\t".$reads_per_rmt{$rmts[$i]}."\t".$offsets_per_rmt{$rmts[$i]}."\n";
}
close(OUT_DEBUG_RMTS);

open(OFFSETS_PER_RMT, ">$noffsets_per_rmt_fn") || die "ERROR: cannot open file $noffsets_per_rmt_fn to read.\n";
print OFFSETS_PER_RMT "noffsets_per_rmt\tspikes\tspike_neg_ctrl\tgenes\tgenes_neg_ctrl\n";
for (my $i=0;$i<=20;$i++){
  print OFFSETS_PER_RMT $i."\t".(0+	$noffsets_per_rmt_hash{"1\t1\t$i"})."\t".(0+$noffsets_per_rmt_hash{"1\t0\t$i"})."\t".(0+$noffsets_per_rmt_hash{"0\t1\t$i"})."\t".(0+$noffsets_per_rmt_hash{"0\t0\t$i"})."\n";
}
close(OFFSETS_PER_RMT);

open(READS_PER_RMT, ">$nreads_per_rmt_fn") || die "ERROR: cannot open file $nreads_per_rmt_fn to read.\n";
print READS_PER_RMT "noffsets_per_rmt\tspikes\tspike_neg_ctrl\tgenes\tgenes_neg_ctrl\n";
for (my $i=0;$i<=100;$i++){
 print READS_PER_RMT $i."\t".(0+	$nreads_per_rmt_hash{"1\t1\t$i"})."\t".(0+$nreads_per_rmt_hash{"1\t0\t$i"})."\t".(0+$nreads_per_rmt_hash{"0\t1\t$i"})."\t".(0+$nreads_per_rmt_hash{"0\t0\t$i"})."\n";
}
close(READS_PER_RMT);


# filt1 RMT error
# filt2 BARCODE error
# filt3 lonely offset error
# filt4 read-poor offset error
