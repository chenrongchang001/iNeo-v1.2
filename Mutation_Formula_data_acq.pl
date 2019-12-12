#!/usr/bin perl -w
use strict;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

my$file1 = $ARGV[0]; # The final table 
my$file3 = $ARGV[1]; # The dbsnp.txt file
my$cancer_type = $ARGV[2];
my$sample_id = $ARGV[3];

die "\ninput: <the major table><*.dbsnp.txt> <cancer type> <ID>\n\n" if @ARGV < 4;

my($header,@col_names,%candidate_pos,%HLA_I_aff,%HLA_II_aff,%HLA_I_aff_rank,%HLA_II_aff_rank,%AC,$sample_work_dir,%BR_2_HLA);
my$max_len = 31;
my$min_len = 15;
my$col_num;
$sample_work_dir = `dirname $file1`; chomp $sample_work_dir;
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
#print "Input file is $opts{i}\nOutput file is $opts{o}\n";
#print "Database file is $opts{d}\n" if defined($opts{d});
#=================================READ INFO============================================
#my@pep_netMHCpan_BR_I = `cat $sample_work_dir/../../../pipeline/8.MHC/MHC_I/*.netMHCpan.raw.result | sed 's/<= WB//g' | sed 's/<= SB//g' | sed 's/<---->//g' | perl -lane '\@temp=split /\t/;print \$temp[2]' | perl -lane '\@temp=split / /;(\$hla,\$pep)=(\$1,\$2) if /(HLA-\\w\\*\\d+:\\d+)\\s+(\\w+)/;\@number;while(\$_=~ /(\\d+\\.\\d+)/g){push \@number,\$1};print "\$hla\t\$pep\t\$number[\$#number]"' | sort -u`;
#chomp @pep_netMHCpan_BR_I;
#my@pep_netMHCpan_BR_II = `cat $sample_work_dir/../../../pipeline/8.MHC/MHC_II/*.netMHCIIpan.raw.result | sed 's/<=WB//g' | sed 's/<=SB//g' | sed 's/<---->//g' | perl -lane '\@temp=split /\t/;print \$temp[2]' | perl -lane '(\$hla,\$pep)=(\$1,\$2) if /(D\\w+\\d+)\\s+(\\w+)/;\@number;while(\$_=~ /(\\d+\\.\\d+)/g){push \@number,\$1};print "\$hla\t\$pep\t\$number[\$#number-1]"'| sort -u`;
#chomp @pep_netMHCpan_BR_II;
#my(%pep_netMHCpan_BR_I,%pep_netMHCpan_BR_II);
#for(@pep_netMHCpan_BR_I){
#	my@temp=split /\t/;
#	$pep_netMHCpan_BR_I{$temp[1]} = $temp[2] if !exists $pep_netMHCpan_BR_I{$temp[1]};
#	$pep_netMHCpan_BR_I{$temp[1]} = $temp[2] if $temp[2] < $pep_netMHCpan_BR_I{$temp[1]};
#	$BR_2_HLA{$temp[1]} = $temp[0] if !exists $BR_2_HLA{$temp[1]} or $temp[2] < $pep_netMHCpan_BR_I{$temp[1]};
#}
#for(@pep_netMHCpan_BR_II){
#        my@temp=split /\t/;
#        $pep_netMHCpan_BR_II{$temp[1]} = $temp[2] if !exists $pep_netMHCpan_BR_II{$temp[1]};
#        $pep_netMHCpan_BR_II{$temp[1]} = $temp[2] if $temp[2] < $pep_netMHCpan_BR_II{$temp[1]};
#	$BR_2_HLA{$temp[1]} = $temp[0] if !exists $BR_2_HLA{$temp[1]} or $temp[2] < $pep_netMHCpan_BR_II{$temp[1]};
#}
#for(keys%BR_2_HLA){
#	print "$_\t$BR_2_HLA{$_}\n";
#}
open IN,$file1;
while(<IN>){
	if(/^Chr/){
		$header = $_;		
		@col_names = ("Chr","Position","Mutation_type","Gene_name","Gene_biotype","ENSG","AA_change","AF_ave","Mutation_level","BamCount_WES_normal_AF","BamCount_WES_tumor_AF","RNA_GATK_AF","RNA_Mutation_level","Bamcount_RNA_AF","Bamcount_RNA_cov","RPKM","HLA_I_type","wild_affinity_I(pssmhcpan_ic50-pssmhcpan_type:netMHC_aff-netMHC_bind:netMHCpan_aff-netMHCpan_bind:pickpocket_logaff-pickpocket_bind)","mut_pep_I","mut_affinity_I(pssmhcpan_ic50-pssmhcpan_type:netMHC_aff-netMHC_bind:netMHCpan_aff-netMHCpan_bind:pickpocket_logaff-pickpocket_bind)","wild_pep_all_I","mut_pep_all_I","wild_pep_all_II","mut_pep_all_II","NeoantigenQuality(1.5)","Wild_immunogenicity","Mut_immunogenicity","HLA_II_type","wild_affinity_II(iedb_rank:netMHCII_aff-netMHCII_bind:netMHCIIpan_aff-netMHCIIpan_bind)","mut_pep_II","mut_affinity_II(iedb_rank:netMHCII_aff-netMHCII_bind:netMHCIIpan_aff-netMHCIIpan_bind)");
		$col_num =  find_col_num($header,\@col_names);	
		next;
	}
	chomp;
	my($chr,$pos,$mut_type,$gene_name,$gene_biotype,$ensg,$aa_change,$AF,$Mute_level,$BamCount_WES_normal_AF,$BamCount_WES_tumor_AF,$RNA_GATK_AF,$RNA_Mutation_level,$Bamcount_RNA_AF,$Bamcount_RNA_cov,$rpkm,$HLA_I,$wild_aff_i,$mut_pep_I,$mut_aff_i,$wild_pep_all_I,$mut_pep_all_I,$wild_pep_all_II,$mut_pep_all_II,$neo_qual,$wild_IG,$mut_IG,$HLA_II,$wild_aff_ii,$mut_pep_II,$mut_aff_ii) = (split /\t/)[$$col_num[0],$$col_num[1],$$col_num[2],$$col_num[3],$$col_num[4],$$col_num[5],$$col_num[6],$$col_num[7],$$col_num[8],$$col_num[9],$$col_num[10],$$col_num[11],$$col_num[12],$$col_num[13],$$col_num[14],$$col_num[15],$$col_num[16],$$col_num[17],$$col_num[18],$$col_num[19],$$col_num[20],$$col_num[21],$$col_num[22],$$col_num[23],$$col_num[24],$$col_num[25],$$col_num[26],$$col_num[27],$$col_num[28],$$col_num[29],$$col_num[30]];	
	next if $mut_type =~ /Synonymous/ or $mut_type =~ /Splice/;
        next if $gene_biotype =~ /IG_/;
        next if $gene_biotype =~ /TR_/;

	my$enst = $1 if $_=~ /(ENST\d+)/;
	$Mute_level = 0.9*0.978 if $Mute_level eq "A";
	$Mute_level = 0.5*0.171 if $Mute_level eq "B";
	$Mute_level = 0.1*0.01 if $Mute_level eq "C";

#	my$hla_i_aff = process_aff_i($mut_aff_i) if $mut_aff_i ne ".";	
#	my$hla_ii_aff = process_aff_ii($mut_aff_ii) if $mut_aff_ii ne ".";

#	$HLA_I_aff{"$chr\t$pos\t$HLA_I\t$mut_aff_i"} = $hla_i_aff if $hla_i_aff;
#	$HLA_II_aff{"$chr\t$pos\t$HLA_II\t$mut_aff_ii"} = $hla_ii_aff if $hla_ii_aff;

	my$mutated_peptide = ($mut_pep_I eq ".")? $mut_pep_II : $mut_pep_I;
	my$pep_netMHCpan_BR = 1;
	my$wild_pep_BR = 0;
	$wild_pep_BR = 1 if $wild_aff_i =~ /WB|SB/ or $wild_aff_ii =~ /WB|SB/;
#	$pep_netMHCpan_BR = $pep_netMHCpan_BR_I{$mutated_peptide} if exists $pep_netMHCpan_BR_I{$mutated_peptide} and length($mutated_peptide) <= 11;
#	$pep_netMHCpan_BR = $pep_netMHCpan_BR_II{$mutated_peptide} if exists $pep_netMHCpan_BR_II{$mutated_peptide} and length($mutated_peptide) > 11;
	chomp $pep_netMHCpan_BR;

	if(length($mutated_peptide) <= 11 and $mut_aff_i =~ /SB/){$pep_netMHCpan_BR = 1}else{$pep_netMHCpan_BR = 1}
        if(length($mutated_peptide) > 11 and $mut_aff_i =~ /SB/){$pep_netMHCpan_BR = 1}else{$pep_netMHCpan_BR = 1}

#	$HLA_I = $BR_2_HLA{$mutated_peptide} if exists $BR_2_HLA{$mutated_peptide} and length($mutated_peptide) <= 11;
#	$HLA_II = $BR_2_HLA{$mutated_peptide} if exists $BR_2_HLA{$mutated_peptide} and length($mutated_peptide) > 11;


	my$M = 1;$M = 1 if $mut_type =~ /nonsynonymous/;
	if($mut_type =~ /nonframeshift/ and $aa_change =~ /\d+_\d+/){
		my($t1,$t2) = ($1,$2) if $aa_change =~ /(\d+)_(\d+)/;
		$M = $t2 -$t1 + 1;
	}
	if($mut_type =~ /frameshift/ and $mut_type !~ /non/){
		$M = length($mut_pep_all_I) - length($wild_pep_all_I) if $mut_pep_all_I ne ".";
		$M = length($mut_pep_all_II) - length($wild_pep_all_II) if $mut_pep_all_II ne ".";
	}
	my$TCGA_Exp= "";	
	$TCGA_Exp = `grep -w $cancer_type /data2/database/TCGA_download_expression/batch2/TCGA.exp.file | grep -w $ensg | cut -f 4 | cut -d " " -f 3`;$TCGA_Exp =~ s/\n//g;
	my$Pro_atlas_Exp = "";
	$Pro_atlas_Exp = `grep -w $chr $file3 | grep -w $pos | cut -f 33`;$Pro_atlas_Exp =~ s/\n//g;	
	my$temp = $1 if $Pro_atlas_Exp =~ /\((.*)\)\/\d+/;
	my@temp = split /,/,$temp if $temp;
	my$Pro_atlas_val = ".";
	$Pro_atlas_val = 10*$temp[0] + 6*$temp[1] + 1*$temp[2] if $temp;

	my$Driver_Gene_bonus = `grep -w $gene_name /data2/database/Cancer_Driver_Gene/*`;
	$Driver_Gene_bonus = 0 if length($Driver_Gene_bonus) < 2;
	$Driver_Gene_bonus = 1 if length($Driver_Gene_bonus) >= 2;
	push @{$candidate_pos{"$chr\t$pos\t$gene_name\t$Driver_Gene_bonus\t$aa_change\t$AF\t$BamCount_WES_normal_AF\t$BamCount_WES_tumor_AF\t$RNA_GATK_AF\t$Bamcount_RNA_AF\t$Bamcount_RNA_cov\t$M\t$Mute_level\t$RNA_Mutation_level\t$rpkm\t$TCGA_Exp\t$Pro_atlas_val\t$mutated_peptide\t"}},"$HLA_I\t$wild_aff_i\t$mut_aff_i\t$wild_pep_BR\t$pep_netMHCpan_BR\t$neo_qual\t$wild_IG\t$mut_IG\t" if $HLA_I ne "." and $HLA_I !~ /D/ and length($mutated_peptide) <= 11;	
	push @{$candidate_pos{"$chr\t$pos\t$gene_name\t$Driver_Gene_bonus\t$aa_change\t$AF\t$BamCount_WES_normal_AF\t$BamCount_WES_tumor_AF\t$RNA_GATK_AF\t$Bamcount_RNA_AF\t$Bamcount_RNA_cov\t$M\t$Mute_level\t$RNA_Mutation_level\t$rpkm\t$TCGA_Exp\t$Pro_atlas_val\t$mutated_peptide\t"}},"$HLA_II\t$wild_aff_ii\t$mut_aff_ii\t$wild_pep_BR\t$pep_netMHCpan_BR\t" if $HLA_II ne "." and $HLA_II=~ /D/ and length($mutated_peptide) > 11;
	
}
close IN;

open IN,$file3;
my$col_num_t;
while(<IN>){
	if(/Chr/){
		my$header = $_;
		my@col_names = ("pos_mhc1_stat_all(A,B,C,D)","pos_mhc2_stat_all(A,B,C,D)","Chr","Position");
		$col_num_t =  find_col_num($header,\@col_names);
		next;	
	}
	chomp;
	my($hla_i,$hla_ii,$chr,$pos) = (split /\t/)[$$col_num_t[0],$$col_num_t[1],$$col_num_t[2],$$col_num_t[3]];

	my@hla_i = split /\|/,$hla_i;
	my@hla_ii = split /\|/,$hla_ii;
	my$ac_i = 0;
	my$ac_ii = 0;
	my$ac = 0;
	if($hla_i ne "."){
		for(@hla_i){
			my@temp_i = split /,/;
			$ac_i += 5*$temp_i[0] + 3*$temp_i[1] + 1*$temp_i[2] + 0.5*$temp_i[3];	
		}
		$ac_i = $ac_i/@hla_i;
	}
	if($hla_ii ne "."){
                for(@hla_ii){
                        my@temp_ii = split /,/;
                        $ac_ii += 2.5*$temp_ii[0] + 1.5*$temp_ii[1] + 0.5*$temp_ii[2] + 0.25*$temp_ii[3];   
                }
                $ac_ii += $ac_ii/@hla_ii;
        }
	$ac = $ac_i + $ac_ii;
	$AC{"$chr\t$pos"} = $ac;
}
close IN;

#my@HLA_I_aff = values%HLA_I_aff;
#my@HLA_II_aff = values%HLA_II_aff;
#my@HLA_aff_all = (@HLA_I_aff,@HLA_II_aff);
#@HLA_aff_all = sort {$a<=>$b} @HLA_aff_all;

#for my$k(keys%HLA_I_aff){
#	for(my$i=0;$i<@HLA_aff_all;$i++){
#		$HLA_I_aff_rank{$k} = $i/@HLA_aff_all if $HLA_aff_all[$i] eq $HLA_I_aff{$k};		
#	}
#}
#for my$k(keys%HLA_II_aff){
#        for(my$i=0;$i<@HLA_aff_all;$i++){
#                $HLA_II_aff_rank{$k} = $i/@HLA_aff_all if $HLA_aff_all[$i] eq $HLA_II_aff{$k};
#        }
#}
#================================================================OUTPUT RESULT====================================================================================
open OUT,"> $sample_work_dir/$sample_id.formula.input";
print OUT "Chr\tPos\tGene\tDrive_Gene\tAA_change\tAF\tBamCount_WES_normal_AF\tBamCount_WES_tumor_AF\tRNA_GATK_AF\tBamcount_RNA_AF\tBamcount_RNA_cov\tM\tMut_level\tRNA_Mutation_level\tRPKM\tTCGA_Exp\tPro_atlas_Exp\tmutated_peptide\tHLA\tWild_aff\tMut_aff\tWild_BR\tNetMHCpan_BR\tNeo_qual\twild_IG\tmut_IG\tBR\tAC\tHLA_AB\tDNA_Mut_level\n";
for my$r(keys%candidate_pos){
	my@r = split /\t/,$r;
	for(@{$candidate_pos{$r}}){
		my@temp = split /\t/,$_;
#		print OUT "$r$_".$HLA_I_aff_rank{"$r[0]\t$r[1]\t$temp[0]\t$temp[2]"}."\t".$AC{"$r[0]\t$r[1]"}."\n" if $temp[0] =~ /HLA-[ABC]/;	
#		print OUT "$r$_.\t.\t.\t".$HLA_II_aff_rank{"$r[0]\t$r[1]\t$temp[0]\t$temp[2]"}."\t".$AC{"$r[0]\t$r[1]"}."\n" if $temp[0] =~ /D/;	
		$r =~ s/_/-/g;
		$_ =~ s/_/-/g;
		$AC{"$r[0]\t$r[1]"} =~ s/_/-/g;
		print OUT "$r$_"."none"."\t".$AC{"$r[0]\t$r[1]"}."\n" if $temp[0] =~ /HLA-[ABC]/;
		print OUT "$r$_.\t.\t.\t"."none"."\t".$AC{"$r[0]\t$r[1]"}."\n" if $temp[0] =~ /D/;
	}	
	print OUT "\n";
}


sub find_col_num{
        die "\nsub input error:<headers with \\t> <col_names in @> in find_col_num\n" if @_<2;
        my$header = shift;
        my$col_names = shift;
        my@header = split /\t/,$header;
        my@result;
        for my$c(@{$col_names}){
                for(my$i=0;$i<@header;$i++){
                        if($c eq $header[$i]){
                                push @result,$i;
                        }
                }
                }
        return \@result;
}

sub process_aff_i{
	my$input = shift;
	my@input = split /:/,$input;
	my@score;
	my$sum = 0;
	my$ave = 0;
	for(1..3){
		my$score = (split /-/,$input[$_])[0];
		$score = 2.718**((1-$score)*10.82) if $_ == 3;
		push @score,$score if $score =~ /\d/;
		$sum += $score if $score =~ /\d/;
	}		
	$ave = $sum/@score;		
	return $ave;
}

sub process_aff_ii{
        my$input = shift;
        my@input = split /:/,$input;
        my@score;
        my$sum = 0;
        my$ave = 0;
        for(@input){
		my$score = 0;
                $score = (split /-/,$_)[0] if $_ ne ".";
#                $score = 2.718**((1-$score)*10.82) if $_ == 3;
                push @score,$score if $score =~ /\d/;
                $sum += $score if $score =~ /\d/;
        }
        $ave = $sum/@score;
        return $ave;
}
