#!/usr/bin/perl

#################################################################################################################

# Copyright (C) 2015  Xiaowei Wang (email: xwang@radonc.wustl.edu), Washington University in St. Louis
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#################################################################################################################

use strict;
use warnings;
use Cwd;
use Sys::Hostname;

use File::Path;

my $version =         'V0.9';
my $file_dir =        './';
my $result_dir =      './result';
my $libsvm_dir =      './libsvm-2.82';
my $classifier_dir =  "./SVM_model";
my $inputFile;
my $outputFile = "WU-CRISPR_$version"."_prediction_result.xls";
my $minLength= 24;

my $result_file =     "$result_dir/sgOligo_$version" . "_prediction_result.xls";
my $feature_file =    "$classifier_dir/feature_list_v0.9.txt";

# Optional parameters for selection of sgRNA. Defaults to no restrictions.
my $scaffold_seq = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTT";
my $min_pos = 0;    
my $max_pos = 100000000; 
my $retained_portion = 1; 
my $scoreFilter = 1; 

################################### USER INPUTS ############################################################

my @inputs = @ARGV;
my $option = shift(@inputs);
&helpText if !@inputs;

my @sequences;
print "Welcome to WU-CRISPR.\n\n";
if ($option eq '-e' or $option eq '--example'){
        my $sampleSelection = shift @inputs;
        
        die ("Please select a valid sample option.\n") unless $sampleSelection eq 'short' or $sampleSelection eq 'long' or $sampleSelection eq 'multiple';
        $inputFile = "./samples/test_sequence_short.fasta" if $sampleSelection eq 'short';
        $inputFile = "./samples/test_sequence_long.fasta" if $sampleSelection eq 'long';
        $inputFile = "./samples/test_sequence_multiple.fasta" if $sampleSelection eq 'multiple';
        print "Selected file: $inputFile\n";
        @sequences = importFasta($inputFile);        
}
elsif($option eq '-f' or $option eq '--file'){
        $inputFile = shift @inputs;
        print "Selected file: $inputFile\n";
        if (-e $inputFile && -r $inputFile && -f $inputFile && -T $inputFile){ #check to ensure file exists, readable, plain text
                open (INPUTCHECK,$inputFile);
                while (<INPUTCHECK>){
                        s/\s+$//;
                        die ("Please ensure that the file is in FASTA format.\n") if $_ !~ /^>/;
                        last;
                }
                close (INPUTCHECK);
                # MAX - changing output file names
                $result_dir = $inputFile.".outDir";
                $outputFile = $inputFile.".outTab";
                # -- MAX
                @sequences = importFasta($inputFile);  
        }else{
                print "Error: Please check to make sure the file \"$inputFile\" exists and is a readable plain text file.\n";
                exit;       
        }
}
elsif ($option eq '-s' or $option eq '--sequence'){
        my $submission = shift @inputs;
        ${$sequences[0]}{'seq'} = $submission;
        my $seqLength = length $submission;
        ${$sequences[0]}{'id'} = "submittedSequence|length_$seqLength";
}
else{
        &helpText;
}
my $startTime = time();


print "\n******************** Genome-Wide sgOligo Version $version Standard Output **************************\n";

unlink $result_file if -e $result_file;
mkdir $result_dir if !(-e $result_dir);

my ($dH_ref, $dS_ref, $dG_ref, $dHi, $dSi, $dGi) = rna_dna_param();
my %dG = %{$dG_ref};

$scaffold_seq =~ tr/ATCGU/atcgu/;  $scaffold_seq =~ tr/t/u/;

my %p_val = import_p($feature_file);

my @feature;    
my @annotation; 
my %oligo_pos; 
my $gene_count = 0;
my $line_count = 0;
my $total_line_count = 0;
my $submittedSeq;
my $id;
my %id2Sequence;
foreach (@sequences){
        my $sequence = ${$_}{'seq'};
        $id = ${$_}{'id'};
        $id2Sequence{$id}=$sequence; 
        $submittedSeq = $sequence;
        $submittedSeq =~ tr/ATCGU/atcgu/; $submittedSeq =~ tr/u/t/;

        print "Error: Sequence contains bases other than A, T, C, G, or U. \n\tWU-CRISPR will now proceed to the next sequence.\n\n" and next if $submittedSeq =~/[^atcg]/i;
        print "Error: Sequence is shorter than $minLength bases. \n\tWU-CRISPR will now proceed to the next sequence.\n\n" and next if length ($submittedSeq)<$minLength;
        print "Error: Sequence is longer than 100,000 bases. \n\tWU-CRISPR will now now proceed to the next sequence.\n\n" and next if length ($submittedSeq)>100000;

        my $submittedSeq_rc = dnaComplement($submittedSeq);

        generate_feature($submittedSeq,"sense");
        generate_feature($submittedSeq_rc,"antisense");
}
predict(\@feature, \@annotation);

print "\n************* sgOligo selection process is done. Program completed successfully. *******************\n";

open(RESULT, $result_file) or die $!;
my %resultSeqs;
my %scoreList;
my $seqId;
while (<RESULT>) {
        s/\s+$//;
        next if $_ =~/^labels/;
        my @inline = split /\t/, $_;
        my $seq = substr($inline[8],0,20);
        my $seqSearch = $seq;
        my $orient = $inline[4];
        $seqId = $inline[3];
        $resultSeqs{$seq}{'orient'} = $orient;
        my $oligoSearch = $id2Sequence{$seqId};
        my $oligoLoc;
        $seqSearch = reverse($seq) and $seqSearch =~ tr/atcg/tagc/ if $orient eq 'antisense';
        my @pos1based;
        while ($oligoSearch =~ /$seqSearch/g) {
                $oligoLoc = pos($oligoSearch)-length($seqSearch)+1;
                push @pos1based,$oligoLoc;
        }
        
        my $location = join(", ",sort{$a<=>$b} @pos1based);
        $resultSeqs{$seq}{'location'} = $location;
        my $score= int(100*$inline[2]+0.5);
        push @{$scoreList{$seqId}}, join("\t",$score,$seq);
}

open(OUT, ">$outputFile") or die "$outputFile could not be opened for writing\n";

print OUT "seqId\tScore\tSequence\tOrientation\tPosition\n";
foreach my $sequenceId (sort keys %scoreList){
        
        print $result_file;
        foreach my $line (sort {$b cmp $a} @{$scoreList{$sequenceId}}){                my ($score,$seq) = split /\t/, $line;
                next if $score<50 and $scoreFilter ==1;
                my $direction = $resultSeqs{$seq}{'orient'};
                my $position = $resultSeqs{$seq}{'location'};
                print OUT "$sequenceId\t$score\t$seq\t$direction\t$position\n";
        }
        
}
close OUT;

my $endTime = time();
my $finalTime = $endTime-$startTime;
print"\nResults have been printed in $outputFile. Program completed in $finalTime seconds.\n";

########################################################################################################################################################################
sub generate_feature {

    my ($exon_plus,$orientation) = @_;
    
    my $dummyString = 'n'x50;
    $exon_plus = join("",$dummyString,$exon_plus,$dummyString);

    my $exon = substr($exon_plus, 50, length($exon_plus) - 100);

    my $seqStrand = $exon_plus;
    $seqStrand = dnaComplement($exon_plus) if $orientation eq 'antisense'; #return to the original sequence

    my @gg_pos_list;
    for (my $i =0; $i<length($exon);$i++){
        my $dinuc = substr($exon,$i,2);
        push @gg_pos_list, $i if $dinuc eq 'gg';
   }
    
    foreach my $gg_pos (@gg_pos_list) {               
            my $oligo = substr($exon_plus, $gg_pos + 50 - 21, 24);
            my $matched_bases = $oligo;
            $matched_bases = dnaComplement($oligo) if $orientation eq 'antisense';
            next if $matched_bases =~ /n/;
                        
            my $cds_pos = index($seqStrand,$matched_bases)-50+1;
            next unless $cds_pos > $min_pos and $cds_pos < $max_pos and $cds_pos / length($exon) < $retained_portion;
            my $output = generate_feature_line($oligo); # return '' if failed the pre-filters

            if ($output) {
                 $feature[$line_count] = $output;

                 my $exon_pos = $gg_pos;
                 $exon_pos = length($exon) - $gg_pos if $orientation eq "antisense";
                 $exon_pos += 1; # 1-based index position

                 $annotation[$line_count] = "$id\t$orientation\t$exon_pos\t$cds_pos\t".length($exon_plus)."\t$oligo";

                 $line_count++;
            }
            $total_line_count++;
    }
}

sub generate_feature_line {
    my $extSeq = shift;
    $extSeq =~ tr/t/u/;
    my $seq = substr($extSeq, 0, 20);

    my $output = '';
    my $feature_count = 0;

    # ************* Position-specific base composition *******************************************************************

    my $offset = -1; # use the extended sequence to include the PAM sequence
    my $exclude_pos = 0;

    for (my $indx = 0; $indx < length($extSeq); $indx++) {

          next if $indx == $offset + 1 or $indx == $offset + 2 or $indx == $offset + 4 or $indx == $offset + 5
                   or $indx == $offset + 6 or $indx == $offset + 7 or $indx == $offset + 8 or $indx == $offset + 13 or $indx == $offset + 22 or $indx == $offset + 23;

          my $base = substr($extSeq, $indx, 1);
          my $a_base = 0;
          my $c_base = 0;
          my $g_base = 0;
          my $ut_base = 0;

          if ($base eq 'a') {
             $a_base = 1;
          }
          elsif ($base eq 'c') {
                $c_base = 1;
          }
          elsif ($base eq 'g') {
                $g_base = 1;
          }
          elsif ($base eq 'u') {
                $ut_base = 1;
          }

          $output .= (++$feature_count) . ":$a_base " . (++$feature_count) . ":$c_base ";
          $output .= (++$feature_count) . ":$g_base " . (++$feature_count) . ":$ut_base ";

          $exclude_pos = 1 if $indx == 18 and $ut_base == 1;
          $exclude_pos = 1 if $indx == 19 and $ut_base == 1;
          $exclude_pos = 1 if $indx == 19 and $c_base == 1;

    }
    return ''  if $exclude_pos == 1;

    # ****************************************************************************************************************

    # ******************UUU (within last six bases); PolyX in the gRNA***********************************************

    my $polyA_reject_length = 5;
    my $polyC_reject_length = 5;
    my $polyG_reject_length = 4;
    my $polyU_reject_length = 4; 

    return '' if substr($seq, -6, 6) =~ /uuu/;
    return '' if $seq =~ /a{$polyA_reject_length}/ || $seq =~ /c{$polyC_reject_length}/ || $seq =~ /g{$polyG_reject_length}/ || $seq =~ /u{$polyU_reject_length}/;

    # ****************************************************************************************************************

    # ************** GC Content ********************************************************************************
    my $gc_content = gcPercent($seq);

    return '' if $gc_content > 0.8;

    $output .= (++$feature_count) . ":$gc_content ";

    # **********************************************************************************************************

    # ************** Duplex binding stability ******************************************************************
    my $binding_flag = 0;
    for (my $start = 0; $start < length($seq) - 4; $start++) {

          next unless $start == 0 or $start == 7;

          my $dG_binding = dG_binding(substr($seq, $start, length($seq)-$start));

          $output .= (++$feature_count) . ":$dG_binding ";

          $binding_flag = 1 if ($start == 0 and $dG_binding >= -18) or ($start == 7 and $dG_binding < -22) or ($start == 7 and $dG_binding > -9);
    }

    return '' if $binding_flag == 1;

    # ************************************************************************************************************

    # ************** Base accessibility (alignment) ******************************************************************
    my $dG_folding = foldingdG($seq);
        
    return '' if $dG_folding < -8;

    $output .= (++$feature_count) . ":$dG_folding ";

    my $gRNA = $seq . $scaffold_seq;
    my ($dG_folding1, $alignment) = RNA_fold($gRNA);

    my $ext_stem = "(((((((((.((((....))))...)))))))";
    my $aligned_stem = substr($alignment, 18, length($ext_stem));

    return '' if $aligned_stem eq $ext_stem;

    $alignment =~ tr/\.\(\)/011/;
    my @aligned = split "", $alignment;
    for (my $i = 0; $i <= $#aligned; $i++) {

          next unless $i == 14 or $i == 17 or $i == 18 or $i == 19 or $i == 20 or $i == 50 or $i == 51 or $i == 52;

          $output .= (++$feature_count) . ":" . $aligned[$i] . " ";
    }
    # ************************************************************************************************************

    # **************mono-, di- and tri-nucleotide compoisition *********************************************************
    my %base_num = (0 => 'a', 1 => 'u', 2 => 'c', 3 => 'g');
    my %dimer_count;
    my %trimer_count;
    my %tetramer_count;

    my @base_ary = $seq =~ /a/g; $output .= (++$feature_count) . ':' . scalar(@base_ary) . ' ';
    @base_ary = $seq =~ /u/g; $output .= (++$feature_count) . ':' . scalar(@base_ary) . ' ';
    @base_ary = $seq =~ /c/g; $output .= (++$feature_count) . ':' . scalar(@base_ary) . ' ';
    @base_ary = $seq =~ /g/g; $output .= (++$feature_count) . ':' . scalar(@base_ary) . ' ';

    for (my $pos = 0; $pos < length($seq)-1; $pos++) {
         my $dimer = substr($seq, $pos, 2);
         $dimer_count{$dimer}++;
    }
    for (my $i = 0; $i <= 3; $i++) {
         for (my $j = 0; $j <= 3; $j++) {
             my $dimer = $base_num{$i} . $base_num{$j};
             $dimer_count{$dimer} = 0 if !exists $dimer_count{$dimer};
             $output .= (++$feature_count) . ':' . $dimer_count{$dimer} . ' ';
         }
    }

    for (my $pos = 0; $pos < length($seq)-2; $pos++) {
         my $trimer = substr($seq, $pos, 3);
         $trimer_count{$trimer}++;
    }
    for (my $i = 0; $i <= 3; $i++) {
         for (my $j = 0; $j <= 3; $j++) {
              for (my $m = 0; $m <= 3; $m++) {
                        my $trimer = $base_num{$i} . $base_num{$j} . $base_num{$m};
                        $trimer_count{$trimer} = 0 if !exists $trimer_count{$trimer};
                        $output .= (++$feature_count) . ':' . $trimer_count{$trimer} . ' ';
             }
         }
    }

    # **********************************************************************************************************
    # All the features have been generated. Now it is time to generate the output line
    # **********************************************************************************************************

    $output = -1 . " " . $output;

    return $output;
}

sub predict {
    my ($feature_ref, $annotation_ref) = @_;
    my @feature = @{$feature_ref};

    open(PREDICT, ">$result_dir/sgOligo_predict.txt");
    for (my $j = 0; $j <= $#feature; $j++) {
         print PREDICT p_cutoff($feature[$j]) . "\n";
    }
    close(PREDICT);
    libsvm($annotation_ref);
}


sub import_p {
    my $file = shift;
    my %p;
    open(IN, $file) or die "Cannot open file for reading: $!\n";
    while (<IN>) {
         $_ =~ s/\s+$//;
         my @line = split /\t/, $_;
         $p{$line[0]} = $line[2];# if $line[2] < 0.01;  
    }
    close(IN);
    return %p;
}

sub p_cutoff {
    my $feat_line = shift;
    my @feature = split / /, $feat_line;

    my @output;
    for (my $i = 1; $i <= $#feature; $i++) {
         my ($index, $value) = split /:/, $feature[$i];
         push @output, $value if exists $p_val{$index};
    }
    my $out_string = $feature[0] . ' ';   # first include the class label
    for (my $j = 0; $j <= $#output; $j++) {
         $out_string .= ($j+1) . ':' . $output[$j] . ' ';
    }
    return $out_string;

}


sub libsvm {
    my ($annotation_ref) = @_;

    my $param_c;
    my $param_g;

    my $command1 = "$libsvm_dir/svm-scale -r $classifier_dir/sgOligo_train.range $result_dir/sgOligo_predict.txt > $result_dir/sgOligo_predict.scale";
    system($command1);

    my $command2 = "$libsvm_dir/svm-predict -b 1 $result_dir/sgOligo_predict.scale $classifier_dir/sgOligo_train.scale.model $result_dir/sgOligo.predict.xls";
    system($command2);

    my @prediction;
    my @annotation = @{$annotation_ref};

    my $reverse_order = 0;
    open(SVM, "$result_dir/sgOligo.predict.xls") or die "cannot open file1 for reading: $!\n";
    while(<SVM>){
           $_ =~ s/\s+$//;
           my @line = split / /, $_;

           if ($line[0] eq 'labels') {   
                if ($line[2] == -1) {
                     $reverse_order = 1;
                }
           }
           my $new_line = $_;
           $new_line = join (" ", ($line[0], $line[2], $line[1])) if $reverse_order == 1;
           push @prediction, $new_line;

    }
    close(SVM);

    open(OUT, ">$result_file");

    my @svm_header = split / /, $prediction[0];
    print OUT join("\t", @svm_header), "\tsequenceID\tOrientation\tPosition in Exon\tPosition in CDS\tCDS Length\tOligo Sequence\n";

    for (my $i = 1; $i <= $#prediction; $i++) {
         my @pred = split / /, $prediction[$i];
         if ($pred[1] >= 0) {
                print OUT  $pred[0] . "\t" . $pred[1] . "\t" . $pred[2] . "\t" . $annotation[$i-1] . "\n";
         }
    }
    close(OUT);
}


sub dG_binding {
# The overall deltaG of the target binding duplex. This is similar to GC content
    my $oligo = shift;
    $oligo =~ s/u/t/g;
    my $binding_dG = 0;
    for (my $indx = 0; $indx < length($oligo) - 1; $indx++) {
        next if substr($oligo, $indx, 2) =~ /n/;
         $binding_dG += $dG{substr($oligo, $indx, 2)};
    }
    $binding_dG += $dGi;

    return $binding_dG;
}

sub foldingdG {
   my $sequence = shift;
   my $tempSeq = "tempSeq$$";
   my $host = hostname();
   my $tempOUT = "$host.tempOUT$$";

   $sequence =~ s/[5|3|'|\-|\s+]//g;
   $sequence =~ tr/Tt/Uu/;

   open(OLIGO, ">$tempSeq") or die "can not open $tempSeq for writing: $!\n";
   print OLIGO $sequence;
   close OLIGO;

   my $dG;
   system("./RNAfold < $tempSeq > $tempOUT");
   open(RESULT, "$tempOUT") or die "Cannot open $tempOUT for reading $!\n";
   while(my $line = <RESULT>){
      if($line =~ /([\-|\d][\.|\d]+)\)/){
         $dG = $1;
         last;
      }
   }
   close(RESULT);
   unlink $tempSeq;
   unlink $tempOUT;
   return $dG;
}

sub RNA_fold {
   my $sequence = shift;
   my $tempSeq = "tempSeq$$";
   my $host = hostname();
   my $tempOUT = "$host.tempOUT$$";
   $sequence =~ s/[5|3|'|\-|\s+]//g;
   $sequence =~ tr/Tt/Uu/;
   
   open(OLIGO, ">$tempSeq") or die "can not open $tempSeq for writing: $!\n";
   print OLIGO $sequence;
   close OLIGO;

   my ($dG, $align);
   system("./RNAfold < $tempSeq > $tempOUT");
   open(RESULT, "$tempOUT") or die "Cannot open $tempOUT for reading $!\n";
   while(my $line = <RESULT>){
         if($line =~ /([\-|\d][\.|\d]+)\)/){
             $dG = $1;
             if($line =~ /^([\.\(\)]+)\s+/){
                 $align = $1;
             }
             else {
                   print $line, "\tUm... RNAfold alignment is empty!\n"; exit;
             }
             last;
         }
   }
   close(RESULT);
   unlink $tempSeq;
   unlink $tempOUT;
   return ($dG, $align);
}

sub dnaComplement {
    my ($sequence) = @_;
    $sequence =~ tr/atcgATCG/tagcTAGC/;
    $sequence = reverse($sequence);
    return $sequence;
}

sub gcPercent {
    my ($sequence) = @_;
    $sequence =~ tr/AUTCG/autcg/;
    my $length = length($sequence);
    my $gcCount= 0;

    for (my $j = 0; $j < $length; $j++) {
         if (substr($sequence, $j, 1) eq 'c' || substr($sequence, $j, 1) eq 'g') {
              $gcCount++;
         }
    }
    return sprintf("%.2f", $gcCount / $length);
}

sub rna_dna_param {
    my %dH = ();
    my %dS = ();
    my %dG = ();

    ($dH{'aa'}, $dS{'aa'}, $dG{'aa'}) = (-11.5, -36.4, -0.2);
    ($dH{'tt'}, $dS{'tt'}, $dG{'tt'}) = (-7.8, -21.9, -1.0);
    ($dH{'at'}, $dS{'at'}, $dG{'at'}) = (-8.3, -23.9, -0.9);
    ($dH{'ta'}, $dS{'ta'}, $dG{'ta'}) = (-7.8, -23.2, -0.6);
    ($dH{'ca'}, $dS{'ca'}, $dG{'ca'}) = (-10.4, -28.4, -1.6);
    ($dH{'tg'}, $dS{'tg'}, $dG{'tg'}) = (-9.0, -26.1, -0.9);
    ($dH{'ct'}, $dS{'ct'}, $dG{'ct'}) = (-9.1, -23.5, -1.8);
    ($dH{'ag'}, $dS{'ag'}, $dG{'ag'}) = (-7.0, -19.7, -0.9);
    ($dH{'ga'}, $dS{'ga'}, $dG{'ga'}) = (-8.6, -22.9, -1.5);
    ($dH{'tc'}, $dS{'tc'}, $dG{'tc'}) = (-5.5, -13.5, -1.3);
    ($dH{'gt'}, $dS{'gt'}, $dG{'gt'}) = (-5.9, -12.3, -2.1);
    ($dH{'ac'}, $dS{'ac'}, $dG{'ac'}) = (-7.8, -21.6, -1.1);
    ($dH{'cg'}, $dS{'cg'}, $dG{'cg'}) = (-16.3, -47.1, -1.7);
    ($dH{'gc'}, $dS{'gc'}, $dG{'gc'}) = (-8.0, -17.1, -2.7);
    ($dH{'gg'}, $dS{'gg'}, $dG{'gg'}) = (-9.3, -23.2, -2.1);
    ($dH{'cc'}, $dS{'cc'}, $dG{'cc'}) = (-12.8, -31.9, -2.9);

    my ($dHi, $dSi, $dGi) = (1.9, -3.9, 3.1);
    return (\%dH, \%dS, \%dG, $dHi, $dSi, $dGi);
}


sub fastaToTab {
     my ($fastaFile, $tabFile) = @_;
     my $id = "";
     my $dna= "";
     my $lastLine = "";
     open (IN, "$fastaFile") || die("Can not open $fastaFile file for reading in fastaToTab sub!\n");
     open (OUT, ">$tabFile") || die("Can not open $tabFile file for writing!\n");

     while (<IN>) {
          s/\s+$//;
          next if ($_ !~ /\S/);
          if ($_ =~ /^\>/) {
               $id = $_;
               $id =~ s/^\>//;
               if ($lastLine =~ /^\>/) {   
                    $id .= $_;
               }
               else {
                    print OUT "\n" if ($dna ne ""); 
                    print OUT "$id\t";
                    $id = "";
               }
          }
          else {
               $_ =~ s/\s//g;
               $dna = $_;
               print OUT $dna;
          }
          $lastLine = $_;
     }
     close(IN);
     close(OUT);
}

sub importTabSeq {
     my ($tabFile) = @_;
     my @sequence = ();
     my $index = 0;
     open (IN, "$tabFile") || die("Cannot open $tabFile file for reading in importTab sub!\n");
     while (<IN>) {
          s/\s+$//;
          my ($id, $sequence) = split /\t/, $_;
          $sequence[$index]{'id'} = $id;
          $sequence =~ tr/A-Z/a-z/;
          $sequence[$index]{'seq'} = $sequence;
          $index++;
     }
     close(IN);
     return @sequence;
}

sub importFasta {
    my ($fastaFile) = @_;
    my $tabFile = "$fastaFile $$.tab";
    fastaToTab($fastaFile, $tabFile);
    my @seq = importTabSeq($tabFile);
    unlink $tabFile if -e $tabFile;
    return @seq;
}

sub helpText{
        print "\n";
        
        print "USAGE:\n\tperl wu-crispr.pl [option] [path]\n\tperl wu-crispr.pl [option] [sequence]\n\n";
        print "SEQUENCE SUBMISSION:\n\t-s|--sequence <sequence>\n\t\t";
        print "Identifies sgRNA oligos from a single submitted sequence \n\t\tand provides a score for all potential active oligos. \n\t\tResults for submitted sequences will be printed to a \n\t\ttab-delimited text file.\n";
        print "\n\t\tExample: perl wu-crispr.pl -s attcatagagacacagagagaggt\n\n";
        
        print "FILE SUBMISSION:\n\t-f|--file <file>\n\t\t";
        print "Imports a FASTA file of sequences from <file> and \n\t\tidentifies potential sgRNA oligos for each submitted\n\t\tsequence. Resulting oligos are available in a tab-\n\t\tdelimited text file.\n";
        print "\n\t\tExample: perl wu-crispr.pl -f mySampleFile.fasta\n\n";
        
        print "EXAMPLE FILE SUBMISSION:\n\t-e|--example <short|long|multiple>\n\t\t";
        print "Uses the short, long, or multiple sample sequence in the \n\t\tsamples directory to generate sgRNA oligos.\n";
        print "\n\t\tExample: perl wu-crispr.pl -e short\n\n";
        
        print "HELP SCREEN:\n\t-h|--help \n\t\t";
        print "Brings up this help menu";
        print "\n\t\tExample: perl wu-crispr.pl -h\n\n";
        
        exit;
}
