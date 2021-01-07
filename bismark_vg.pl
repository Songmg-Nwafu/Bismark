#!/usr/bin/env perl
use strict;
use warnings;
use IO::Handle;
use Cwd;
$|++;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";

my $parent_dir = getcwd();
my $bismark_version = 'v0.22.3';
my $start_run = time();
my $command_line = join (" ",@ARGV);

### before processing the command line we will replace --solexa1.3-quals with --phred64-quals as the '.' in the option name will cause Getopt::Long to fail
foreach my $arg (@ARGV){
	if ($arg eq '--solexa1.3-quals'){
		$arg = '--phred64-quals';
	}
}

my @filenames;

my @fhs;         # stores alignment process names, bisulfite index location, bowtie filehandles and the number of times sequences produced an alignment
my %chromosomes; # stores the chromosome sequences of the mouse genome
my %SQ_order;    # stores the order of sequences in the reference. This is to produce SAM/BAM files with a known order of chromosomes
my %counting;    # counting various events
my $final_output_filename; # required for the nucleotide coverage report
my @pids; # storing the process IDs of child processes in parallel mode

my ($genome_folder,$CT_index_basename,$GA_index_basename,$path_to_vg,$aligner_options,$directional,$phred64,$output_dir,$vg,$bam,$samtools_path,$temp_dir,$multicore) =  process_command_line();

my $seqID_contains_tabs;
my $verbose = 0;

if ($multicore > 1){
	warn "Running Bismark Parallel version. Number of parallel instances to be spawned: $multicore\n\n";
}
##command line

foreach my $filename (@filenames){
	my $original_filename = $filename;
	my $original_filename_1;
	my $original_filename_2;
	
	chdir $parent_dir or die "Unable to move to initial working directory '$parent_dir' $!\n";
	
	### resetting the counting hash and fhs
	reset_counters_and_fhs($filename);
	@pid = ();
	$seqID_contains_tabs = 0;

	### if 2 or more files are provided we can hold the genome in memory and don't need to read it in a second time
	unless (%chromosomes){
		my $cwd = getcwd(); # storing the path of the current working directory
		warn "Current working directory is: $cwd\n\n";
		read_genome_into_memory($cwd);
	}
	
	# get general settings (also for single-threaded use)
	my ($pid,$offset) = multi_process_handling ();
	my ($single_end,$paired_end);
	
	### PAIRED-END ALIGNMENT
	if ($filename =~ ','){
		
		$single_end = 0;
		$paired_end = 1;
		
		my ($C_to_T_infile_1,$G_to_A_infile_1);
		$fhs[0]->{name} = 'CTread1GAread2CTgenome';
		$fhs[1]->{name} = 'GAread1CTread2GAgenome';
		$fhs[2]->{name} = 'GAread1CTread2CTgenome';
		$fhs[3]->{name} = 'CTread1GAread2GAgenome';
		warn "\nPaired-end alignments will be performed\n",'='x39,"\n\n";

		my ($filename_1,$filename_2) = (split (/,/,$filename));
		$original_filename_1 = $filename_1;
		$original_filename_2 = $filename_2;

		warn "The provided filenames for paired-end alignments are $filename_1 and $filename_2\n";

		### subsetting the input file(s)
		unless ($multicore == 1){ # not needed in single-core mode
			# warn "My PID: $pid\nMy offset: $offset\n";
			my $temp_filename_1 = subset_input_file_FastQ($filename_1,$pid,$offset);
			warn "Using the subset file >${temp_dir}$temp_filename_1< as new in-file 1 (instead of >$filename_1<)\n";

			my $temp_filename_2 = subset_input_file_FastQ($filename_2,$pid,$offset);
			warn "Using the subset file >${temp_dir}$temp_filename_2< as new in-file 2 (instead of >$filename_2<)\n";
			$filename_2 = "${temp_dir}$temp_filename_2";
		}

		### additional variables only for paired-end alignments
		my ($C_to_T_infile_2,$G_to_A_infile_2); # to be made from mate2 file

		my $read1_count; # to see if R1 and R2 have the same length
		my $read2_count;
		
		
		warn "Input files are in FastQ format\n";
		if ($directional){
			
			($C_to_T_infile_1,$read1_count) = biTransformFastQFiles_paired_end ($filename_1,1); # also passing the read number
			($G_to_A_infile_2,$read2_count) = biTransformFastQFiles_paired_end ($filename_2,2);
				
			$fhs[0]->{inputfile_1} = $C_to_T_infile_1;
			$fhs[0]->{inputfile_2} = $G_to_A_infile_2;
			$fhs[1]->{inputfile_1} = undef;
			$fhs[1]->{inputfile_2} = undef;
			$fhs[2]->{inputfile_1} = undef;
			$fhs[2]->{inputfile_2} = undef;
			$fhs[3]->{inputfile_1} = $C_to_T_infile_1;
			$fhs[3]->{inputfile_2} = $G_to_A_infile_2;
		}
		else{ # non-directional
			
			($C_to_T_infile_1,$G_to_A_infile_1,$read1_count) = biTransformFastQFiles_paired_end ($filename_1,1); # also passing the read number
			($C_to_T_infile_2,$G_to_A_infile_2,$read2_count) = biTransformFastQFiles_paired_end ($filename_2,2);

			$fhs[0]->{inputfile_1} = $C_to_T_infile_1;
			$fhs[0]->{inputfile_2} = $G_to_A_infile_2;
			$fhs[1]->{inputfile_1} = $G_to_A_infile_1;
			$fhs[1]->{inputfile_2} = $C_to_T_infile_2;
			$fhs[2]->{inputfile_1} = $G_to_A_infile_1;
			$fhs[2]->{inputfile_2} = $C_to_T_infile_2;
			$fhs[3]->{inputfile_1} = $C_to_T_infile_1;
			$fhs[3]->{inputfile_2} = $G_to_A_infile_2;
		}
			
		unless ($read1_count eq $read2_count){
			die "[FATAL ERROR]:\tNumber of bisulfite transformed reads are not equal between Read 1 (\#$read1_count) and Read 2 (\#$read2_count).\nPossible causes: file truncation, or as a result of specifying read pairs that do not belong to each other?! Please re-specify file names! Exiting...\n\n";
		}
						
		if ($bowtie2){
			paired_end_align_fragments_to_bisulfite_genome_fastQ_bowtie2 ($C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2);
		}
		else{ #
			paired_end_align_fragments_to_bisulfite_genome_fastQ_hisat2 ($C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2);
		}
	}
	
	### Else we are performing SINGLE-END ALIGNMENTS
	else{
		warn "\nSingle-end alignments will be performed\n",'='x39,"\n\n";

		$single_end = 1;
		$paired_end = 0;

		### subsetting the input file(s)
		unless ($multicore == 1){ # not needed in single-core mode
			# warn "My PID: $pid\nMy offset: $offset\n";
			# FastQ format, default
			my $temp_filename = subset_input_file_FastQ($filename,$pid,$offset);
			warn "Using the subset file >${temp_dir}$temp_filename< as new in-file (instead of >$filename<)\n";
			$filename = "${temp_dir}$temp_filename";
		}

		### Initialising bisulfite conversion filenames
		my ($C_to_T_infile,$G_to_A_infile);

		warn "Input file is in FastQ format\n";
		if ($directional){
		
			($C_to_T_infile) = biTransformFastQFiles ($filename);
			
			$fhs[0]->{inputfile} = $fhs[1]->{inputfile} = $C_to_T_infile;
		}	
		else{
			($C_to_T_infile,$G_to_A_infile) = biTransformFastQFiles ($filename);
			
			$fhs[0]->{inputfile} = $fhs[1]->{inputfile} = $C_to_T_infile;
			$fhs[2]->{inputfile} = $fhs[3]->{inputfile} = $G_to_A_infile;
		}

		### Creating up to 4 different filehandles and storing the first entry
		if ($pbat){
			if ($bowtie2){ # as of version 0.10.2 we also support PBAT alignments for Bowtie 2
				single_end_align_fragments_to_bisulfite_genome_fastQ_bowtie2 (undef,$G_to_A_infile);
			}
			else{ # HISAT2
				single_end_align_fragments_to_bisulfite_genome_fastQ_hisat2 (undef,$G_to_A_infile);
			}
		}
		elsif ($bowtie2){
			single_end_align_fragments_to_bisulfite_genome_fastQ_bowtie2 ($C_to_T_infile,$G_to_A_infile);
		}
		else{ # HISAT2
			single_end_align_fragments_to_bisulfite_genome_fastQ_hisat2 ($C_to_T_infile,$G_to_A_infile);
		}
		start_methylation_call_procedure_single_ends($filename,$C_to_T_infile,$G_to_A_infile,$pid);

	}
	### MERGING AND DELETING TEMP FILES // TIDYING UP AFTER A MULTICORE PROCESS

	if ($pid){ # only performing this for the parent process

		if ($multicore > 1){

			warn "Now waiting for all child processes to complete\n";

			### we need to ensure that we wait for all child processes to be finished before continuing
			# warn "here are the child IDs: @pids\n";
			# warn "Looping through the child process IDs:\n";

			foreach my $id (@pids){
				# print "$id\t";
				my $kid = waitpid ($id,0);
				# print "Returned: $kid\nExit status: $?\n";
				unless ($? == 0){
					warn "\nChild process terminated with exit signal: '$?'\n\n";
				}
			}

			# regenerating names for temporary files
			my @temp_input;
			my @temp_output;
			my @temp_reports;
			my @temp_unmapped_1;  # will store single end reads or R1 of paired-end
			my @temp_unmapped_2;
			my @temp_ambiguous_1; # will store single end reads or R1 of paired-end
			my @temp_ambiguous_2;
			my @temp_ambig_bam;

			for (1..$offset){

				# Temp Input Files
				if ($single_end){
					
					push @temp_input, "${original_filename}.temp.${_}";
					
				}
				elsif($paired_end){
					
					push @temp_input, "${original_filename_1}.temp.${_}";
					push @temp_input, "${original_filename_2}.temp.${_}";
					
				}

				# if files had a prefix we need to specify it
				my $add_prefix;
				if (defined $prefix){
					$add_prefix = "${prefix}.";
				}
				else{
					$add_prefix = '';
				}

				if ($single_end){
			 
					# Temp Output Files
					my $pathless_filename = ${original_filename}; # 10 Jan 2017
					$pathless_filename =~ s/.*\///; # deleting path information    
				
					push @temp_output,     "${output_dir}${add_prefix}${pathless_filename}.temp.${_}._bismark_vg.bam";
					push @temp_reports,    "${output_dir}${add_prefix}${pathless_filename}.temp.${_}._bismark_vg_SE_report.txt";
					
				}
				elsif($paired_end){
				
					# Temp Output Files
					my $pathless_filename_1 = ${original_filename_1}; # 10 Jan 2017
					my $pathless_filename_2 = ${original_filename_2};
					$pathless_filename_1 =~ s/.*\///; # deleting path information 
					$pathless_filename_2 =~ s/.*\///; # deleting path information 
					
					push @temp_output,     "${output_dir}${add_prefix}${pathless_filename}.temp.${_}._bismark_vg.bam";
					push @temp_reports,    "${output_dir}${add_prefix}${pathless_filename}.temp.${_}._bismark_vg_PE_report.txt";
							
				}
			}

			warn "\n\nRight, cleaning up now...\n\n";

			# deleting temp files;
			warn "Deleting temporary sequence files...\n";
			foreach my $temp (@temp_input){
				#print "$temp\t";
				$temp =~ s/.*\///; # deleting path information
				print "${temp_dir}${temp}\t";
				unlink "${temp_dir}${temp}" or warn "Failed to delete temporary FastQ file ${temp_dir}$temp: $!\n";
			}
			print "\n\n";

			# merging temp BAM files
			if ($single_end){
				merge_individual_BAM_files(\@temp_output,$original_filename,$single_end);
			}
			else{
				merge_individual_BAM_files(\@temp_output,$original_filename_1,$single_end);
			}

			# deleting temp BAM files
			warn "Deleting temporary BAM files...\n";
			foreach my $temp (@temp_output){
				# print "$temp\t";
				$temp =~ s/.*\///; # deleting path information
				print "${output_dir}${temp}\t";
				unlink "${output_dir}${temp}" or warn "Failed to delete temporary BAM file ${output_dir}${temp}: $!\n";
			}
			print "\n\n";

			# resetting the counters once more so we can add all data from all temporary reports
			reset_counters_and_fhs($original_filename);

			### Merging the Bismark mapping report files
			if ($single_end){
				merge_individual_mapping_reports(\@temp_reports,$original_filename,$single_end);
				print_final_analysis_report_single_end('mock_file1','mock_file_2','mock_pid','mergeThis');
			}
			else{
				merge_individual_mapping_reports(\@temp_reports,$original_filename_1,$single_end,$original_filename_2);
				print_final_analysis_report_paired_ends('mock_file1','mock_file_2','mock_file3','mock_file_4','mock_pid','mergeThis');
			}

			# deleting temp report files
			warn "Deleting temporary report files...\n";
			foreach my $temp (@temp_reports){
				print "$temp\t";
				unlink "${output_dir}${temp}" or warn "Failed to delete temporary report file $output_dir$temp: $!\n";
			}
			print "\n\n";
		}
	}

	if ($pid){ # only for the Parent
	
		### Produce Run Time
		my $end_run = time();
		my $run_time = $end_run - $start_run;
		my $days  = int($run_time/(24*60*60));
		my $hours = ($run_time/(60*60))%24;
		my $mins  = ($run_time/60)%60;
		my $secs  = $run_time%60;

		warn "Bismark completed in ${days}d ${hours}h ${mins}m ${secs}s\n";
		print REPORT "Bismark completed in ${days}d ${hours}h ${mins}m ${secs}s\n";
	
		warn "\n====================\nBismark run complete\n====================\n\n";

	}
	else{
		# If multiple files were supplied as the command line, like so:
		# -1 R1.fastq,simulated_1.fastq,ZZZ_R1.fastq -2 R2.fastq,simulated_2.fastq,ZZZ_R2.fastq --multicore 4
		# we need to exit from the child processes if we don't want a steady increase of new Bismark instances! Fixed 30 10 2017
		# warn "Terminating Child process\n\n"; 
		exit;
	}
}

sub reset_counters_and_fhs{
  my $filename = shift;
  %counting=(
	     total_meCHH_count => 0,
	     total_meCHG_count => 0,
	     total_meCpG_count => 0,
	     total_meC_unknown_count => 0,
	     total_unmethylated_CHH_count => 0,
	     total_unmethylated_CHG_count => 0,
	     total_unmethylated_CpG_count => 0,
	     total_unmethylated_C_unknown_count => 0,
	     sequences_count => 0,
	     no_single_alignment_found => 0,
	     unsuitable_sequence_count => 0,
	     genomic_sequence_could_not_be_extracted_count => 0,
	     unique_best_alignment_count => 0,
	     low_complexity_alignments_overruled_count => 0,
	     CT_CT_count => 0, #(CT read/CT genome, original top strand)
	     CT_GA_count => 0, #(CT read/GA genome, original bottom strand)
	     GA_CT_count => 0, #(GA read/CT genome, complementary to original top strand)
	     GA_GA_count => 0, #(GA read/GA genome, complementary to original bottom strand)
	     CT_GA_CT_count => 0, #(CT read1/GA read2/CT genome, original top strand)
	     GA_CT_GA_count => 0, #(GA read1/CT read2/GA genome, complementary to original bottom strand)
	     GA_CT_CT_count => 0, #(GA read1/CT read2/CT genome, complementary to original top strand)
	     CT_GA_GA_count => 0, #(CT read1/GA read2/GA genome, original bottom strand)
	     alignments_rejected_count => 0, # only relevant if --directional was specified
	    );

  if ($directional){
    if ($filename =~ ','){ # paired-end files
      @fhs=(
	    { name => 'CTreadCTgenome',
	      strand_identity => 'con ori forward',
	      bisulfiteIndex => $CT_index_basename,
	      seen => 0,
	      wrong_strand => 0,
	    },
	    { name => 'CTreadGAgenome',
	      strand_identity => 'con ori reverse',
	      bisulfiteIndex => $GA_index_basename,
	      seen => 0,
	      wrong_strand => 0,
	    },
	    { name => 'GAreadCTgenome',
	      strand_identity => 'compl ori con forward',
	      bisulfiteIndex => $CT_index_basename,
	      seen => 0,
	      wrong_strand => 0,
	    },
	    { name => 'GAreadGAgenome',
	    strand_identity => 'compl ori con reverse',
	      bisulfiteIndex => $GA_index_basename,
	      seen => 0,
	      wrong_strand => 0,
	    },
	   );
    }
    else{ # single-end files
      @fhs=(
	    { name => 'CTreadCTgenome',
	      strand_identity => 'con ori forward',
	      bisulfiteIndex => $CT_index_basename,
	      seen => 0,
	      wrong_strand => 0,
	    },
	    { name => 'CTreadGAgenome',
	      strand_identity => 'con ori reverse',
	      bisulfiteIndex => $GA_index_basename,
	      seen => 0,
	      wrong_strand => 0,
	    },
	   );
    }
  }
  else{
    @fhs=(
	  { name => 'CTreadCTgenome',
	    strand_identity => 'con ori forward',
	    bisulfiteIndex => $CT_index_basename,
	    seen => 0,
	    wrong_strand => 0,
	  },
	  { name => 'CTreadGAgenome',
	    strand_identity => 'con ori reverse',
	    bisulfiteIndex => $GA_index_basename,
	    seen => 0,
	    wrong_strand => 0,
	  },
	  { name => 'GAreadCTgenome',
	    strand_identity => 'compl ori con forward',
	    bisulfiteIndex => $CT_index_basename,
	    seen => 0,
	    wrong_strand => 0,
	  },
	  { name => 'GAreadGAgenome',
	    strand_identity => 'compl ori con reverse',
	    bisulfiteIndex => $GA_index_basename,
	    seen => 0,
	    wrong_strand => 0,
	  },
	 );
  }
}

sub read_genome_into_memory{

    ## working directoy
    my $cwd = shift;

    ## reading in and storing the specified genome in the %chromosomes hash
    chdir ($genome_folder) or die "Can't move to $genome_folder: $!";
    warn "Now reading in and storing sequence information of the genome specified in: $genome_folder\n\n";

    my @chromosome_filenames =  <*.fa>;

    ### if there aren't any genomic files with the extension .fa we will look for files with the extension .fa.gz
    unless (@chromosome_filenames){
	@chromosome_filenames =  <*.fa.gz>;
    }
    
    ### if there aren't any genomic files with the extension .fa or .fa.gz we will look for files with the extension .fasta
    unless (@chromosome_filenames){
	@chromosome_filenames =  <*.fasta>;
    }
    
    ### if there aren't any genomic files with the extension .fa or .fa.gz or .fasta we will look for files with the extension .fasta.gz
    unless (@chromosome_filenames){
	@chromosome_filenames =  <*.fasta.gz>;
    }
    
    unless (@chromosome_filenames){
	die "The specified genome folder $genome_folder does not contain any sequence files in FastA format (with .fa or .fasta file extensions, with or w/o .gz extension)\n";
    }

    my $SQ_count = 0;

    foreach my $chromosome_filename (@chromosome_filenames){
	# warn "Now processing: $chromosome_filename\n";
	if ($chromosome_filename =~ /\.gz$/){
	    open (CHR_IN,"gunzip -c $chromosome_filename |") or die "Failed to read from sequence file $chromosome_filename $!\n";
	}
	else{
	    open (CHR_IN,$chromosome_filename) or die "Failed to read from sequence file $chromosome_filename $!\n";
	}

	### first line needs to be a fastA header
	my $first_line = <CHR_IN>;
	chomp $first_line;
	$first_line =~ s/\r//;
	### Extracting chromosome name from the FastA header
	my $chromosome_name = extract_chromosome_name($first_line);
	if ($chromosome_name eq ''){ # should prevent chromosome name with spaces at the start such as > chr1, > chr2
	    die "Chromosome names must not be empty! Please check that there are no spaces at the start of the FastA header(s) and try again\n\n";
	}
	my $sequence;

	while (<CHR_IN>){
	  chomp;
	  $_ =~ s/\r//; # removing carriage returns if present
	  if ($_ =~ /^>/){

	    ### storing the previous chromosome in the %chromosomes hash, only relevant for Multi-Fasta-Files (MFA)
	    if (exists $chromosomes{$chromosome_name}){
	      print "chr $chromosome_name (",length $sequence ," bp)\n";
	      die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name!\n";
	    }
	    else {
	      if (length($sequence) == 0){
		warn "Chromosome $chromosome_name in the multi-fasta file $chromosome_filename did not contain any sequence information!\n";
	      }
	      print "chr $chromosome_name (",length $sequence ," bp)\n";
	      $chromosomes{$chromosome_name} = $sequence;
	      # warn "$SQ_count\t$chromosome_name\n";
	      ++$SQ_count;
	      $SQ_order{$SQ_count} = $chromosome_name;
	    }
	    ### resetting the sequence variable
	    $sequence = '';
	    ### setting new chromosome name
	    $chromosome_name = extract_chromosome_name($_);
	    if ($chromosome_name eq ''){ # should prevent chromosome name with spaces at the start such as > chr1, > chr2
		die "Chromosome names must not be empty! Please check that there are no spaces at the start of the FastA header(s) and try again.\n\n";
	    }
	  }
	  else{
	    $sequence .= uc$_;
	  }
	}

 	### Processing last chromosome of a multi Fasta File or the only entry in case of single entry FastA files

	if (exists $chromosomes{$chromosome_name}){
	    print "chr $chromosome_name (",length $sequence ," bp)\t";
	    die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name.\n";
	}
	else{
	    if (length($sequence) == 0){
		warn "Chromosome $chromosome_name in the file $chromosome_filename did not contain any sequence information!\n";
	    }

	    ++$SQ_count;
	    # warn "$SQ_count\t$chromosome_name\n";
	    print "chr $chromosome_name (",length $sequence ," bp)\n";
	    $chromosomes{$chromosome_name} = $sequence;
	    $SQ_order{$SQ_count} = $chromosome_name;
	}
    }
    print "\n";
    chdir $cwd or die "Failed to move to directory $cwd\n";

}

sub extract_chromosome_name {
    ## Bowtie seems to extract the first string after the inition > in the FASTA file, so we are doing this as well
    my $fasta_header = shift;
    if ($fasta_header =~ s/^>//){
			my ($chromosome_name) = split (/\s+/,$fasta_header);
			return $chromosome_name;
    }
    else{
			die "The specified chromosome ($fasta_header) file doesn't seem to be in FASTA format as required!\n";
    }
}

sub muti_process_handling{
	my $offset = 1;
	my $process_id;
	if ($multicore > 1){

		until ($offset == $multicore){
		# warn "multicore: $multicore\noffset: $offset\n";
			my $fork = fork;

			if (defined $fork){
				if ($fork != 0){
					$process_id = $fork;
					push @pids, $process_id;
					if ($offset < $multicore){
						++$offset;
					# warn "I am the parent process, child pid: $fork\nIncrementing offset counter to: $offset\n\n";
					}
					else{
						# warn "Reached the number of maximum multicores. Proceeeding to processing...\n";
					}
				}
				elsif ($fork == 0){
					# warn "I am a child process, pid: $fork\nOffset counter is: $offset\nProceeding to processing...\n";
					$process_id = $fork;
					last;
				}
			}
			else{
				die "[FATAL ERROR]: Forking unsuccessful. This normally means that something is fundamentally not working with the fork command. Please run again without the --parallel option, or ask your system admin to look into this.\n";
			}	
		}

		# warn "\nThe Thread Identity\n===================\n";
		if ($process_id){
			# print "I am the parent process. My children are called:\n";
			# print join ("\t",@pids),"\n";
			# print "I am going to process the following line count: $offset\n\n";
		}
		elsif($process_id == 0){
			# warn "I am a child process: Process ID: $process_id\n";
			# warn "I am going to process the following line count: $offset\n\n";
		}
		else{
			die "Process ID was: '$process_id'\n";
		}
	}
	else{
		warn "Single-core mode: setting pid to 1\n";
		$process_id = 1;
	}

	return ($process_id,$offset);
}

sub process_command_line{
	my @aligner_options;
	my $help;
	my $mates1;
	my $mates2;
	my $phred64;
	my $phred33;
	my $fastq;
	my $mismatches;
	my $seed_length;
	my $sequence_format;
	my $output_dir;
	my $parallel;
	my $multicore;
	my $vg;
	my $path_to_vg;
	my $samtools_path;
	my $genome_folder;
	my $non_directional;
	my $bam;
	
	my $command_line = Getoptions(
			'help|man'				  => \$help,
			'1=s'					  => \$mates1,
			'2=s'					  => \$mates2,
			'phred33-quals'           => \$phred33,
		    'phred64-quals'           => \$phred64,
			'n|seedmms=i'             => \$mismatches,
			'l|seedlen=i'             => \$seed_length,
			'parallel|multicore=i'    => \$multicore,
			'samtools_path=s'         => \$samtools_path,
			'non_directional'         => \$non_directional,
			'vg'					  => \$vg;
			'path_to_vg'			  => \$path_to_vg;	
			'o|output_dir=s'          => \$output_dir,
	);
	
	unless ($command_line){
		die "Please respecify command line options\n";
	}
	
	### HELPFILE
	if ($help){
		print_helpfile();
		exit;
	}
	
	### PROCESSING OPTIONS
	
	if ($vg){
		warn "You select the vg to align\n";
	}
	else {
		die "You may not select --vg. Make your pick! \n";
	}
	
	### Path to vg
	if (defined $path_to_vg){
		unless ($path_to_vg =~ /\/$/){
			$path_to_vg =~ s/$/\//;
		}
		if (-d $path_to_vg){
			$path_to_vg = "${path_to_vg}vg";
		}
		else{
			die "The path to vg provided ($path_to_vg) is invalid (not a directory)!\n";
		}
	}
	else {
		$path_to_vg = 'vg';
		warn "Path to vg specified as: $path_to_vg\n";
	}
	
	### BAM format
	$bam = 1;
	warn "Output format is BAM\n";
	
	### OUTPUT REQUESTED AS BAM FILE (default)
	if ($bam){
	
		### PATH TO SAMTOOLS
		if (defined $samtools_path){
			# if Samtools was specified as full command
			if ($samtools_path =~ /samtools$/){
				if (-e $samtools_path){
					# Samtools executable found
				}
				else{
					die "Could not find an installation of Samtools at the location $samtools_path. Please respecify\n";
				}	
			}
			else{
				unless ($samtools_path =~ /\/$/){
					$samtools_path =~ s/$/\//;
				}
				$samtools_path .= 'samtools';
				if (-e $samtools_path){
					# Samtools executable found
				}
				else{
					die "Could not find an installation of Samtools at the location $samtools_path. Please respecify\n";
				}
			}

			warn "Alignments will be written out in BAM format. Samtools path provided as: '$samtools_path'\n";
			$bam = 1;
		}
		# Check whether Samtools is in the PATH if no path was supplied by the user
		else{
			if (!system "which samtools >/dev/null 2>&1"){ # STDOUT is binned, STDERR is redirected to STDOUT. Returns 0 if samtools is in the PATH
				$samtools_path = `which samtools`;
				chomp $samtools_path;
				warn "Alignments will be written out in BAM format. Samtools found here: '$samtools_path'\n";
				$bam = 1;
			}
		}

		unless (defined $samtools_path){
			$bam = 2;
			warn "Did not find Samtools on the system. Alignments will be compressed with GZIP instead (.sam.gz)\n";
		}
		sleep (1);
	}
	
	### PROCESSING ARGUMENTS
	
	### GENOME FOLDER
	if (defined $genome_folder){
		# warn "Genome folder specified with --genome_folder $genome_folder\n";
	}
	else{
		$genome_folder = shift @ARGV; # mandatory
	}
	unless ($genome_folder){
		warn "Genome folder was not specified!\n";
		print_helpfile();
		exit;
	}

	### checking that the genome folder, all subfolders and the required bowtie index files exist
	unless ($genome_folder =~/\/$/){
		$genome_folder =~ s/$/\//;
	}

	if (chdir $genome_folder){
		my $absolute_genome_folder = getcwd(); ## making the genome folder path absolute
		unless ($absolute_genome_folder =~/\/$/){
			$absolute_genome_folder =~ s/$/\//;
		}
		warn "Reference genome folder provided is $genome_folder\t(absolute path is '$absolute_genome_folder)'\n";
		$genome_folder = $absolute_genome_folder;
	}
	else{
		die "Failed to move to $genome_folder: $!\nUSAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>} [<hits>]    (--help for more details)\n";
	}
	
	### index directory
	my $CT_dir = "${genome_folder}Bisulfite_Genome/CT_conversion/";
	my $GA_dir = "${genome_folder}Bisulfite_Genome/GA_conversion/";
	
	
	
	my $CT_index_basename = "${CT_dir}BS_CT";
	my $GA_index_basename = "${GA_dir}BS_GA";
	
	### INPUT OPTIONS
	$fastq = 1;
	warn "FastQ format assumed (by default)\n";
	$sequence_format = 'FASTQ';
	push @aligner_options, '-q';
	
	### QUALITY VALUES
	if ($phred33 and $phred64){
		die "You can only specify one type of quality value at a time! (--phred33 or --phred64)";
	}
	if ($phred33){ ## if nothing else is specified $phred33 will be used as default by both Bowtie 1 and 2.
		# Phred quality values work only when -q is specified
		unless ($fastq){
			die "Phred quality values works only when -q (FASTQ) is specified\n";
		}
		push @aligner_options,"--phred33";
	}

	if ($phred64){
		# Phred quality values work only when -q is specified
		unless ($fastq){
			die "Phred quality values work only when -q (FASTQ) is specified\n";
		}
		push @aligner_options,"--phred64";
	}
	else{
		$phred64 = 0;
	}
	
	### ALIGNMENT OPTIONS

	### MISMATCHES
	if (defined $mismatches){
		if ($mismatches == 0 or $mismatches == 1){
			push @aligner_options,"-N $mismatches";
		}
		else{
			die "Please set the number of multiseed mismatches with '-N <int>' (where <int> can be 0 or 1)\n";
		}
	}

	### SEED LENGTH
	if (defined $seed_length){
		push @aligner_options,"-L $seed_length";
	}
	
	
	### vg PARALLELIZATION OPTIONS
	
	
	### REPORTING OPTIONS

	push @aligner_options,'--ignore-quals'; ## All mismatches will receive penalty for mismatches as if they were of high quality, which is 6 by default

	### PAIRED-END MAPPING
	if ($mates1){

		if (defined $singles){ # if --single_end has been set explicitely
			die "You cannot set --single_end and supply files in paired-end format (-1 <mates1> -2 <mates2>). Please respecify!\n";
		}

		my @mates1 = (split (/,/,$mates1));
		die "Paired-end mapping requires the format: -1 <mates1> -2 <mates2>, please respecify!\n" unless ($mates2);
		my @mates2 = (split(/,/,$mates2));
		unless (scalar @mates1 == scalar @mates2){
			die "Paired-end mapping requires the same amounnt of mate1 and mate2 files, please respecify! (format: -1 <mates1> -2 <mates2>)\n";
		}
		while (1){
			my $mate1 = shift @mates1;
			my $mate2 = shift @mates2;
			last unless ($mate1 and $mate2);

			if ($mate1 eq $mate2){
				die "\n[FATAL ERROR]: Read 1 ($mate1) and Read 2 ($mate2) files were specified as the exact same file, which is almost certainly unintentional (and wrong). Please re-specify!\n\n";
			}
			push @filenames,"$mate1,$mate2";
		}
	}
	elsif ($mates2){
		die "Paired-end mapping requires the format: -1 <mates1> -2 <mates2>, please respecify!\n";
	}

	chdir $parent_dir or die "Failed to move back to parent directory: $!\n\n";
	my $current = getcwd();

	### SINGLE-END MAPPING
	# Single-end mapping will be performed if no mate pairs for paired-end mapping have been specified

	unless ($mates1 and $mates2){
		if (defined $singles){ # if --single_end has been set explicitely
			warn "Mapping set to single-end mode (user defined). File names need to be separated by commas [,] or colons [:]! Supplied file names are: $singles\n";
			$singles =~ s/:/,/g; # replacing colons (:) with commas
		}
		else{
			$singles = join (',',@ARGV);
			unless ($singles){
				die "\nNo filename supplied! Please specify one or more files for single-end Bismark mapping!\n";
			}
			$singles =~ s/\s/,/g; # replacing spaces with commas
		}
		@filenames = (split(/,/,$singles));
	}

	warn "\nInput files to be analysed (in current folder '$current'):\n";
	# checking if files exist and bail if they don't
	foreach my $f (@filenames){
		if ($mates1 and $mates2){
			my ($f1,$f2) = (split (/,/,$f));
			# mate 1
			if (-e $f1){
				warn "$f1\n";
			}
			else{
				die "Supplied filename '$f1' does not exist, please respecify!\n\n";
			}
			# mate 2
			if (-e $f2){
				warn "$f2\n";
			}
			else{
				die "Supplied filename '$f2' does not exist, please respecify!\n\n";
			}
		}
		else{
			if (-e $f){
				warn "$f\n";
			}
			else{
				die "Supplied filename '$f' does not exist, please respecify!\n\n";
			}
		}	
	}
	
	### STRAND-SPECIFIC LIBRARIES
	my $directional;
	if ($non_directional){
		die "A library can only be specified to be either non-directional or a PBAT-Seq library. Please respecify!\n\n" if ($pbat);
		warn "Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported\n";
		sleep (1);
		$directional = 0;
	}
	else{
		warn "Library is assumed to be strand-specific (directional), alignments to strands complementary to the original top or bottom strands will be ignored (i.e. not performed!)\n";
		sleep (1);
		$directional = 1; # default behaviour
	}
	
	### OUTPUT DIRECTORY
	
	chdir $parent_dir or die "Failed to move back to current working directory\n";
	if ($output_dir){
		unless ($output_dir =~ /\/$/){
			$output_dir =~ s/$/\//;
		}

		if (chdir $output_dir){
			$output_dir = getcwd; #  making the path absolute
			unless ($output_dir =~ /\/$/){
				$output_dir =~ s/$/\//;
			}
		}
		else{
			mkdir $output_dir or die "Unable to create directory $output_dir $!\n";
			warn "Created output directory $output_dir!\n\n";
			chdir $output_dir or die "Failed to move to $output_dir\n";
			$output_dir = getcwd; #  making the path absolute
			unless ($output_dir =~ /\/$/){
				$output_dir =~ s/$/\//;
			}
		}
		warn "Output will be written into the directory: $output_dir\n";
	}
	else{
		$output_dir = '';
	}

	### TEMPORARY DIRECTORY for C->T and G->A transcribed files

	chdir $parent_dir or die "Failed to move back to current working directory\n";
	if ($temp_dir){
		warn "\nUsing temp directory: $temp_dir\n";
		unless ($temp_dir =~ /\/$/){
			$temp_dir =~ s/$/\//;
		}

		if (chdir $temp_dir){
			$temp_dir = getcwd; #  making the path absolute
			unless ($temp_dir =~ /\/$/){
				$temp_dir =~ s/$/\//;
			}
		}
		else{
			mkdir $temp_dir or die "Unable to create directory $temp_dir $!\n";
			warn "Created temporary directory $temp_dir!\n\n";
			chdir $temp_dir or die "Failed to move to $temp_dir\n";
			$temp_dir = getcwd; #  making the path absolute
			unless ($temp_dir =~ /\/$/){
				$temp_dir =~ s/$/\//;
			}
		}
		warn "Temporary files will be written into the directory: $temp_dir\n";
	}
	else{
		$temp_dir = '';
	}
	
	my $aligner_options = join (' ',@aligner_options);
	warn "Summary of all aligner options:\t$aligner_options\n"; sleep(2);
	
	return ($genome_folder,$CT_index_basename,$GA_index_basename,$path_to_vg,$aligner_options,$directional,$phred64,$output_dir,$vg,$bam,$samtools_path,$temp_dir,$multicore);
}

sub print_helpfile{
	print << "HOW_TO";

DESCRIPTION
	
ARGUMENTS:

<genome_folder>          The path to the folder containing the unmodified reference genome
                         as well as the subfolders created by the Bismark_Genome_Preparation
                         script (/Bisulfite_Genome/CT_conversion/ and /Bisulfite_Genome/GA_conversion/).
                         Bismark expects one or more fastA files in this folder (file extension: .fa, .fa.gz
                         or .fasta or .fasta.gz). The path can be relative or absolute. The path may also be set
                         as '--genome_folder /path/to/genome/folder/'.

-1 <mates1>              Comma-separated list of files containing the #1 mates (filename usually includes
                         "_1"), e.g. flyA_1.fq,flyB_1.fq). Sequences specified with this option must
                         correspond file-for-file and read-for-read with those specified in <mates2>.
                         Reads may be a mix of different lengths. Bismark will produce one mapping result
                         and one report file per paired-end input file pair.

-2 <mates2>              Comma-separated list of files containing the #2 mates (filename usually includes
                         "_2"), e.g. flyA_2.fq,flyB_2.fq). Sequences specified with this option must
                         correspond file-for-file and read-for-read with those specified in <mates1>.
                         Reads may be a mix of different lengths.

<singles>                A comma- or space-separated list of files containing the reads to be aligned (e.g.
                         lane1.fq,lane2.fq lane3.fq). Reads may be a mix of different lengths. Bismark will
                         produce one mapping result and one report file per input file. Please note that
                         one should supply a list of files in conjunction with --basename as the output files
                         will constantly overwrite each other...

OPTIONS:

Input:

--se/--single_end <list> Sets single-end mapping mode explicitly giving a list of file names as <list>.
                         The filenames may be provided as a comma [,] or colon [:] separated list.

-q/--fastq               The query input files (specified as <mate1>,<mate2> or <singles> are FASTQ
                         files (usually having extension .fg or .fastq). This is the default. See also
                         --solexa-quals.

-f/--fasta               The query input files (specified as <mate1>,<mate2> or <singles> are FASTA
                         files (usually having extensions .fa, .mfa, .fna or similar). All quality values
                         are assumed to be 40 on the Phred scale. FASTA files are expected to contain both
                         the read name and the sequence on a single line (and not spread over several lines).

-s/--skip <int>          Skip (i.e. do not align) the first <int> reads or read pairs from the input.

-u/--upto <int>          Only aligns the first <int> reads or read pairs from the input. Default: no limit.

--phred33-quals          FASTQ qualities are ASCII chars equal to the Phred quality plus 33. Default: ON.

--phred64-quals          FASTQ qualities are ASCII chars equal to the Phred quality plus 64. Default: off.

--path_to_bowtie2        The full path </../../> to the Bowtie 2 installation on your system. If not
                         specified it is assumed that Bowtie 2 is in the PATH.

--path_to_hisat2         The full path </../../> to the HISAT2 installation on your system. If not
                         specified it is assumed that HISAT2 is in the PATH.

	
	
Alignment:


-N <int>                 Sets the number of mismatches to allowed in a seed alignment during multiseed alignment.
                         Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower)
                         but increases sensitivity. Default: 0. This option is only available for Bowtie 2 (for
                         Bowtie 1 see -n).

-L <int>                 Sets the length of the seed substrings to align during multiseed alignment. Smaller values
                         make alignment slower but more senstive. Default: the --sensitive preset of Bowtie 2 is
                         used by default, which sets -L to 20. maximum of L can be set to 32. The length of the seed
                         would effect the alignment speed dramatically while the larger L, the faster the aligment.
                         This option is only available for Bowtie 2 (for Bowtie 1 see -l).

--parallel <int>         (May also be --multicore <int>) Sets the number of parallel instances of Bismark to be run concurrently.
                         This forks the Bismark alignment step very early on so that each individual Spawn of Bismark processes
                         only every n-th sequence (n being set by --parallel). Once all processes have completed,
                         the individual BAM files, mapping reports, unmapped or ambiguous FastQ files are merged
                         into single files in very much the same way as they would have been generated running Bismark
                         conventionally with only a single instance.

                         If system resources are plentiful this is a viable option to speed up the alignment process
                         (we observed a near linear speed increase for up to --parallel 8 tested). However, please note
                         that a typical Bismark run will use several cores already (Bismark itself, 2 or 4 threads of
                         Bowtie2/HISAT2, Samtools, gzip etc...) and ~10-16GB of memory depending on the choice of aligner
                         and genome. WARNING: Bismark Parallel (BP?) is resource hungry! Each value of --parallel specified
                         will effectively lead to a linear increase in compute and memory requirements, so --parallel 4 for
                         e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~40GB or RAM, but at the same time
                         reduce the alignment time to ~25-30%. You have been warned.

	
Output:

-o/--output_dir <dir>    Write all output files into this directory. By default the output files will be written into
                         the same folder as the input file(s). If the specified folder does not exist, Bismark will attempt
                         to create it first. The path to the output folder can be either relative or absolute.

--temp_dir <dir>         Write temporary files to this directory instead of into the same directory as the input files. If
                         the specified folder does not exist, Bismark will attempt to create it first. The path to the
                         temporary folder can be either relative or absolute.

Bismark BAM/SAM OUTPUT (default):

 (1) QNAME  (seq-ID)
 (2) FLAG   (this flag tries to take the strand a bisulfite read originated from into account (this is different from ordinary DNA alignment flags!))
 (3) RNAME  (chromosome)
 (4) POS    (start position)
 (5) MAPQ   (always 255 for use with Bowtie)
 (6) CIGAR
 (7) RNEXT
 (8) PNEXT
 (9) TLEN
(10) SEQ
(11) QUAL   (Phred33 scale)
(12) NM-tag (edit distance to the reference)
(13) MD-tag (base-by-base mismatches to the reference (handles indels)
(14) XM-tag (methylation call string)
(15) XR-tag (read conversion state for the alignment)
(16) XG-tag (genome conversion state for the alignment)
(17) XA/XB-tag (non-bisulfite mismatches) (optional!)

Each read of paired-end alignments is written out in a separate line in the above format.


Last modified on 26 July 2019
HOW_TO
}

sub subset_input_file_FastQ{

	my ($filename,$process_id,$offset) = @_;

	if ($filename =~ /gz$/){
		open (OFFSET,"gunzip -c $filename |") or die "Couldn't read from file '$filename': $!\n";
	}
	else{
		open (OFFSET,$filename) or die "Couldn't read from file '$filename': $!\n";
	}

	# warn "offset is $offset\n";
	my $temp = $filename;
	$temp .= ".temp.$offset";
	$temp =~ s/^.*\///; # replacing everything upto and including the last /, i.e. removing file path information

	open (TEMPFQ,'>',"${temp_dir}${temp}") or die "Failed to write output ${temp_dir}${temp}: $!\n";

	my $line_count = 0;
	my $seqs_processed = 0;
	
	while (1){
		my $l1 = <OFFSET>;
		my $l2 = <OFFSET>;
		my $l3 = <OFFSET>;
		my $l4 = <OFFSET>;

		last unless ($l4);
		++$line_count;

		if ( ($line_count - $offset)%$multicore == 0){
			# warn "line count: $line_count\noffset: $offset\n";
			# warn "Modulus: ",($line_count - $offset)%$multicore,"\n";
			# warn "processing this line $line_count (processID: $process_id with \$offset $offset)\n";
			print TEMPFQ "$l1$l2$l3$l4";
			$seqs_processed++;
		}
		else{
			# warn "skipping line $line_count for processID: $process_id with \$offset $offset)\n";
			next;
		}
	}

	close OFFSET; # or warn $!;
	close TEMPFQ or warn "Failed to close file handle TEMPFQ: $!\n";

	warn "Finished subdividing $filename for PID: $process_id and offset $offset (sequences written out: $seqs_processed)\n\n";

	return ($temp); # returning the subset filename
}

sub biTransformFastQFiles_paired_end {
	my ($file,$read_number) = @_;
	my ($dir,$filename);

	if ($file =~ /\//){
		($dir,$filename) = $file =~ m/(.*\/)(.*)$/;
	}
	else{
		$filename = $file;
	}

	### gzipped version of the infile
	if ($file =~ /\.gz$/){
		open (IN,"gunzip -c $file |") or die "Couldn't read from file $file: $!\n";
	}
	else{
		open (IN,$file) or die "Couldn't read from file $file: $!\n";
	}

		
	my $C_to_T_infile = my $G_to_A_infile = $filename;

	
	$C_to_T_infile =~ s/$/_C_to_T.fastq/;
	$G_to_A_infile =~ s/$/_G_to_A.fastq/;

	

	if ($directional){
		if ($read_number == 1){
			warn "Writing a C -> T converted version of the input file $filename to $temp_dir$C_to_T_infile\n";
			
			open (CTOT,'>',"$temp_dir$C_to_T_infile") or die "Couldn't write to file $!\n";
		}
		elsif ($read_number == 2){
			warn "Writing a G -> A converted version of the input file $filename to $temp_dir$G_to_A_infile\n";
			
			open (GTOA,'>',"$temp_dir$G_to_A_infile") or die "Couldn't write to file $!\n";
			
		}
		else{
			die "Read number needs to be 1 or 2, but was $read_number!\n\n";
		}
	}
	else{ # non-directional
		warn "Writing a C -> T converted version of the input file $filename to $temp_dir$C_to_T_infile\n";
		warn "Writing a G -> A converted version of the input file $filename to $temp_dir$G_to_A_infile\n";
		
		open (CTOT,'>',"$temp_dir$C_to_T_infile") or die "Couldn't write to file $!\n";
		open (GTOA,'>',"$temp_dir$G_to_A_infile") or die "Couldn't write to file $!\n";
			
	}	
	

	my $count = 0;
	while (1){
		my $identifier = <IN>;
		my $sequence = <IN>;
		my $identifier2 = <IN>;
		my $quality_score = <IN>;
		last unless ($identifier and $sequence and $identifier2 and $quality_score);
		++$count;

		chomp $identifier;
		$identifier = fix_IDs($identifier); # this is to avoid problems with truncated read ID when they contain white spaces
		$identifier .= "\n";
		
		if ($skip){
			next unless ($count > $skip);
		}
		if ($upto){
			last if ($count > $upto);
		}

		$sequence= uc$sequence; # make input file case insensitive

		## small check if the sequence file appears to be a FastQ file
		if ($count == 1){
			if ($identifier !~ /^\@/ or $identifier2 !~ /^\+/){
				die "Input file doesn't seem to be in FastQ format at sequence $count: $!\n";
			}	
		}
		my $sequence_C_to_T = my $sequence_G_to_A = $sequence;

		if ($read_number == 1){
			$identifier =~ s/$/\/1\/1/;
		}
		elsif ($read_number == 2){
			$identifier =~ s/$/\/2\/2/;
		}
		else{
			die "Read number needs to be 1 or 2\n";
		}

		
		if ($directional){
			if ($read_number == 1){
				$sequence_C_to_T =~ tr/C/T/;
				print CTOT join ('',$identifier,$sequence_C_to_T,$identifier2,$quality_score);
			}
			else{
				$sequence_G_to_A =~ tr/G/A/;
				print GTOA join ('',$identifier,$sequence_G_to_A,$identifier2,$quality_score);
			}	
		}
		else{
			$sequence_C_to_T =~ tr/C/T/;
			$sequence_G_to_A =~ tr/G/A/;			
			print CTOT join ('',$identifier,$sequence_C_to_T,$identifier2,$quality_score);
			print GTOA join ('',$identifier,$sequence_G_to_A,$identifier2,$quality_score);
		}	
	}

	if ($directional){
		if ($read_number == 1){
			warn "\nCreated C -> T converted version of the FastQ file $filename ($count sequences in total)\n\n";
			close CTOT or die "Failed to close filehandle $!\n";
			return ($C_to_T_infile,$count); # passing back the number of transliterated sequences so we can make sure R1 and R2 files have the same length
		}
		else{
			warn "\nCreated G -> A converted version of the FastQ file $filename ($count sequences in total)\n\n";
			close GTOA or die "Failed to close filehandle $!\n";
			return ($G_to_A_infile,$count);
		}
	}
	else{
		warn "\nCreated C -> T as well as G -> A converted versions of the FastQ file $filename ($count sequences in total)\n\n";
		close CTOT or die "Failed to close filehandle $!\n";
		close GTOA or die "Failed to close filehandle $!\n";
		return ($C_to_T_infile,$G_to_A_infile,$count);  # passing back the number of transliterated sequences so we can make sure R1 and R2 files have the same length
	}
}

sub paired_end_align_fragments_to_bisulfite_genome_fastQ_hisat2 {

	my ($C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2) = @_;
	if ($directional){
		warn "Input files are $C_to_T_infile_1 and $G_to_A_infile_2 (FastQ)\n";
	}
	elsif ($pbat){
		warn "Input files are $G_to_A_infile_1 and $C_to_T_infile_2 (FastQ)\n";
	}
	else{
		warn "Input files are $C_to_T_infile_1 and $G_to_A_infile_1 and $C_to_T_infile_2 and $G_to_A_infile_2 (FastQ)\n";
	}

	## Now starting up 4 instances of HISAT2 feeding in the converted sequence files and reading in the first line of the output, and storing it in the
	## data structure above
	if ($directional or $pbat){
		warn "Now running 2 instances of HISAT2 against the bisulfite genome of $genome_folder with the specified options: $aligner_options\n\n";
	}
	else{
		warn "Now running 4 individual instances of HISAT2 against the bisulfite genome of $genome_folder with the specified options: $aligner_options\n\n";
	}

	foreach my $fh (@fhs) {

		if ($directional or $pbat){ # skipping unwanted filehandles
			unless ($fh->{inputfile_1}){
				$fh->{last_seq_id} = undef;
				$fh->{last_line_1} = undef;
				$fh->{last_line_2} = undef;
				next;
			}
		}

		my $hisat2_options = $aligner_options;
		if ($fh->{name} eq 'CTread1GAread2CTgenome' or $fh->{name} eq 'GAread1CTread2GAgenome'){
			$hisat2_options .= ' --norc'; ### ensuring the alignments are only reported in a sensible manner
		}
		else {
			$hisat2_options .= ' --nofw';
		}

		warn "Now starting a HISAT2 paired-end alignment for $fh->{name} (reading in sequences from $temp_dir$fh->{inputfile_1} and $temp_dir$fh->{inputfile_2}, with the options: $hisat2_options))\n";
		open ($fh->{fh},"$path_to_hisat2 $hisat2_options -x $fh->{bisulfiteIndex} -1 $temp_dir$fh->{inputfile_1} -2 $temp_dir$fh->{inputfile_2} |") or die "Can't open pipe to HISAT2: $!";

		### HISAT2 outputs out SAM format, so we need to skip everything until the first sequence
		while (1){
			$_ = $fh->{fh}->getline();
			if ($_) {
				last unless ($_ =~ /^\@/); # SAM headers start with @
			}
			else{
				last; # no alignment output
			}
		}

		my $line_1 = $_;
		my $line_2 = $fh->{fh}->getline();

		# if HISAT2 produces an alignment we store the first line of the output
		if ($line_1 and $line_2) {
			chomp $line_1;
			chomp $line_2;
			
			# Need to make sure HISAT2 does the same
			### Bowtie always reports the alignment with the smaller chromosomal position first. This can be either sequence 1 or sequence 2.
			### We will thus identify which sequence was read 1 and store this ID as last_seq_id

			my $id_1 = (split(/\t/,$line_1))[0]; # this is the first element of the first HISAT2 output line (= the sequence identifier)
			my $id_2 = (split(/\t/,$line_2))[0]; # this is the first element of the second HISAT2 output line
			
			if ($id_1 =~ s/\/1$//){ # removing the read 1 tag if present (remember that HISAT2 clips off /1 or /2 line endings itself, so we added /1/1 or /2/2 to start with
				$fh->{last_seq_id} = $id_1;
			}
			elsif ($id_2 =~ s/\/1$//){ # removing the read 1 tag if present
				$fh->{last_seq_id} = $id_2;
			}
			else{
				die "Either the first or the second id need to be read 1! ID1 was: $id_1; ID2 was: $id_2\n";
			}
			
			$fh->{last_line_1} = $line_1; # this contains read 1 or read 2
			$fh->{last_line_2} = $line_2; # this contains read 1 or read 2
			# warn "Set last_seq_id as: $fh->{last_seq_id}\n";
			warn "Found first alignment:\n$fh->{last_line_1}\n$fh->{last_line_2}\n";
		}
		# otherwise we just initialise last_seq_id and last_lines as undefined
		else {
			warn "Found no alignment, assigning undef to last_seq_id and last_lines\n";
			$fh->{last_seq_id} = undef;
			$fh->{last_line_1} = undef;
			$fh->{last_line_2} = undef;
		}
	}
}

sub start_methylation_call_procedure_paired_ends {
	my ($sequence_file_1,$sequence_file_2,$C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2,$pid) = @_;
	my ($dir_1,$filename_1);

	if ($sequence_file_1 =~ /\//){
		($dir_1,$filename_1) = $sequence_file_1 =~ m/(.*\/)(.*)$/;
	}
	else{
		$filename_1 = $sequence_file_1;
	}

	my ($dir_2,$filename_2);

	if  ($sequence_file_2 =~ /\//){
		($dir_2,$filename_2) = $sequence_file_2 =~ m/(.*\/)(.*)$/;
	}
	else{
		$filename_2 = $sequence_file_2;
	}

	### printing all alignments to a results file
	my $outfile = $filename_1;
	# warn "Outfile: $outfile\n";
	$outfile =~ s/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//; # attempting to remove fastq.gz etc to make filename a little shorter
	# warn "Outfile: $outfile\n";sleep(5);

	if ($prefix){
		$outfile = "$prefix.$outfile";
	}
	if ($bowtie2){ # SAM format is the default Bowtie 2 output
		$outfile =~ s/$/_bismark_bt2_pe.sam/;
	}
	else{ # SAM format is the default for HISAT2
		$outfile =~ s/$/_bismark_hisat2_pe.sam/;
	}

	if ($basename){ # Output file basename is set using the -B argument
		$outfile = "${basename}_pe.sam";
	}

	if ($ambig_bam){
		my $ambig_bam_out = $outfile;
		$ambig_bam_out =~ s/sam$/ambig.bam/;
		warn "Ambiguous BAM output: $ambig_bam_out\n";
		open (AMBIBAM,"| $samtools_path view -bSh 2>/dev/null - > $output_dir$ambig_bam_out") or die "Failed to write to $ambig_bam_out: $!\n";
	}

	$bam = 0 unless (defined $bam);

	if ($cram){ ### Samtools is installed, writing out CRAM directly. This qill require Samtools version 1.2 or higher!
		if ($multicore > 1){
			$outfile =~ s/sam$/bam/;
			open (OUT,"| $samtools_path view -bSh 2>/dev/null - > $output_dir$outfile") or die "Failed to write to $outfile: $!\n";
		}
		else{ # single-core mode
			$outfile =~ s/sam$/cram/;
			$final_output_filename = "${output_dir}${outfile}";
			open (OUT,"| $samtools_path view -h -C -T $cram_ref 2>/dev/null - > $output_dir$outfile") or die "Failed to write to CRAM file $outfile: $!\nPlease note that this option requires Samtools version 1.2 or higher!\n\n";
		}
	}
	elsif ($bam == 1){ ### Samtools is installed, writing out BAM directly
		$outfile =~ s/sam$/bam/;
		$final_output_filename = "${output_dir}${outfile}";
		open (OUT,"| $samtools_path view -bSh 2>/dev/null - > $output_dir$outfile") or die "Failed to write to $outfile: $!\n";
	}
	elsif($bam == 2){ ### no Samtools found on system. Using GZIP compression instead
		$outfile .= '.gz';
		open (OUT,"| gzip -c - > $output_dir$outfile") or die "Failed to write to $outfile: $!\n";
	}
	else{ # uncompressed ouput, default
		open (OUT,'>',"$output_dir$outfile") or die "Failed to write to $outfile: $!\n";
	}

	warn "\n>>> Writing bisulfite mapping results to $outfile <<<\n\n";
	sleep(1);

  
	### printing alignment and methylation call summary to a report file
	my $reportfile = $filename_1;
	$reportfile =~ s/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//; # attempting to remove fastq.gz etc to make filename a little shorter

	$reportfile =~ s/$/_bismark_vg_PE_report.txt/;

	if ($basename){ # Output file basename is set using the -B argument
		$reportfile = "${basename}_PE_report.txt";
	}

	open (REPORT,'>',"$output_dir$reportfile") or die "Failed to write to $reportfile: $!\n";
	print REPORT "Bismark report for: $sequence_file_1 and $sequence_file_2 (version: $bismark_version)\n";


	print REPORT "Bismark was run with HISAT2 against the bisulfite genome of $genome_folder with the specified options: $aligner_options\n";
	
	if ($directional){
		print REPORT "Option '--directional' specified (default mode): alignments to complementary strands (CTOT, CTOB) were ignored (i.e. not performed)\n\n";
	}
	else{
		print REPORT "Option '--non_directional' specified: alignments to all strands were being performed (OT, OB, CTOT, CTOB)\n\n";
	}

	unless ($sam_no_hd){
		generate_SAM_header();
	}


	### Input files are in FastQ format
	process_fastQ_files_for_paired_end_methylation_calls($sequence_file_1,$sequence_file_2,$C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2,$pid);

}

sub generate_SAM_header{

    print OUT "\@HD\tVN:1.0\tSO:unsorted\n";          # @HD = header, VN = version, SO = sort order

    #  Unordered printing of @SQ headers
    #  foreach my $chr (keys %chromosomes){
    #    my $length = length ($chromosomes{$chr});
    #    print "\@SQ\tSN:$chr\tLN:$length\n";
    #    print OUT "\@SQ\tSN:$chr\tLN:$length\n";        # @SQ = sequence, SN = seq name, LN = length
    #  }

    foreach my $chr (sort {$a<=>$b} keys %SQ_order){
		# warn "$chr\t$SQ_order{$chr}\n";
		my $length = length ($chromosomes{$SQ_order{$chr}});
		print OUT "\@SQ\tSN:$SQ_order{$chr}\tLN:$length\n"; # @SQ = sequence, SN = seq name, LN = length
		if ($ambig_bam){
			print AMBIBAM "\@SQ\tSN:$SQ_order{$chr}\tLN:$length\n";
		}
    }
    
    # 18 11 2015: Added @RG as a header line if --rg_tag or --rg_id/--rg_sample were set as well

    print OUT "\@PG\tID:Bismark\tVN:$bismark_version\tCL:\"bismark $command_line\"\n";        # @PG = program, ID = unique identifier, PN = program name name, VN = program version
}

sub process_fastQ_files_for_paired_end_methylation_calls{
	
	my ($sequence_file_1,$sequence_file_2,$C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2,$pid) = @_;
	### Processing the two Illumina sequence files; we need the actual sequence of both reads to compare them against the genomic sequence in order to
	### make a methylation call. The sequence identifier per definition needs to be same for a sequence pair used for paired-end alignments.
	### Now reading in the sequence files sequence by sequence and see if the current sequences produced a paired-end alignment to one (or both)
	### of the converted genomes (either C->T or G->A version)

	### gzipped version of the infiles
	if ($sequence_file_1 =~ /\.gz$/ and $sequence_file_2 =~ /\.gz$/){
		open (IN1,"gunzip -c $sequence_file_1 |") or die "Failed to open gunzip -c pipe to $sequence_file_1 $!\n";
		open (IN2,"gunzip -c $sequence_file_2 |") or die "Failed to open gunzip -c pipe to $sequence_file_2 $!\n";
	}
	else{
		open (IN1,$sequence_file_1) or die $!;
		open (IN2,$sequence_file_2) or die $!;
	}	

	my $count = 0;

	warn "\nReading in the sequence files $sequence_file_1 and $sequence_file_2\n";
	### Both files are required to have the exact same number of sequences, therefore we can process the sequences jointly one by one
	while (1) {
		# reading from the first input file
		my $identifier_1 = <IN1>;
		my $sequence_1 = <IN1>;
		my $ident_1 = <IN1>;         # not needed
		my $quality_value_1 = <IN1>; # not needed
		# reading from the second input file
		my $identifier_2 = <IN2>;
		my $sequence_2 = <IN2>;
		my $ident_2 = <IN2>;         # not needed
		my $quality_value_2 = <IN2>; # not needed
		last unless ($identifier_1 and $sequence_1 and $quality_value_1 and $identifier_2 and $sequence_2 and $quality_value_2);
	
		chomp $sequence_1;
		chomp $identifier_1;
		chomp $sequence_2;
		chomp $identifier_2;
		chomp $quality_value_1;
		chomp $quality_value_2;
	
		$identifier_1 = fix_IDs($identifier_1); # this is to avoid problems with truncated read ID when they contain white spaces
		$identifier_2 = fix_IDs($identifier_2);

		++$count;

		if ($skip){
			next unless ($count > $skip);
		}
		if ($upto){
			last if ($count > $upto);
		}

		$counting{sequences_count}++;
		if ($counting{sequences_count}%1000000==0) {
			warn "Processed $counting{sequences_count} sequence pairs so far\n";
		}

		my $orig_identifier_1 = $identifier_1;
		my $orig_identifier_2 = $identifier_2;

		$identifier_1 =~ s/^\@//;  # deletes the @ at the beginning of the FastQ ID

		my $return = check_results_paired_end (uc$sequence_1,uc$sequence_2,$identifier_1,$quality_value_1,$quality_value_2);
		
		unless ($return){
			$return = 0;
		}
	
	}

	warn "Processed $counting{sequences_count} sequences in total\n\n";

	close OUT or warn "Failed to close filehandle OUT: $!\n\n";

	print_final_analysis_report_paired_ends($C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2,$pid);

}

sub fix_IDs{
	my $id = shift;
	# warn "Got $id\n";
	if ($icpc){ # added 13 03 2019; see https://github.com/FelixKrueger/Bismark/issues/236
		$id =~ s/[ \t].*$//g; # truncating read ID at the first space or tab
	}
	else{
		$id =~ s/[ \t]+/_/g; # replace spaces or tabs with underscores
	}
	# warn "returning $id\n"; sleep(1);
	return $id;
}

sub check_results_paired_end{

    my ($sequence_1,$sequence_2,$identifier,$quality_value_1,$quality_value_2) = @_;
    
    ### quality values are not given for FastA files, so they are initialised with a Phred quality of 40
    unless ($quality_value_1){
		$quality_value_1 = 'I'x(length$sequence_1);
    }

    unless ($quality_value_2){
		$quality_value_2 = 'I'x(length$sequence_2);
    }
    # print "Read ID:$identifier\nLast ID [0]: $fhs[0]->{last_seq_id}\nLast ID [1]: $fhs[1]->{last_seq_id}\nLast ID [2]: $fhs[2]->{last_seq_id}\nLast ID [3]: $fhs[3]->{last_seq_id}\n\n"; sleep(1);
    
    my %alignments;
    my $alignment_ambiguous = 0;
    
    my $first_ambig_alignment_line1; # storing the first ambiguous alignment so it can be written out in case '--ambig_bam' was specified R1
    my $first_ambig_alignment_line2; # R2
    
    my $best_AS_so_far;   ## we need to keep a memory of the best alignment score so far
    my $amb_same_thread = 0;   ## if a read's primary and secondary alignments have the same alignment score we set this to true.

    ### reading from the Bowtie 2 output filehandles
    
    ### for paired end reads we are reporting alignments to the OT strand first (index 0), then the OB strand (index 3!!), similiar to the single end way.
    ### alignments to the complementary strands are reported afterwards (CTOT got index 1, and CTOB got index 2).
    ### This is needed so that alignments which either contain no single C or G or reads which contain only protected Cs are reported to the original strands (OT and OB)
    ### Before the complementary strands. Remember that it does not make any difference for the methylation calls, but it will matter if alignments to the complementary
    ### strands are not being reported when '--directional' is specified

    foreach my $index (0,3,1,2){
		### skipping this index if the last alignment has been set to undefined already (i.e. end of bowtie output)
		next unless ($fhs[$index]->{last_line_1} and $fhs[$index]->{last_line_2} and defined $fhs[$index]->{last_seq_id});

		### if the sequence pair we are currently looking at produced an alignment we are doing various things with it
		if ($fhs[$index]->{last_seq_id} eq $identifier) {
			my ($id_1,$flag_1,$mapped_chromosome_1,$position_1,$mapping_quality_1,$cigar_1,$bowtie_sequence_1,$qual_1) = (split (/\t/,$fhs[$index]->{last_line_1}))[0,1,2,3,4,5,9,10];
			my ($id_2,$flag_2,$mapped_chromosome_2,$position_2,$mapping_quality_2,$cigar_2,$bowtie_sequence_2,$qual_2) = (split (/\t/,$fhs[$index]->{last_line_2}))[0,1,2,3,4,5,9,10];
			# print "Index: $index\t$fhs[$index]->{last_line_1}\n";
			# print "Index: $index\t$fhs[$index]->{last_line_2}\n";
			# print join ("\t",$id_1,$flag_1,$mapped_chromosome_1,$position_1,$mapping_quality_1,$cigar_1,$bowtie_sequence_1,$qual_1),"\n";
			# print join ("\t",$id_2,$flag_2,$mapped_chromosome_2,$position_2,$mapping_quality_2,$cigar_2,$bowtie_sequence_2,$qual_2),"\n";
			$id_1 =~ s/\/1$//;
			$id_2 =~ s/\/2$//;

			### If a sequence has no reported alignments there will be a single output line per sequence with a bit-wise flag value of 77 for read 1 (1+4+8+64), or 141 for read 2 (1+4+8+128).
			### We can store the next alignment and move on to the next Bowtie 2 instance
			if ($flag_1 == 77 and $flag_2 == 141){

				## reading in the next alignment, which must be the next sequence
				my $newline_1 = $fhs[$index]->{fh}-> getline();
				my $newline_2 = $fhs[$index]->{fh}-> getline();
		
				if ($newline_1 and $newline_2){
					chomp $newline_1;
					chomp $newline_2;
					my ($seq_id_1) = split (/\t/,$newline_1);
					my ($seq_id_2) = split (/\t/,$newline_2);
					$seq_id_1 =~ s/\/1$//;
					$seq_id_2 =~ s/\/2$//;
					$fhs[$index]->{last_seq_id} = $seq_id_1;
					$fhs[$index]->{last_line_1} = $newline_1;
					$fhs[$index]->{last_line_2} = $newline_2;
		    
					#  print "current sequence ($identifier) did not map, reading in next sequence\n";
					#  print "$index\t$fhs[$index]->{last_seq_id}\n";
					#  print "$index\t$fhs[$index]->{last_line_1}\n";
					#  print "$index\t$fhs[$index]->{last_line_2}\n";
					next; # next instance
				}
				else{
					# assigning undef to last_seq_id and last_line and jumping to the next index (end of Bowtie 2 output)
					$fhs[$index]->{last_seq_id} = undef;
					$fhs[$index]->{last_line_1} = undef;
					$fhs[$index]->{last_line_2} = undef;
					next;
				}
			}
	    
			### If there are one or more proper alignments we can extract the chromosome number
			my ($chromosome_1,$chromosome_2);
			if ($mapped_chromosome_1 =~ s/_(CT|GA)_converted$//){
				$chromosome_1 = $mapped_chromosome_1;
			}
			else{
				die "Chromosome number extraction failed for $mapped_chromosome_1\n";
			}
			if ($mapped_chromosome_2 =~ s/_(CT|GA)_converted$//){
				$chromosome_2 = $mapped_chromosome_2;
			}
			else{
				die "Chromosome number extraction failed for $mapped_chromosome_2\n";
			}

			die "Paired-end alignments need to be on the same chromosome\n" unless ($chromosome_1 eq $chromosome_2);

			### We will use the optional fields to determine the best alignments. Later on we extract the number of mismatches and/or indels from the CIGAR string
			my ($alignment_score_1,$alignment_score_2,$second_best_1,$second_best_2,$MD_tag_1,$MD_tag_2);

			my @fields_1 = split (/\t/,$fhs[$index]->{last_line_1});
			my @fields_2 = split (/\t/,$fhs[$index]->{last_line_2});

			foreach (11..$#fields_1){
				if ($fields_1[$_] =~ /AS:i:(.*)/){
					$alignment_score_1 = $1;
				}
				elsif ($fields_1[$_] =~ /XS:i:(.*)/){
					$second_best_1 = $1;
				}
				elsif ($fields_1[$_] =~ /MD:Z:(.*)/){
					$MD_tag_1 = $1;
				}
			}

			foreach (11..$#fields_2){
				if ($fields_2[$_] =~ /AS:i:(.*)/){
					$alignment_score_2 = $1;
				}
				elsif ($fields_2[$_] =~ /MD:Z:(.*)/){
					$MD_tag_2 = $1;
				}
				else{
					if ($bowtie2){
						if ($fields_2[$_] =~ /XS:i:(.*)/){
							$second_best_2 = $1;
						}
					}	
					else{ # HISAT2 uses the ZS tag instead
						if($fields_2[$_] =~ /ZS:i:(.*)/){ 
							$second_best_2 = $1;
						}
					}
				}
			}

			die "Failed to extract alignment score 1 ($alignment_score_1) and MD tag ($MD_tag_1)!\nlast alignment 1: $fhs[$index]->{last_line_1}\nlast alignment 2: $fhs[$index]->{last_line_2}\n" unless (defined $alignment_score_1 and defined $MD_tag_1);
			die "Failed to extract alignment score 2 ($alignment_score_2) and MD tag ($MD_tag_2)!\nlast alignment 1: $fhs[$index]->{last_line_1}\nlast alignment 2: $fhs[$index]->{last_line_2}\n" unless (defined $alignment_score_2 and defined $MD_tag_2);

			# warn "First read 1 alignment score is: '$alignment_score_1'\n";
			# warn "First read 2 alignment score is: '$alignment_score_2'\n";
			# warn "XS/ZS tag 1 is: '$second_best_1'\n";
			# warn "XS/ZS tag 2 is: '$second_best_2'\n";
			# warn "MD tag 1 is: '$MD_tag_1'\n";
			# warn "MD tag 2 is: '$MD_tag_2'\n";

			### To decide whether a sequence pair has a unique best alignment we will look at the highest sum of alignment scores from both alignments
			my $sum_of_alignment_scores_1 = $alignment_score_1 + $alignment_score_2 ;
			# warn "sum of alignment scores: $sum_of_alignment_scores_1\n\n"; sleep(1);

			my $overwrite = 0; # If there are 2 alternative alignments to the same position, e.g. OT with 50 mismatches and CTOB with 0 mismatches, the CTOB one trumps the OT one.
			# introduced 13 April 2016 as a suggestion by Sylvain Foret, ANU Canberra

			if (!defined $best_AS_so_far){
				$overwrite = 1;
				$best_AS_so_far = $sum_of_alignment_scores_1;
				# warn "First alignment score, setting \$best_AS_so_far to $best_AS_so_far\n";
				if ($ambig_bam){ # also setting the first_ambig_alignment
					# Read 1
					$first_ambig_alignment_line1 = $fhs[$index]->{last_line_1};
					$first_ambig_alignment_line1 =~ s/_(CT|GA)_converted//;
					# Read 2
					$first_ambig_alignment_line2 = $fhs[$index]->{last_line_2};
					$first_ambig_alignment_line2 =~ s/_(CT|GA)_converted//;
					# warn "$first_ambig_alignment_line1\n$first_ambig_alignment_line2\n\n"; sleep(1);
				}
			}
			else{
				if ($sum_of_alignment_scores_1 >= $best_AS_so_far){ # AS are generally negative with a maximum of 0
					# 19 07 2016 Changed to >= so that equally good alignments to different positions get added as well. Ambiguous alignments are identified and removed later.
					$overwrite = 1;

					# resetting the ambiguous within thread memory (if applicable at all) only if the current alignment is really better than the previous one.
					# 22 07 2016: ambiguous score within same thread only resets if the current alignment is really better than the previous one
					if ($sum_of_alignment_scores_1 > $best_AS_so_far){
						# warn "Resetting amb within thread value to 0\n";
						$amb_same_thread = 0;
		    
					}
					$best_AS_so_far = $sum_of_alignment_scores_1; # moved this down so that $amb_within_thread gets a chance to be reset
					# warn "Found better or equal sum of alignment scores ($sum_of_alignment_scores_1), setting \$best_AS_so_far to $best_AS_so_far\n";
				}
				else{
					# warn "current alignment (AS $sum_of_alignment_scores) isn't better than the best so far ($best_AS_so_far). Not changing anything\n";
				}		
			}
	    
			# If either of the reads has a second best alignment but the other one doesn't we assign a value of the best alignment, i.e. the AS score
			if (defined $second_best_1 or defined $second_best_2){    
				unless (defined $second_best_1){
					$second_best_1 = $alignment_score_1;
				}
				unless (defined $second_best_2){
					$second_best_2 = $alignment_score_2;
				}
				# warn "\n\n#############################\n\nXS Read 1: $second_best_1\nXS Read 2: $second_best_2\n\n##########################\n\n";
			}
	    
			if (defined $second_best_1 and defined $second_best_2){
				my $sum_of_alignment_scores_second_best = $second_best_1 + $second_best_2;
				# warn "Second best alignment_score_1 is: '$second_best_1'\n";
				# warn "Second best alignment_score_2 is: '$second_best_2'\n";
				# warn "Second best alignment sum of alignment scores is: '$sum_of_alignment_scores_second_best'\n";
	    
				# If the first alignment score for the first read pair is the same as the alignment score of the second best hit we we keep a memory of this
				if ($sum_of_alignment_scores_1 == $sum_of_alignment_scores_second_best){
					# checking to see if this read pair produced the best alignment
					if ($sum_of_alignment_scores_1 == $best_AS_so_far){  # yes this is the best read pair so far, either within the thread or between threads, however it is ambiguous
						#warn "Read pair is ambiguous within the same thread, or otherwise as good as the best one so far. Setting \$amb_same_thread to 1 for currently best AS: $best_AS_so_far\n";
						$amb_same_thread = 1;
					}
					else{
						# warn "This read pair has a worse alignment score than the best alignment so far and will be ignored even though it is ambiguous in itself\n";
					}

					### if there is a better alignment later on -> fine. If not, the read will get booted altogether one way or another

					## need to read and discard all additional ambiguous reads until we reach the next sequence
					until ($fhs[$index]->{last_seq_id} ne $identifier){
						my $newline_1 = $fhs[$index]->{fh}-> getline();
						my $newline_2 = $fhs[$index]->{fh}-> getline();
						if ($newline_1 and $newline_2){
							chomp $newline_1;
							chomp $newline_2;
							my ($seq_id_1) = split (/\t/,$newline_1);
							my ($seq_id_2) = split (/\t/,$newline_2);
							$seq_id_1 =~ s/\/1$//;
							$seq_id_2 =~ s/\/2$//;
							# print "New Seq IDs:\t$seq_id_1\t$seq_id_2\n";

							$fhs[$index]->{last_seq_id} = $seq_id_1;
							$fhs[$index]->{last_line_1} = $newline_1;
							$fhs[$index]->{last_line_2} = $newline_2;
						}
						else{
							# assigning undef to last_seq_id and last_line and jumping to the next index (end of Bowtie 2 output)
							$fhs[$index]->{last_seq_id} = undef;
							$fhs[$index]->{last_line_1} = undef;
							$fhs[$index]->{last_line_2} = undef;
							last; # break free if the end of the alignment output was reached
						}
					}
				}
				#  if ($fhs[$index]->{last_seq_id}){
				#    warn "Index: $index\tThis Seq-ID is $identifier, skipped all ambiguous sequences until the next ID which is: $fhs[$index]->{last_seq_id}\n";
				#  }
			    else{ # the next best alignment has a lower alignment score than the current read, so we can safely store the current alignment

					my $alignment_location;
					if ($position_1 <= $position_2){
						$alignment_location = join(":",$chromosome_1,$position_1,$position_2);
					}
					elsif($position_2 < $position_1){
						$alignment_location = join(":",$chromosome_1,$position_2,$position_1);
					}

					### If a sequence aligns to exactly the same location twice the sequence does either not contain any C or G, or all the Cs (or Gs on the reverse
					### strand) were methylated and therefore protected. Alternatively it will align better in one condition than in the other. In any case, it is not needed to overwrite
					### the same positional entry with a second entry for the same location, as the genomic sequence extraction and methylation call would not be affected by this. The only
					### thing which would change is the index number for the found alignment). We will continue to assign these alignments to the first indexes 0 and 3, i.e. OT and OB

					if ($overwrite){ # see comment above at "my $overwrite = ..."
						#unless (exists $alignments{$alignment_location}){
						$alignments{$alignment_location}->{seq_id} = $id_1;
						$alignments{$alignment_location}->{alignment_score_1} = $alignment_score_1;
						$alignments{$alignment_location}->{alignment_score_2} = $alignment_score_2;
						$alignments{$alignment_location}->{sum_of_alignment_scores} = $sum_of_alignment_scores_1;
						$alignments{$alignment_location}->{sum_of_alignment_scores_second_best} = $sum_of_alignment_scores_second_best;
						$alignments{$alignment_location}->{bowtie_sequence_1} = $bowtie_sequence_1;
						$alignments{$alignment_location}->{bowtie_sequence_2} = $bowtie_sequence_2;
						$alignments{$alignment_location}->{index} = $index;
						$alignments{$alignment_location}->{chromosome} = $chromosome_1; # either is fine
						$alignments{$alignment_location}->{position_1} = $position_1;
						$alignments{$alignment_location}->{position_2} = $position_2;
						$alignments{$alignment_location}->{mismatch_info_1} = $MD_tag_1;
						$alignments{$alignment_location}->{mismatch_info_2} = $MD_tag_2;
						$alignments{$alignment_location}->{CIGAR_1} = $cigar_1;
						$alignments{$alignment_location}->{CIGAR_2} = $cigar_2;
						$alignments{$alignment_location}->{flag_1} = $flag_1;
						$alignments{$alignment_location}->{flag_2} = $flag_2;
						# warn "added best of several alignments to \%alignments hash\n";
					}

					### now reading and discarding all (inferior) alignments of this read pair until we hit the next sequence
					until ($fhs[$index]->{last_seq_id} ne $identifier){
						my $newline_1 = $fhs[$index]->{fh}-> getline();
						my $newline_2 = $fhs[$index]->{fh}-> getline();
						if ($newline_1 and $newline_2){
							chomp $newline_1;
							chomp $newline_2;
							my ($seq_id_1) = split (/\t/,$newline_1);
							my ($seq_id_2) = split (/\t/,$newline_2);
							$seq_id_1 =~ s/\/1$//;
							$seq_id_2 =~ s/\/2$//;
							# print "New Seq IDs:\t$seq_id_1\t$seq_id_2\n";

							$fhs[$index]->{last_seq_id} = $seq_id_1;
							$fhs[$index]->{last_line_1} = $newline_1;
							$fhs[$index]->{last_line_2} = $newline_2;
						}
						else{
							# assigning undef to last_seq_id and last_line_1 and _2 and jumping to the next index (end of Bowtie 2 output)
							$fhs[$index]->{last_seq_id} = undef;
							$fhs[$index]->{last_line_1} = undef;
							$fhs[$index]->{last_line_2} = undef;
							last; # break free if the end of the alignment output was reached
						}
					}
			# if($fhs[$index]->{last_seq_id}){
			#   warn "Index: $index\tThis Seq-ID is $identifier, skipped all other alignments until the next ID was reached which is: $fhs[$index]->{last_seq_id}\n";
			# }
				}
			}	
			else{ # there is no second best hit, so we can just store this one and read in the next sequence

				my $alignment_location = join(":",$chromosome_1,$position_1,$position_2);
					# print "$alignment_location\n";
				### If a sequence aligns to exactly the same location with a perfect match twice the sequence does either not contain any C or G, or all the Cs (or Gs on the reverse
				### strand) were methylated and therefore protected. Alternatively it will align better in one condition than in the other. In any case, it is not needed to overwrite
				### the same positional entry with a second entry for the same location, as the genomic sequence extraction and methylation call would not be affected by this. The only
				### thing which would change is the index number for the found alignment). We will continue to assign these alignments to the first indexes 0 and 3, i.e. OT and OB

				#unless (exists $alignments{$alignment_location}){ # see comment above at my $overwrite = ...
				if ($overwrite){
					$alignments{$alignment_location}->{seq_id} = $id_1;
					$alignments{$alignment_location}->{alignment_score_1} = $alignment_score_1;
					$alignments{$alignment_location}->{alignment_score_2} = $alignment_score_2;
					$alignments{$alignment_location}->{sum_of_alignment_scores} = $sum_of_alignment_scores_1;
					$alignments{$alignment_location}->{sum_of_alignment_scores_second_best} = undef;
					$alignments{$alignment_location}->{bowtie_sequence_1} = $bowtie_sequence_1;
					$alignments{$alignment_location}->{bowtie_sequence_2} = $bowtie_sequence_2;
					$alignments{$alignment_location}->{index} = $index;
					$alignments{$alignment_location}->{chromosome} = $chromosome_1; # either is fine
					$alignments{$alignment_location}->{position_1} = $position_1;
					$alignments{$alignment_location}->{position_2} = $position_2;
					$alignments{$alignment_location}->{mismatch_info_1} = $MD_tag_1;
					$alignments{$alignment_location}->{mismatch_info_2} = $MD_tag_2;
					$alignments{$alignment_location}->{CIGAR_1} = $cigar_1;
					$alignments{$alignment_location}->{CIGAR_2} = $cigar_2;
					$alignments{$alignment_location}->{flag_1} = $flag_1;
					$alignments{$alignment_location}->{flag_2} = $flag_2;
					# warn "added unique alignment to \%alignments hash\n";
				}

				# Now reading and storing the next read pair
				until ($fhs[$index]->{last_seq_id} ne $identifier){
					my $newline_1 = $fhs[$index]->{fh}-> getline();
					my $newline_2 = $fhs[$index]->{fh}-> getline();
					if ($newline_1 and $newline_2){
						chomp $newline_1;
						chomp $newline_2;
						# print "$newline_1\n";
						# print "$newline_2\n";
						my ($seq_id_1) = split (/\t/,$newline_1);
						my ($seq_id_2) = split (/\t/,$newline_2);
						$seq_id_1 =~ s/\/1$//;
						$seq_id_2 =~ s/\/2$//;
						# print "New Seq IDs:\t$seq_id_1\t$seq_id_2\n";

						$fhs[$index]->{last_seq_id} = $seq_id_1;
						$fhs[$index]->{last_line_1} = $newline_1;
						$fhs[$index]->{last_line_2} = $newline_2;
					}
					else{
						# assigning undef to last_seq_id and last_line_1 and _2 and jumping to the next index (end of Bowtie 2 output)
						$fhs[$index]->{last_seq_id} = undef;
						$fhs[$index]->{last_line_1} = undef;
						$fhs[$index]->{last_line_2} = undef;
						last; # break free if the end of the alignment output was reached
					}
				}
			}
		}
	}

	### If there were several equally good alignments for the best alignment score we will boot the read
	if ($amb_same_thread){
		# warn "\$alignment_ambiguous now: $alignment_ambiguous\n";
		$alignment_ambiguous = 1;
		# warn "\$alignment_ambiguous now: $alignment_ambiguous\n";
	}
	else{
		# warn "alignment won't be considered ambiguous. This time....\n";
	}


	### if the read produced several ambiguous alignments for a single instance of Bowtie 2 we can return already now. If --ambiguous was specified the read sequence will be printed out in FastQ format
	if ($alignment_ambiguous == 1){
		$counting{unsuitable_sequence_count}++;
		### report that the sequence pair has multiple hits with bitwise flag 256. We can print the sequence to the result file straight away and skip everything else
		#  my $ambiguous_read_1 = join("\t",$identifier.'/1','256','*','0','0','*','*','0','0',$sequence_1,$quality_value_1);
		#  my $ambiguous_read_2 = join("\t",$identifier.'/2','256','*','0','0','*','*','0','0',$sequence_2,$quality_value_2);
		#  print "$ambiguous_read_1\n";
		#  print "$ambiguous_read_2\n";

		if ($ambig_bam){
			# warn "Sequence is ambiguous, printing out to ambiguous BAM file:\n";
			# replacing the first /1\t in the ID of R1
			# warn "Was\n$first_ambig_alignment_line1\n$first_ambig_alignment_line2\n";
			$first_ambig_alignment_line1 =~ s/\/1\t/\t/;
			$first_ambig_alignment_line2 =~ s/\/2\t/\t/;
			# warn "Is:\n$first_ambig_alignment_line1\n$first_ambig_alignment_line2\n\n";

			print AMBIBAM "$first_ambig_alignment_line1\n$first_ambig_alignment_line2\n";
			# print "$first_ambig_alignment_line1\n$first_ambig_alignment_line2\n";
		}

		if ($ambiguous){
			return 2; # => exits to next sequence pair, and prints it out to _ambiguous_reads_1.txt and _ambiguous_reads_2.txt if '--ambiguous' was specified
		}
		elsif ($unmapped){
			return 1; # => exits to next sequence pair, and prints it out to _unmapped_reads_1.txt and _unmapped_reads_2.txt if '--unmapped' but not '--ambiguous' was specified
		}
		else{
			return 0;
		}
	}

	### if no alignment was found for a certain sequence at all we continue with the next sequence in the sequence file
	unless (%alignments){
		$counting{no_single_alignment_found}++;

		# my $unmapped_read_1 = join("\t",$identifier.'/1','77','*','0','0','*','*','0','0',$sequence_1,$quality_value_1);
		# my $unmapped_read_2 = join("\t",$identifier.'/2','141','*','0','0','*','*','0','0',$sequence_2,$quality_value_2);
		# print "$unmapped_read_1\n";
		# print "$unmapped_read_2\n";
		if ($unmapped){
			return 1; # => exits to next sequence pair, and prints it out to _unmapped_reads_1.txt and _unmapped_read_2.txt if '--unmapped' was specified
		}
		else{
			return 0;
		}
	}

	#######################################################################################################################################################
	
	### If the sequence pair was not rejected so far we are now looking if there is a unique best alignment among all alignment instances. If there is only one
	### single best position we are going to store the alignment information in the $meth_call variable. If there are multiple hits with the same (highest)
	### alignment score we are discarding the sequence pair altogether.
	### For end-to-end alignments the maximum alignment score is 0, each mismatch receives a penalty of 6, and each gap receives penalties for opening (5)
	### and extending (3 per bp) the gap.

	#######################################################################################################################################################

	### Declaring an empty hash reference which will store all information we need for the methylation call
	my $methylation_call_params; # hash reference
	my $sequence_pair_fails = 0; # using $sequence_pair_fails as a 'memory' if a sequence could not be aligned uniquely (set to 1 then)

	  ### print contents of %alignments for debugging
	  ##  if (scalar keys %alignments >= 1){
	  #     print "\n******\n";
	  #     foreach my $alignment_location (sort {$a cmp $b} keys %alignments){
	  #       print "Loc:  $alignment_location\n";
	  #       print "ID:      $alignments{$alignment_location}->{seq_id}\n";
	  #       print "AS_1:    $alignments{$alignment_location}->{alignment_score_1}\n";
	  #       print "AS_2:    $alignments{$alignment_location}->{alignment_score_2}\n";
	  #       print "Seq_1:   $alignments{$alignment_location}->{bowtie_sequence_1}\n";
	  #       print "Seq_2:   $alignments{$alignment_location}->{bowtie_sequence_2}\n";
	  #       print "Index    $alignments{$alignment_location}->{index}\n";
	  #       print "Chr:     $alignments{$alignment_location}->{chromosome}\n";
	  #       print "Pos_1:   $alignments{$alignment_location}->{position_1}\n";
	  #       print "Pos_2:   $alignments{$alignment_location}->{position_2}\n";
	  #       print "CIGAR_1: $alignments{$alignment_location}->{CIGAR_1}\n";
	  #       print "CIGAR_2: $alignments{$alignment_location}->{CIGAR_2}\n";
	  #       print "MD_1:    $alignments{$alignment_location}->{mismatch_info_1}\n";
	  #       print "MD_2:    $alignments{$alignment_location}->{mismatch_info_2}\n";
	  #       print "Flag 1:  $alignments{$alignment_location}->{flag_1}\n";
	  #       print "Flag 2:  $alignments{$alignment_location}->{flag_2}\n";
	  #    }
	  #    print "\n******\n";
	  #  }

	## if there is only 1 entry in the %alignments hash we accept it as the best alignment
	if (scalar keys %alignments == 1){
		for my $unique_best_alignment (keys %alignments){
			$methylation_call_params->{$identifier}->{bowtie_sequence_1} = $alignments{$unique_best_alignment}->{bowtie_sequence_1};
			$methylation_call_params->{$identifier}->{bowtie_sequence_2} = $alignments{$unique_best_alignment}->{bowtie_sequence_2};
			$methylation_call_params->{$identifier}->{chromosome}        = $alignments{$unique_best_alignment}->{chromosome};
			$methylation_call_params->{$identifier}->{position_1}        = $alignments{$unique_best_alignment}->{position_1};
			$methylation_call_params->{$identifier}->{position_2}        = $alignments{$unique_best_alignment}->{position_2};
			$methylation_call_params->{$identifier}->{index}             = $alignments{$unique_best_alignment}->{index};
			$methylation_call_params->{$identifier}->{alignment_score_1} = $alignments{$unique_best_alignment}->{alignment_score_1};
			$methylation_call_params->{$identifier}->{alignment_score_2} = $alignments{$unique_best_alignment}->{alignment_score_2};
			$methylation_call_params->{$identifier}->{sum_of_alignment_scores} = $alignments{$unique_best_alignment}->{sum_of_alignment_scores};
			$methylation_call_params->{$identifier}->{sum_of_alignment_scores_second_best} = $alignments{$unique_best_alignment}->{sum_of_alignment_scores_second_best};
			$methylation_call_params->{$identifier}->{mismatch_info_1}   = $alignments{$unique_best_alignment}->{mismatch_info_1};
			$methylation_call_params->{$identifier}->{mismatch_info_2}   = $alignments{$unique_best_alignment}->{mismatch_info_2};
			$methylation_call_params->{$identifier}->{CIGAR_1}           = $alignments{$unique_best_alignment}->{CIGAR_1};
			$methylation_call_params->{$identifier}->{CIGAR_2}           = $alignments{$unique_best_alignment}->{CIGAR_2};
			$methylation_call_params->{$identifier}->{flag_1}            = $alignments{$unique_best_alignment}->{flag_1};
			$methylation_call_params->{$identifier}->{flag_2}            = $alignments{$unique_best_alignment}->{flag_2};
		}
	}

	### otherwise we are going to find out if there is a best match among the multiple alignments, or whether there are 2 or more equally good alignments (in which case
	### we boot the sequence pair altogether)
	elsif (scalar keys %alignments >= 2  and scalar keys %alignments <= 4){
		my $best_sum_of_alignment_scores;
		my $best_alignment_location;
		foreach my $alignment_location (sort {$alignments{$b}->{sum_of_alignment_scores} <=> $alignments{$a}->{sum_of_alignment_scores}} keys %alignments){

			# warn "$alignments{$alignment_location}->{sum_of_alignment_scores}\n"; sleep(1);

			unless (defined $best_sum_of_alignment_scores){
				$best_sum_of_alignment_scores = $alignments{$alignment_location}->{sum_of_alignment_scores};
				$best_alignment_location = $alignment_location;
				# print "setting best alignment score to: $best_sum_of_alignment_scores\n";
			}
			else{
				### if the second best alignment has the same sum of alignment scores as the first one, the sequence pair will get booted
				if ($alignments{$alignment_location}->{sum_of_alignment_scores} == $best_sum_of_alignment_scores){
					# warn "Same sum of alignment scores for 2 different alignments, the sequence pair will get booted!\n";
					$sequence_pair_fails = 1;
					last; # exiting since we know that the sequence has ambiguous alignments
				}
				### else we are going to store the best alignment for further processing
				else{
					$methylation_call_params->{$identifier}->{bowtie_sequence_1} = $alignments{$best_alignment_location}->{bowtie_sequence_1};
					$methylation_call_params->{$identifier}->{bowtie_sequence_2} = $alignments{$best_alignment_location}->{bowtie_sequence_2};
					$methylation_call_params->{$identifier}->{chromosome}        = $alignments{$best_alignment_location}->{chromosome};
					$methylation_call_params->{$identifier}->{position_1}        = $alignments{$best_alignment_location}->{position_1};
					$methylation_call_params->{$identifier}->{position_2}        = $alignments{$best_alignment_location}->{position_2};
					$methylation_call_params->{$identifier}->{index}             = $alignments{$best_alignment_location}->{index};
					$methylation_call_params->{$identifier}->{alignment_score_1} = $alignments{$best_alignment_location}->{alignment_score_1};
					$methylation_call_params->{$identifier}->{alignment_score_2} = $alignments{$best_alignment_location}->{alignment_score_2};
					$methylation_call_params->{$identifier}->{sum_of_alignment_scores} = $alignments{$best_alignment_location}->{sum_of_alignment_scores};
					$methylation_call_params->{$identifier}->{mismatch_info_1}   = $alignments{$best_alignment_location}->{mismatch_info_1};
					$methylation_call_params->{$identifier}->{mismatch_info_2}   = $alignments{$best_alignment_location}->{mismatch_info_2};
					$methylation_call_params->{$identifier}->{CIGAR_1}           = $alignments{$best_alignment_location}->{CIGAR_1};
					$methylation_call_params->{$identifier}->{CIGAR_2}           = $alignments{$best_alignment_location}->{CIGAR_2};
					$methylation_call_params->{$identifier}->{flag_1}            = $alignments{$best_alignment_location}->{flag_1};
					$methylation_call_params->{$identifier}->{flag_2}            = $alignments{$best_alignment_location}->{flag_2};

					if (defined $alignments{$best_alignment_location}->{sum_of_alignment_scores_second_best} and ( $alignments{$best_alignment_location}->{sum_of_alignment_scores_second_best} > $alignments{$alignment_location}->{sum_of_alignment_scores} )) {
						$methylation_call_params->{$identifier}->{sum_of_alignment_scores_second_best} = $alignments{$best_alignment_location}->{sum_of_alignment_scores_second_best};
					}
					else {
						$methylation_call_params->{$identifier}->{sum_of_alignment_scores_second_best} = $alignments{$alignment_location}->{sum_of_alignment_scores};
					}

					last; # exiting since the sequence produced a unique best alignment
				}
			}
		}
	}
	else{
		die "There are too many potential hits for this sequence pair (1-4 expected, but found: '",scalar keys %alignments,"')\n";;
	}

	### skipping the sequence completely if there were multiple alignments with the same best sum of alignment scores at different positions
	if ($sequence_pair_fails == 1){
		$counting{unsuitable_sequence_count}++;

		### report that the sequence has multiple hits with bitwise flag 256. We can print the sequence to the result file straight away and skip everything else
		# my $ambiguous_read_1 = join("\t",$identifier.'/1','256','*','0','0','*','*','0','0',$sequence_1,$quality_value_1);
		# my $ambiguous_read_2 = join("\t",$identifier.'/2','256','*','0','0','*','*','0','0',$sequence_2,$quality_value_2);
		# warn "$ambiguous_read_1\n";
		# warn "$ambiguous_read_2\n";

		if ($ambiguous){
			return 2; # => exits to next sequence pair, and prints it out (in FastQ format) to _ambiguous_reads_1.txt and _ambiguous_reads_2.txt if '--ambiguous' was specified
		}
		elsif ($unmapped){
			return 1; # => exits to next sequence pair, and prints it out (in FastQ format) to _unmapped_reads_1.txt and _unmapped_reads_2.txt if '--unmapped' but not '--ambiguous' was specified
		}
		else{
			return 0; # => exits to next sequence pair (default)
		}
	}

	### --DIRECTIONAL
	### If the option --directional has been specified the user wants to consider only alignments to the original top strand or the original bottom strand. We will therefore
	### discard all alignments to strands complementary to the original strands, as they should not exist in reality due to the library preparation protocol
	if ($directional){
		if ( ($methylation_call_params->{$identifier}->{index} == 1) or ($methylation_call_params->{$identifier}->{index} == 2) ){
			#    warn "Alignment rejected! (index was: $methylation_call_params->{$identifier}->{index})\n";
			$counting{alignments_rejected_count}++;
			return 0;
		}
	}

	### If the sequence pair has not been rejected so far it does have a unique best alignment
	$counting{unique_best_alignment_count}++;
	extract_corresponding_genomic_sequence_paired_end($identifier,$methylation_call_params);

	### check to see if the genomic sequences we extracted has the same length as the observed sequences +2, and only then we perform the methylation call
	if (length($methylation_call_params->{$identifier}->{unmodified_genomic_sequence_1}) != length($sequence_1)+2){
		warn "Chromosomal sequence could not be extracted for\t$identifier\t$methylation_call_params->{$identifier}->{chromosome}\t$methylation_call_params->{$identifier}->{position_1}\n";
		$counting{genomic_sequence_could_not_be_extracted_count}++;
		return 0;
	}
	if (length($methylation_call_params->{$identifier}->{unmodified_genomic_sequence_2}) != length($sequence_2)+2){
		warn "Chromosomal sequence could not be extracted for\t$identifier\t$methylation_call_params->{$identifier}->{chromosome}\t$methylation_call_params->{$identifier}->{position_2}\n";
		$counting{genomic_sequence_could_not_be_extracted_count}++;
		return 0;
	}

	### Compute MAPQ value
	$methylation_call_params->{$identifier}->{mapq} = calc_mapq (length($sequence_1), length($sequence_2),
                                                                           $methylation_call_params->{$identifier}->{sum_of_alignment_scores},
                                                                           $methylation_call_params->{$identifier}->{sum_of_alignment_scores_second_best});


	### now we are set to perform the actual methylation call
	if ($slam){
		$methylation_call_params->{$identifier}->{methylation_call_1} = methylation_call_slam($identifier,$sequence_1,$methylation_call_params->{$identifier}->{unmodified_genomic_sequence_1},$methylation_call_params->{$identifier}->{read_conversion_1});
		$methylation_call_params->{$identifier}->{methylation_call_2} = methylation_call_slam($identifier,$sequence_2,$methylation_call_params->{$identifier}->{unmodified_genomic_sequence_2},$methylation_call_params->{$identifier}->{read_conversion_2});
	}
	else{
		$methylation_call_params->{$identifier}->{methylation_call_1} = methylation_call($identifier,$sequence_1,$methylation_call_params->{$identifier}->{unmodified_genomic_sequence_1},$methylation_call_params->{$identifier}->{read_conversion_1});
		$methylation_call_params->{$identifier}->{methylation_call_2} = methylation_call($identifier,$sequence_2,$methylation_call_params->{$identifier}->{unmodified_genomic_sequence_2},$methylation_call_params->{$identifier}->{read_conversion_2});
	}
	# warn "$methylation_call_params->{$identifier}->{read_conversion_2}\n";
	# warn "  $sequence_2\n";
	# warn "$methylation_call_params->{$identifier}->{unmodified_genomic_sequence_2}\n";
	# warn "  $methylation_call_params->{$identifier}->{methylation_call_2}\n";

	print_bisulfite_mapping_results_paired_ends($identifier,$sequence_1,$sequence_2,$methylation_call_params,$quality_value_1,$quality_value_2);
	return 0; ## otherwise 1 will be returned by default, which would print the sequence pair to unmapped_1 and _2
}

### EXTRACT GENOMIC SEQUENCE | PAIRED-END

sub extract_corresponding_genomic_sequence_paired_end{
	my ($sequence_identifier,$methylation_call_params) = @_;
	### A bisulfite sequence pair for 1 location in the genome can theoretically be on any of the 4 possible converted strands. We are also giving the
	### sequence a 'memory' of the conversion we are expecting which we will need later for the methylation call
	# $verbose = 1;
	my $cigar_1 = $methylation_call_params->{$sequence_identifier}->{CIGAR_1};
	my $cigar_2 = $methylation_call_params->{$sequence_identifier}->{CIGAR_2};
	my $flag_1 =  $methylation_call_params->{$sequence_identifier}->{flag_1};
	my $flag_2 =  $methylation_call_params->{$sequence_identifier}->{flag_2};

	my $contains_deletion_1 = 0;
	my $contains_deletion_2 = 0;
	if ($cigar_1 =~ /D/){
		$contains_deletion_1 = 1;
		if ($verbose){ warn "$cigar_1\n$methylation_call_params->{$sequence_identifier}->{mismatch_info_1}\n";}
	}
	if ($cigar_2 =~ /D/){
		$contains_deletion_2 = 1;
		if ($verbose){ warn "$cigar_2\n$methylation_call_params->{$sequence_identifier}->{mismatch_info_2}\n";}
	}

	
	# warn "$cigar_1\t$cigar_2\t$flag_1\t$flag_2\n";
	### We are now extracting the corresponding genomic sequence, +2 extra bases at the end (or start) so that we can also make a CpG methylation call and
	### in addition make differential calls for Cs in CHG or CHH context if the C happens to be at the last (or first)  position of the actually observed sequence

	### the alignment_strand information is needed to determine which strand of the genomic sequence we are comparing the read against,
	### the read_conversion information is needed to know whether we are looking for C->T or G->A substitutions
	my $alignment_read_1;
	my $alignment_read_2;
	my $read_conversion_info_1;
	my $read_conversion_info_2;
	my $genome_conversion;

	### Now extracting the same sequence from the mouse genomic sequence, +2 extra bases at one of the ends so that we can also make a CpG, CHG or CHH methylation call
	### if the C happens to be at the last position of the actually observed sequence
	my $non_bisulfite_sequence_1 = '';
	my $non_bisulfite_sequence_2 = '';
	my $genomic_seq_for_MD_tag_1 = ''; # this sequence contains potential deletions in the genome as well so that we can generate a proper MD tag for the SAM output
	my $genomic_seq_for_MD_tag_2 = '';

	### Positions in SAM format are 1 based, so we need to subract 1 when getting substrings
	my $pos_1 = $methylation_call_params->{$sequence_identifier}->{position_1}-1;
	my $pos_2 = $methylation_call_params->{$sequence_identifier}->{position_2}-1;

	# parsing CIGAR 1 string
	my @len_1 = split (/\D+/,$cigar_1); # storing the length per operation
	my @ops_1 = split (/\d+/,$cigar_1); # storing the operation
	shift @ops_1; # remove the empty first element
	die "CIGAR 1 string contained a non-matching number of lengths and operations\n" unless (scalar @len_1 == scalar @ops_1);
	# parsing CIGAR 2 string
	my @len_2 = split (/\D+/,$cigar_2); # storing the length per operation
	my @ops_2 = split (/\d+/,$cigar_2); # storing the operation
	shift @ops_2; # remove the empty first element
	die "CIGAR 2 string contained a non-matching number of lengths and operations\n" unless (scalar @len_2 == scalar @ops_2);

	my $indels_1 = 0; # adding these to the hemming distance value (needed for the NM field in the final SAM output
	my $indels_2 = 0;

	### Extracting read 1 genomic sequence ###

	# extracting 2 additional bp at the 5' end (read 1)
	if ( ($methylation_call_params->{$sequence_identifier}->{index} == 1) or ($methylation_call_params->{$sequence_identifier}->{index} == 3) ){
		# checking if the substring will be valid or if we can't extract the sequence because we are right at the edge of a chromosome
		unless ( ($pos_1-2) > 0){# exiting with en empty genomic sequence otherwise
			$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_1} = $non_bisulfite_sequence_1;
			$methylation_call_params->{$sequence_identifier}->{genomic_seq_for_MD_tag_1} = $genomic_seq_for_MD_tag_1;
			return;
		}
		$non_bisulfite_sequence_1 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_1-2,2);
	}

	foreach (0..$#len_1){
		if ($ops_1[$_] eq 'M'){
			# extracting genomic sequence
			$non_bisulfite_sequence_1 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_1,$len_1[$_]);
			if ($contains_deletion_1){
				$genomic_seq_for_MD_tag_1 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_1,$len_1[$_]);
			}
			#   warn "$non_bisulfite_sequence_1\n";
			# adjusting position
			$pos_1 += $len_1[$_];
		}
		elsif ($ops_1[$_] eq 'I'){ # insertion in the read sequence
			# we simply add padding Xs instead of finding genomic sequence. This will not be used to infer methylation calls, and we can later ignore it for the generation of the MD;Z: tag
			$non_bisulfite_sequence_1 .= 'X' x $len_1[$_];
			if ($contains_deletion_1){
				$genomic_seq_for_MD_tag_1 .= 'X' x $len_1[$_];
			}
			# warn "$non_bisulfite_sequence_1\n";
			# position doesn't need adjusting

			### 03 06 2014: In fact we don't need to add anything to the hemming distance for insertions since we use padding Xs which will fail a base by base comparison in hemming_dist()
			# indels_1 += $len_1[$_]; # adding to $indels_1 to determine the hemming distance (= single base mismatches, insertions or deletions) for the SAM output
		}
		elsif ($ops_1[$_] eq 'S'){ # soft-clipped read sequence
			# we simply add padding Xs instead of finding genomic sequence. These will not be used to infer methylation calls,
			# and we can later ignore them better during the generation of the MD:Z: tag
			$non_bisulfite_sequence_1 .= 'X' x $len_1[$_];
			if ($contains_deletion_1){
				$genomic_seq_for_MD_tag_1 .= 'X' x $len_1[$_];
			}
			# warn "Soft-clipped Read 1 sequence! $sequence_identifier\n";
			# warn "$non_bisulfite_sequence_1\n";
			# position doesn't need adjusting

			### We don't need to add anything to the hemming distance for soft-clipped bases 
			# since we use padding Xs which will fail a base by base comparison in hemming_dist()
		}
		elsif ($ops_1[$_] eq 'D'){ # deletion in the read sequence
			# we do not add any genomic sequence but only adjust the position
			# we do however need to add the genomic sequence to $genomic_seq_for_MD-tag so we can create a proper MD tag later
			if ($contains_deletion_1){
				$genomic_seq_for_MD_tag_1 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_1,$len_1[$_]);
			}
			#     warn "Just adjusting the position by: ",$len_1[$_],"bp\n";
			$pos_1 += $len_1[$_];
			$indels_1 += $len_1[$_]; # adding to $indels_1 to determine the hemming distance (= single base mismatches, insertions or deletions) for the SAM output
		}
		elsif ($ops_1[$_] eq 'N'){ # skipped region in the read; a splice junction
			# we do not add any genomic sequence but only adjust the position
			#     warn "Just adjusting the position by: ",$len_1[$_],"bp\n";
			$pos_1 += $len_1[$_];
			# not altering the variable needed for the hemming distance
			# $indels_1 += $len_1[$_]; # adding to $indels_1 to determine the hemming distance (= single base mismatches, insertions or deletions) for the SAM output
		}
		elsif($cigar_1 =~ tr/[HPX=]//){ # if these (for standard mapping) illegal characters exist we die
			die "The CIGAR 1 string contained illegal CIGAR operations in addition to 'M', 'I', 'D', 'S' and 'N': $cigar_1";
		}
		else{
			die "The CIGAR 1 string contained undefined CIGAR operations in addition to 'M', 'I', 'D', 'S' and 'N': $cigar_1";
		}	
	}

	### 3' end of read 1
	if ( ($methylation_call_params->{$sequence_identifier}->{index} == 0) or ($methylation_call_params->{$sequence_identifier}->{index} == 2) ){
		## checking if the substring will be valid or if we can't extract the sequence because we are right at the edge of a chromosome
		unless (length($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}}) >= $pos_1+2){# exiting with en empty genomic sequence otherwise
			$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_1} = $non_bisulfite_sequence_1;
			return;
		}

		$non_bisulfite_sequence_1 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_1,2);
	}	


	### Extracting read 2 genomic sequence ###

	### 5' end of read 2
	if ( ($methylation_call_params->{$sequence_identifier}->{index} == 1) or ($methylation_call_params->{$sequence_identifier}->{index} == 3) ){
		## checking if the substring will be valid or if we can't extract the sequence because we are right at the edge of a chromosome
		unless ( ($pos_2-2) >= 0){# exiting with en empty genomic sequence otherwise
			$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_1} = $non_bisulfite_sequence_1;
			$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_2} = $non_bisulfite_sequence_2;
			$methylation_call_params->{$sequence_identifier}->{genomic_seq_for_MD_tag_2} = $genomic_seq_for_MD_tag_2;
			return;
		}
		$non_bisulfite_sequence_2 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_2-2,2);
	}

	foreach (0..$#len_2){
		if ($ops_2[$_] eq 'M'){
			# extracting genomic sequence
			$non_bisulfite_sequence_2 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_2,$len_2[$_]);
			if ($contains_deletion_2){
				$genomic_seq_for_MD_tag_2 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_2,$len_2[$_]);
			}
			# warn "$non_bisulfite_sequence_2\n";
			# adjusting position
			$pos_2 += $len_2[$_];
		}
		elsif ($ops_2[$_] eq 'I'){ # insertion in the read sequence
			# we simply add padding Xs instead of finding genomic sequence. This will not be used to infer methylation calls and we can ignore this later during the generation of the MD:Z: tag
			$non_bisulfite_sequence_2 .= 'X' x $len_2[$_];
			if ($contains_deletion_2){
				$genomic_seq_for_MD_tag_2 .= 'X' x $len_2[$_];
			}
			# warn "$non_bisulfite_sequence_2\n";
			# position doesn't need adjusting

			### 03 06 2014: In fact we don't need to add anything to the hemming distance for insertions since we use padding Xs which will fail a base by base comparison in hemming_dist()
			# $indels_2 += $len_2[$_]; # adding to $indels_1 to determine the hemming distance (= single base mismatches, insertions or deletions) for the SAM output
		}
		elsif ($ops_2[$_] eq 'S'){ # soft-clipped read sequence
			# we simply add padding Xs instead of finding genomic sequence. This will not be used to infer methylation calls
			# and we can ignore this later during the generation of the MD:Z: tag
			$non_bisulfite_sequence_2 .= 'X' x $len_2[$_];
			if ($contains_deletion_2){
				$genomic_seq_for_MD_tag_2 .= 'X' x $len_2[$_];
			}
			# warn "Soft-clipped Read 2sequence! ID: $sequence_identifier\n";
			# warn "$non_bisulfite_sequence_2\n";
			# position doesn't need adjusting

			# We don't need to add anything to the hemming distance for insertions since we use padding Xs which will fail a base by base comparison in hemming_dist()
		}
		elsif ($ops_2[$_] eq 'D'){ # deletion in the read sequence
			# we do not add any genomic sequence but only adjust the position
			# we do however need to add the genomic sequence to $genomic_seq_for_MD-tag so we can create a proper MD tag later
			if ($contains_deletion_2){
				$genomic_seq_for_MD_tag_2 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_2,$len_2[$_]);
			}
			# warn "Just adjusting the position by: ",$len_2[$_],"bp\n";
			$pos_2 += $len_2[$_];
			$indels_2 += $len_2[$_]; # adding to $indels_1 to determine the hemming distance (= single base mismatches, insertions or deletions) for the SAM output
		}
		elsif ($ops_2[$_] eq 'N'){ # deletion in the read sequence
			# we do not add any genomic sequence but only adjust the position
			# warn "Just adjusting the position by: ",$len_2[$_],"bp\n";
			$pos_2 += $len_2[$_];
		}
		elsif($cigar_2 =~ tr/[SHPX=]//){ # if these (for standard mapping) illegal characters exist we die
			die "The CIGAR 2 string contained illegal CIGAR operations in addition to 'M', 'I', 'D', 'S' and 'N': $cigar_2";
		}
		else{
			die "The CIGAR 2 string contained undefined CIGAR operations in addition to 'M', 'I', 'D', 'S' and 'N': $cigar_2";
		}
	}

	### 3' end of read 2
	if ( ($methylation_call_params->{$sequence_identifier}->{index} == 0) or ($methylation_call_params->{$sequence_identifier}->{index} == 2) ){
		## checking if the substring will be valid or if we can't extract the sequence because we are right at the edge of a chromosome
			unless (length($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}}) >= $pos_2+2){# exiting with en empty genomic sequence otherwise
			# need to set read 1 as well now to prevent warning
			#  warn "'$non_bisulfite_sequence_1'\n'$non_bisulfite_sequence_2'\n\n";
			#  sleep(5);
			$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_1} = $non_bisulfite_sequence_1;
			$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_2} = $non_bisulfite_sequence_2;
			return;
		}
		$non_bisulfite_sequence_2 .= substr ($chromosomes{$methylation_call_params->{$sequence_identifier}->{chromosome}},$pos_2,2);
	}

	### all paired-end alignments reported by Bowtie 2 have the Read 1 alignment first and the Read 2 alignment as the second one irrespective of whether read 1 or read 2 was
	### the + alignment. We also read in sequences read 1 then read 2 so they should correspond perfectly

	### results from CT converted read 1 plus GA converted read 2 vs. CT converted genome (+/- orientation alignments are reported only)
	if ($methylation_call_params->{$sequence_identifier}->{index} == 0){
		### [Index 0, sequence originated from (converted) forward strand]
		$counting{CT_GA_CT_count}++;
		$alignment_read_1 = '+';
		$alignment_read_2 = '-';
		$read_conversion_info_1 = 'CT';
		$read_conversion_info_2 = 'GA';
		$genome_conversion = 'CT';
		### Read 1 is always the forward hit
		### Read 2 is will always on the reverse strand, so it needs to be reverse complemented
		$non_bisulfite_sequence_2 = reverse_complement($non_bisulfite_sequence_2);
		if ($contains_deletion_2){
			$genomic_seq_for_MD_tag_2 = reverse_complement($genomic_seq_for_MD_tag_2);
		}
	}

	### results from GA converted read 1 plus CT converted read 2 vs. GA converted genome (+/- orientation alignments are reported only)
	elsif ($methylation_call_params->{$sequence_identifier}->{index} == 1){
		### [Index 1, sequence originated from complementary to (converted) bottom strand]
		$counting{GA_CT_GA_count}++;
		$alignment_read_1 = '+';
		$alignment_read_2 = '-';
		$read_conversion_info_1 = 'GA';
		$read_conversion_info_2 = 'CT';
		$genome_conversion = 'GA';
		### Read 1 is always the forward hit
		### Read 2 is will always on the reverse strand, so it needs to be reverse complemented
		$non_bisulfite_sequence_2 = reverse_complement($non_bisulfite_sequence_2);
		if ($contains_deletion_2){
			$genomic_seq_for_MD_tag_2 = reverse_complement($genomic_seq_for_MD_tag_2);
		}
	}

	### results from GA converted read 1 plus CT converted read 2 vs. CT converted genome (-/+ orientation alignments are reported only)
	elsif ($methylation_call_params->{$sequence_identifier}->{index} == 2){
		### [Index 2, sequence originated from the complementary to (converted) top strand]
		$counting{GA_CT_CT_count}++;
		$alignment_read_1 = '-';
		$alignment_read_2 = '+';
		$read_conversion_info_1 = 'GA';
		$read_conversion_info_2 = 'CT';
		$genome_conversion = 'CT';

		### Read 1 (the reverse strand) genomic sequence needs to be reverse complemented
		$non_bisulfite_sequence_1 = reverse_complement($non_bisulfite_sequence_1);
		if ($contains_deletion_1){
			$genomic_seq_for_MD_tag_1 = reverse_complement($genomic_seq_for_MD_tag_1);
		}
	}

	### results from CT converted read 1 plus GA converted read 2 vs. GA converted genome (-/+ orientation alignments are reported only)
	elsif ($methylation_call_params->{$sequence_identifier}->{index} == 3){
		### [Index 3, sequence originated from the (converted) reverse strand]
		$counting{CT_GA_GA_count}++;
		$alignment_read_1 = '-';
		$alignment_read_2 = '+';
		$read_conversion_info_1 = 'CT';
		$read_conversion_info_2 = 'GA';
		$genome_conversion = 'GA';
		### Read 1 (the reverse strand) genomic sequence needs to be reverse complemented
		$non_bisulfite_sequence_1 = reverse_complement($non_bisulfite_sequence_1);
		if ($contains_deletion_1){
			$genomic_seq_for_MD_tag_1 = reverse_complement($genomic_seq_for_MD_tag_1);
		}
	}
	else{
		die "Too many bowtie result filehandles\n";
	}
	### the alignment_strand information is needed to determine which strand of the genomic sequence we are comparing the read against,
	### the read_conversion information is needed to know whether we are looking for C->T or G->A substitutions

	$methylation_call_params->{$sequence_identifier}->{alignment_read_1} = $alignment_read_1;
	$methylation_call_params->{$sequence_identifier}->{alignment_read_2} = $alignment_read_2;
	$methylation_call_params->{$sequence_identifier}->{genome_conversion} = $genome_conversion;
	$methylation_call_params->{$sequence_identifier}->{read_conversion_1} = $read_conversion_info_1;
	$methylation_call_params->{$sequence_identifier}->{read_conversion_2} = $read_conversion_info_2;
	$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_1} = $non_bisulfite_sequence_1;
	$methylation_call_params->{$sequence_identifier}->{unmodified_genomic_sequence_2} = $non_bisulfite_sequence_2;
	$methylation_call_params->{$sequence_identifier}->{genomic_seq_for_MD_tag_1} = $genomic_seq_for_MD_tag_1;
	$methylation_call_params->{$sequence_identifier}->{genomic_seq_for_MD_tag_2} = $genomic_seq_for_MD_tag_2;

	## the end position of a read is stored in $pos
	$methylation_call_params->{$sequence_identifier}->{end_position_1} = $pos_1;
	$methylation_call_params->{$sequence_identifier}->{end_position_2} = $pos_2;
	$methylation_call_params->{$sequence_identifier}->{indels_1} = $indels_1;
	$methylation_call_params->{$sequence_identifier}->{indels_2} = $indels_2;
}

sub reverse_complement{
	my $sequence = shift;
	$sequence =~ tr/CATG/GTAC/;
	$sequence = reverse($sequence);
	return $sequence;
}

# Compute MAPQ value for a read or read pair as in Bowtie2-2.2.2 (specifically, V2 of the MAPQ calculator: "class BowtieMapq2")
# assuming end-to-end alignment with the default calculation of the minimum alignment score

sub calc_mapq {
	my ($read1Len, $read2Len, $AS_best, $AS_secBest) = @_;

	# Calculate the minimum alignment score either with linear or logarithmic function
	# Bismark hardcodes the expectation that end-to-end alignments will receive a linear score_min function (L,Intercept,Coefficient) while local alignment will receive logarithmic score_min function(G,Intercept,Coefficient)
	# This matches the defaults function forms in bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-score-min
	# For details on scoring functions: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#setting-function-options
	# If this expectation is lifted, the following code will need to account for the user-requested function too.
	# For the moment, we can do:
	my $scMin = $score_min_intercept + $score_min_slope * ($local ? log $read1Len : $read1Len);
	### read2Len is only defined for paired-end reads, so for single-end mode we can just a score min value for read 1
	if (defined $read2Len){
		$scMin += $score_min_intercept + $score_min_slope * ($local ? log $read2Len : $read2Len);
	}

	my $diff = abs$scMin; # scores can vary by up to this much (since max AS is 0 for end-to-end alignment)
	my $bestOver = $AS_best - $scMin;
	#warn "AS_best: $AS_best\n";
	#warn "scMin: $scMin\n";
	#warn "diff: $diff\n";
	#warn "bestOver (AS_best - scMin): $bestOver\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	
	if (!$local){
		# warn "End-to-End alignment\n";
		if (!defined $AS_secBest) {
			if    ($bestOver >= $diff * 0.8) { return 42; }
			elsif ($bestOver >= $diff * 0.7) { return 40; }
			elsif ($bestOver >= $diff * 0.6) { return 24; }
			elsif ($bestOver >= $diff * 0.5) { return 23; }
			elsif ($bestOver >= $diff * 0.4) { return  8; }
			elsif ($bestOver >= $diff * 0.3) { return  3; }
			else                             { return  0; }
		}
		else{
			my $bestDiff = abs(abs($AS_best) - abs($AS_secBest));
			if ($bestDiff >= $diff * 0.9) {
				if ($bestOver == $diff) {
					return 39;
				} 
				else {
					return 33;
				}
			} 
			elsif ($bestDiff >= $diff * 0.8) {
				if ($bestOver == $diff) {
					return 38;
				} 
				else {
					return 27;
				}
			} 
			elsif ($bestDiff >= $diff * 0.7) {
				if ($bestOver == $diff) {
					return 37;
				}
				else {
					return 26;
				}
			} 
			elsif ($bestDiff >= $diff * 0.6) {
				if ($bestOver == $diff) {
					return 36;
				} 
				else {
					return 22;
				}
			}
			elsif ($bestDiff >= $diff * 0.5) {
				if ($bestOver == $diff) {
					return 35;
				} 
				elsif ($bestOver >= $diff * 0.84) {
					return 25;
				}
				elsif ($bestOver >= $diff * 0.68) {
					return 16;
				}
				else {
					return 5;
				}
			} 
			elsif ($bestDiff >= $diff * 0.4) {
				if ($bestOver == $diff) {
					return 34;
				} 
				elsif ($bestOver >= $diff * 0.84) {
					return 21;
				} 
				elsif ($bestOver >= $diff * 0.68) {
					return 14;
				} 
				else {
					return 4;
				}
			} 
			elsif ($bestDiff >= $diff * 0.3) {
				if ($bestOver == $diff) {
					return 32;
				} 
				elsif ($bestOver >= $diff * 0.88) {
					return 18;
				} 
				elsif ($bestOver >= $diff * 0.67) {
					return 15;
				} 
				else {
					return 3;
				}
			} 
			elsif ($bestDiff >= $diff * 0.2) {
				if ($bestOver == $diff) {
					return 31;
				} 
				elsif ($bestOver >= $diff * 0.88) {
					return 17;
				} 
				elsif ($bestOver >= $diff * 0.67) {
					return 11;
				} 
				else {
					return 0;
				}
			} 
			elsif ($bestDiff >= $diff * 0.1) {
				if ($bestOver == $diff) {
					return 30;
				} 
				elsif ($bestOver >= $diff * 0.88) {
					return 12;
				}
				elsif ($bestOver >= $diff * 0.67) {
					return 7;
				} 
				else {
					return 0;
				}
			} 
			elsif ($bestDiff > 0) {
				if ($bestOver >= $diff * 0.67) {
					return 6;
				} 
				else {
					return 2;
				}
			}
			else {
				if ($bestOver >= $diff * 0.67) {
					return 1;
				}
				else {
					return 0;
				}
			}
		}
	}
	else{
		## Local alignment 
		## For more information see here: https://github.com/FelixKrueger/Bismark/issues/260
		if(!defined $AS_secBest) {
				if   ($bestOver >= $diff * 0.8)   { return 44; }
				elsif($bestOver >= $diff * 0.7)   { return 42; }
				elsif($bestOver >= $diff * 0.6)   { return 41; }
				elsif($bestOver >= $diff * 0.5)   { return 36; }
				elsif($bestOver >= $diff * 0.4)   { return 28; }
				elsif($bestOver >= $diff * 0.3)   { return 24; }
				else                              { return 22; }
		}
		else {
			# FK: Not sure what to do about this to be honest secbest = s.paired() ?
			# FK: Not sure what to do about this to be honest s.bestUnchosenCScore().score() : s.bestUnchosenScore(mate1).score();
			my $bestDiff = abs(abs($AS_best) - abs($AS_secBest));
			# warn "bestDiff is: $bestDiff\n";
			if   ($bestDiff >= $diff * 0.9){
				return 40;
			}	
			elsif($bestDiff >= $diff * 0.8){
				return 39;
			}
			elsif($bestDiff >= $diff * 0.7){
				return 38;
			}
			elsif($bestDiff >= $diff * 0.6){
				return 37;
			}
			elsif($bestDiff >= $diff * 0.5){
				if ($bestOver == $diff){
					return 35;
				}
				elsif($bestOver >= $diff * 0.50){
					return 25;
				}
				else{
					return 20;
				}
			} 
			elsif($bestDiff >= $diff * 0.4){
				if ($bestOver == $diff){
					return 34;
				}
				elsif($bestOver >= $diff * 0.50){
					return 21;
				}
				else{
					return 19;
				}
			}
			elsif($bestDiff >= $diff * 0.3){
				if ($bestOver == $diff){
					return 33;
				}
				elsif($bestOver >= $diff * 0.5){
					return 18;
				}
				else{
					return 16;
				}
			}	
			elsif($bestDiff >= $diff * 0.2){
				if ($bestOver == $diff){
					return 32;
				}
				elsif($bestOver >= $diff * 0.5){
					return 17;
				}
				else{
					return 12;
				}
			}
			elsif($bestDiff >= $diff * 0.1){
				if ($bestOver == $diff){ 
					return 31;
				}
				elsif($bestOver >= $diff * 0.5){ 
					return 14; 
				}
				else{ 
					return 9;
				}
			}
			elsif($bestDiff > 0){
				if($bestOver >= $diff * 0.5){
					return 11;
				}
				else{
					return 2;
				}
			} 
			else {
				if($bestOver >= $diff * 0.5){
					return 1;
				}
				else{
					return 0;
				}
			}
		}	
	}
}

sub methylation_call{
	my ($identifier,$sequence_actually_observed,$genomic_sequence,$read_conversion) = @_;
	### splitting both the actually observed sequence and the genomic sequence up into single bases so we can compare them one by one
	my @seq = split(//,$sequence_actually_observed);
	my @genomic = split(//,$genomic_sequence);
	#  print join ("\n",$identifier,$sequence_actually_observed,$genomic_sequence,$read_conversion),"\n";
	### Creating a match-string with different characters for non-cytosine bases (disregarding mismatches here), methyl-Cs or non-methyl Cs in either
	### CpG, CHH or CHG context

  #################################################################
  ### . for bases not involving cytosines                       ###
  ### X for methylated C in CHG context (was protected)         ###
  ### x for not methylated C in CHG context (was converted)     ###
  ### H for methylated C in CHH context (was protected)         ###
  ### h for not methylated C in CHH context (was converted)     ###
  ### Z for methylated C in CpG context (was protected)         ###
  ### z for not methylated C in CpG context (was converted)     ###
  ### U for methylated C in unknown context (was protected)     ###
  ### u for not methylated C in unknwon context (was converted) ###
  #################################################################

  my @match =();
  warn "length of \@seq: ",scalar @seq,"\tlength of \@genomic: ",scalar @genomic,"\n" unless (scalar @seq eq (scalar@genomic-2)); ## CHH changed to -2
  my $methyl_CHH_count = 0;
  my $methyl_CHG_count = 0;
  my $methyl_CpG_count = 0;
  my $methyl_C_unknown_count = 0;
  my $unmethylated_CHH_count = 0;
  my $unmethylated_CHG_count = 0;
  my $unmethylated_CpG_count = 0;
  my $unmethylated_C_unknown_count = 0;

  if ($read_conversion eq 'CT'){
    for my $index (0..$#seq) {
      if ($seq[$index] eq $genomic[$index]) {
	### The residue can only be a C if it was not converted to T, i.e. protected my methylation
	if ($genomic[$index] eq 'C') {
	  ### If the residue is a C we want to know if it was in CpG context or in any other context
	  my $downstream_base = $genomic[$index+1];

	  if ($downstream_base eq 'G'){
	    ++$methyl_CpG_count;
	    push @match,'Z'; # protected C, methylated, in CpG context
	  }
	  elsif ($downstream_base eq 'N' or $downstream_base eq 'X'){ # if the downstream base was an N we cannot really be sure about the sequence context (as it might have been a CG)
	    ++$methyl_C_unknown_count;
	    push @match,'U'; # protected C, methylated, in Unknown context
	  }
	  else {
	    ### C in not in CpG-context, determining the second downstream base context
	    my $second_downstream_base = $genomic[$index+2];

	    if ($second_downstream_base eq 'G'){
	      ++$methyl_CHG_count;
	      push @match,'X'; # protected C, methylated, in CHG context
	    }
	    elsif ($second_downstream_base eq 'N' or $second_downstream_base eq 'X'){
	      ++$methyl_C_unknown_count; # if the second downstream base was an N we cannot really be sure about the sequence context (as it might have been a CHH or CHG)
	      push @match,'U'; # protected C, methylated, in Unknown context
	    }
	    else{
	      ++$methyl_CHH_count;
	      push @match,'H'; # protected C, methylated, in CHH context
	    }
	  }
	}
	else {
	  push @match, '.';
	}
      }
      elsif ($seq[$index] ne $genomic[$index]) {
	### for the methylation call we are only interested in mismatches involving cytosines (in the genomic sequence) which were converted into Ts
	### in the actually observed sequence
	if ($genomic[$index] eq 'C' and $seq[$index] eq 'T') {
	  ### If the residue was converted to T we want to know if it was in CpG, CHG or CHH  context
	  my $downstream_base = $genomic[$index+1];

	  if ($downstream_base eq 'G'){
	    ++$unmethylated_CpG_count;
	    push @match,'z'; # converted C, not methylated, in CpG context
	  }
	  elsif ($downstream_base eq 'N' or $downstream_base eq 'X'){ # if the downstream base was an N we cannot really be sure about the sequence context (as it might have been a CG)
	    ++$unmethylated_C_unknown_count;
	    push @match,'u'; # converted C, not methylated, in Unknown context
	  }
	  else{
	    ### C in not in CpG-context, determining the second downstream base context
	    my $second_downstream_base = $genomic[$index+2];

	    if ($second_downstream_base eq 'G'){
	      ++$unmethylated_CHG_count;
	      push @match,'x'; # converted C, not methylated, in CHG context
	    }
	    elsif ($second_downstream_base eq 'N' or $second_downstream_base eq 'X'){
	      ++$unmethylated_C_unknown_count; # if the second downstream base was an N we cannot really be sure about the sequence context (as it might have been a CHH or CHG)
	      push @match,'u'; # converted C, not methylated, in Unknown context
	    }
	    else{
	      ++$unmethylated_CHH_count;
	      push @match,'h'; # converted C, not methylated, in CHH context
	    }
	  }
	}
	### all other mismatches are not of interest for a methylation call
	else {
	  push @match,'.';
	}
      }
      else{
	die "There can be only 2 possibilities\n";
      }
    }
  }
  elsif ($read_conversion eq 'GA'){
    # print join ("\n",'***',$identifier,$sequence_actually_observed,$genomic_sequence,$read_conversion,'***'),"\n";

    for my $index (0..$#seq) {
      if ($seq[$index] eq $genomic[$index+2]) {
	### The residue can only be a G if the C on the other strand was not converted to T, i.e. protected my methylation
	if ($genomic[$index+2] eq 'G') {
	  ### If the residue is a G we want to know if the C on the other strand was in CpG, CHG or CHH context, therefore we need
	  ### to look if the base upstream is a C

	  my $upstream_base = $genomic[$index+1];

	  if ($upstream_base eq 'C'){
	    ++$methyl_CpG_count;
	    push @match,'Z'; # protected C on opposing strand, methylated, in CpG context
	  }
	  elsif ($upstream_base eq 'N' or $upstream_base eq 'X'){ # if the upstream base was an N we cannot really be sure about the sequence context (as it might have been a CG)
	    ++$methyl_C_unknown_count;
	    push @match,'U'; # protected C on opposing strand, methylated, in Unknown context
	  }
	  else{
	    ### C in not in CpG-context, determining the second upstream base context
	    my $second_upstream_base = $genomic[$index];

	    if ($second_upstream_base eq 'C'){
	      ++$methyl_CHG_count;
	      push @match,'X'; # protected C on opposing strand, methylated, in CHG context
	    }
	    elsif ($second_upstream_base eq 'N' or $second_upstream_base eq 'X'){
	      ++$methyl_C_unknown_count; # if the second upstream base was an N we cannot really be sure about the sequence context (as it might have been a CHH or CHG)
	      push @match,'U'; # protected C, methylated, in Unknown context
	    }
	    else{
	      ++$methyl_CHH_count;
	      push @match,'H'; # protected C on opposing strand, methylated, in CHH context
	    }
	  }
	}
	else{
	  push @match, '.';
	}
      }
      elsif ($seq[$index] ne $genomic[$index+2]) {
	### for the methylation call we are only interested in mismatches involving cytosines (in the genomic sequence) which were converted to Ts
	### on the opposing strand, so G to A conversions in the actually observed sequence
	if ($genomic[$index+2] eq 'G' and $seq[$index] eq 'A') {
	  ### If the C residue on the opposing strand was converted to T then we will see an A in the currently observed sequence. We want to know if
	  ### the C on the opposing strand was it was in CpG, CHG or CHH context, therefore we need to look one (or two) bases upstream!

	  my $upstream_base = $genomic[$index+1];

	  if ($upstream_base eq 'C'){
	    ++$unmethylated_CpG_count;
	    push @match,'z'; # converted C on opposing strand, not methylated, in CpG context
	  }
	  elsif ($upstream_base eq 'N' or $upstream_base eq 'X'){ # if the upstream base was an N we cannot really be sure about the sequence context (as it might have been a CG)
	    ++$unmethylated_C_unknown_count;
	    push @match,'u'; # converted C on opposing strand, not methylated, in Unknown context
	  }
	  else{
	    ### C in not in CpG-context, determining the second upstream base context
	    my $second_upstream_base = $genomic[$index];

	    if ($second_upstream_base eq 'C'){
	      ++$unmethylated_CHG_count;
	      push @match,'x'; # converted C on opposing strand, not methylated, in CHG context
	    }
	    elsif ($second_upstream_base eq 'N' or $second_upstream_base eq 'X'){
	      ++$unmethylated_C_unknown_count; # if the second upstream base was an N we cannot really be sure about the sequence context (as it might have been a CHH or CHG)
	      push @match,'u'; # converted C on opposing strand, not methylated, in Unknown context
	    }
	    else{
	      ++$unmethylated_CHH_count;
	      push @match,'h'; # converted C on opposing strand, not methylated, in CHH context
	    }
	  }
	}
	### all other mismatches are not of interest for a methylation call
	else {
	  push @match,'.';
	}
      }
      else{
	die "There can be only 2 possibilities\n";
      }
    }
  }
  else{
    die "Strand conversion info is required to perform a methylation call\n";
  }

  my $methylation_call = join ("",@match);

  $counting{total_meCHH_count} += $methyl_CHH_count;
  $counting{total_meCHG_count} += $methyl_CHG_count;
  $counting{total_meCpG_count} += $methyl_CpG_count;
  $counting{total_meC_unknown_count} += $methyl_C_unknown_count;
  $counting{total_unmethylated_CHH_count} += $unmethylated_CHH_count;
  $counting{total_unmethylated_CHG_count} += $unmethylated_CHG_count;
  $counting{total_unmethylated_CpG_count} += $unmethylated_CpG_count;
  $counting{total_unmethylated_C_unknown_count} += $unmethylated_C_unknown_count;

  # print "\n$sequence_actually_observed\n$genomic_sequence\n",@match,"\n$read_conversion\n\n";

  return $methylation_call;
}

sub print_bisulfite_mapping_results_paired_ends{
	my ($identifier,$sequence_1,$sequence_2,$methylation_call_params,$quality_value_1,$quality_value_2)= @_;

  ### we will output the FastQ quality in Sanger encoding (Phred 33 scale)
  if ($phred64){
    $quality_value_1 = convert_phred64_quals_to_phred33($quality_value_1);
    $quality_value_2 = convert_phred64_quals_to_phred33($quality_value_2);
  }
 
  ### writing every single aligned read and its methylation call to the output file  (unmapped and ambiguous reads were already printed)
  paired_end_SAM_output($identifier,$sequence_1,$sequence_2,$methylation_call_params,$quality_value_1,$quality_value_2); # at the end of the script

}

sub paired_end_SAM_output{

	my ($id,$actual_seq_1,$actual_seq_2,$methylation_call_params,$qual_1,$qual_2) = @_;
	my $strand_1                = $methylation_call_params->{$id}->{alignment_read_1}; # Bowtie 1 only reports the read 1 alignment strand
	my $strand_2                = $methylation_call_params->{$id}->{alignment_read_2};
	my $chr                     = $methylation_call_params->{$id}->{chromosome};
	my $ref_seq_1               = $methylation_call_params->{$id}->{unmodified_genomic_sequence_1};
	my $ref_seq_2               = $methylation_call_params->{$id}->{unmodified_genomic_sequence_2};
	my $methcall_1              = $methylation_call_params->{$id}->{methylation_call_1};
	my $methcall_2              = $methylation_call_params->{$id}->{methylation_call_2};
	my $read_conversion_1       = $methylation_call_params->{$id}->{read_conversion_1};
	my $read_conversion_2       = $methylation_call_params->{$id}->{read_conversion_2};
	my $genome_conversion       = $methylation_call_params->{$id}->{genome_conversion};

	my $id_1;
	my $id_2;

	if ($old_flag){
		$id_1 = $id.'/1';
		$id_2 = $id.'/2';
	}
	else{
		$id_1 = $id; # appending /1 or /2 confuses some downstream programs such as Picard
		$id_2 = $id;
	}

	# Allows all degenerate nucleotide sequences in reference genome
	# die "Reference sequence ($ref_seq_1) contains invalid nucleotides!\n" if $ref_seq_1 =~ /[^ACTGNRYMKSWBDHVX]/i; # X are padded nucleotides in case of insertions in the read
	# die "Reference sequence ($ref_seq_2) contains invalid nucleotides!\n" if $ref_seq_2 =~ /[^ACTGNRYMKSWBDHVX]/i;

	my $index; # used to store the srand origin of the alignment in a less convoluted way

	if ($read_conversion_1 eq 'CT' and $genome_conversion eq 'CT'){
		$index = 0; ## this is OT   (original top strand)
	}
	elsif ($read_conversion_1 eq 'GA' and $genome_conversion eq 'GA'){
		$index = 1; ## this is CTOB (complementary to OB)
	}
	elsif ($read_conversion_1 eq 'GA' and $genome_conversion eq 'CT'){
		$index = 2; ## this is CTOT (complementary to OT)
	}
	elsif ($read_conversion_1 eq 'CT' and $genome_conversion eq 'GA'){
		$index = 3; ## this is OB   (original bottom)
	}	
	else {
		die "Unexpected combination of read 1 and genome conversion: $read_conversion_1 / $genome_conversion\n";
	}

	my $number_of_mismatches_1  = $methylation_call_params->{$id}->{alignment_score_1}; # only needed for custom allele-specific output, not the default!
	my $number_of_mismatches_2  = $methylation_call_params->{$id}->{alignment_score_2};
	
	### we need to remove 2 bp of the genomic sequence as we were extracting read + 2bp long fragments to make a methylation call at the
	### first or last position.

	if ($index == 0 or $index == 3){ # OT or OB
		$ref_seq_1 = substr($ref_seq_1,0,length($ref_seq_1)-2);
		$ref_seq_2 = substr($ref_seq_2,2,length($ref_seq_2)-2);
	}
	else{ # CTOT or CTOB
		$ref_seq_1 = substr($ref_seq_1,2,length($ref_seq_1)-2);
    	$ref_seq_2 = substr($ref_seq_2,0,length($ref_seq_2)-2);
	}

	
	#####

	# start positions

	my $start_read_1 = $methylation_call_params->{$id}->{position_1};
	my $start_read_2 = $methylation_call_params->{$id}->{position_2};
	
	#####

	# end positions

	my $end_read_1 = $methylation_call_params->{$id}->{end_position_1};
	my $end_read_2 = $methylation_call_params->{$id}->{end_position_2};

	#####

	### This is a description of the bitwise FLAG field which needs to be set for the SAM file taken from: "The SAM Format Specification (v1.4-r985), September 7, 2011"
	## FLAG: bitwise FLAG. Each bit is explained in the following table:
	## Bit    Description                                                Comment                                Value
	## 0x1    template having multiple segments in sequencing            0: single-end 1: paired end            value: 2^^0 (  1)
	## 0x2    each segment properly aligned according to the aligner     true only for paired-end alignments    value: 2^^1 (  2)
	## 0x4    segment unmapped                                           ---                                           ---
	## 0x8    next segment in the template unmapped                      ---                                           ---
	## 0x10   SEQ being reverse complemented                             - strand alignment                     value: 2^^4 ( 16)
	## 0x20   SEQ of the next segment in the template being reversed     + strand alignment                     value: 2^^5 ( 32)
	## 0x40   the first segment in the template                          read 1                                 value: 2^^6 ( 64)
	## 0x80   the last segment in the template                           read 2                                 value: 2^^7 (128)
	## 0x100  secondary alignment                                        ---                                           ---
	## 0x200  not passing quality controls                               ---                                           ---
	## 0x400  PCR or optical duplicate                                   ---                                           ---

	### As the FLAG value do not consider that there might be 4 different bisulfite strands of DNA, we are trying to make FLAG tags which take the strand identity into account

	# strands OT and CTOT will be treated as aligning to the top strand (both sequences are scored as aligning to the top strand)
	# strands OB and CTOB will be treated as aligning to the bottom strand (both sequences are scored as reverse complemented sequences)

	my $flag_1;                                                            # FLAG variable used for SAM format
	my $flag_2;

	### The new default FLAG values were changed on 21 07 2015, so that reads do not ignored as discordant reads by the new SeqMonk BAM import
	### In essence we are going to flip the R1 R2 flags around for CTOT and CTOB reads. We still report the first and second read in the same
	### order and only change the actual FLAG value. This should not affect the methylation extraction in any way

	if ($index == 0){       # OT
		unless ($old_flag){
			$flag_1 = 99;                                                      # Read 1 is on the + strand and Read 2 is reversed  (1+2+32+64)
			$flag_2 = 147;                                                     # Read 2 is reverse complemented but informative for the OT  (1+2+16+128)
		}
		else{
			$flag_1 = 67;                                                      # Read 1 is on the + strand  (1+2+64) (Read 2 is technically reverse-complemented, but we do not score it)
			$flag_2 = 131;                                                     # Read 2 is on - strand but informative for the OT        (1+2+128)
		}
	}
	elsif ($index == 1){    # CTOB
		unless($old_flag){
			$flag_1 = 163;                                               # Read 1 is on the forward strand (CTOB) and Read 2 is reverse complemented but we swap round the FLAG
																		 # for R1 and R2 so that we don't end up with discordant pairs
                                                                         # So Read 1 gets Paired read, mapped in proper pair, mate is reversed and second in pair  (1+2+32+128)
			$flag_2 = 83;                                                      # Read 2 gets Read paired, mapped in proper pair, first in pair and Read 2 is reversed  (1+2+16+64)
		}
		else{
			$flag_1 = 115;                                                     # Read 1 is on the + strand, we score for OB  (1+2+16+32+64)
			$flag_2 = 179;                                                     # Read 2 is on the - strand  (1+2+16+32+128)
		}
	}
	elsif ($index == 2){    # CTOT
		unless ($old_flag){
			$flag_1 = 147;                                                     # Read 1 is reverse complemented (CTOT) and Read 2 is the forward read
                                                                         # but we swap round the FLAG for R1 and R2 so that we do not end up with discordant pairs
                                                                         # So Read 1 gets Read paired, read mapped in proper pair, read reverse complemented and second in pair (1+2+32+128)
			$flag_2 = 99;                                                      # Read 2 gets Read paired, read mapped in proper pair, mate reverse strand and First in Pair (1+2+32+64)
		}
		else{
			$flag_1 = 67;                                                      # Read 1 is on the - strand (CTOT) strand, but we score it for OT (1+2+64)
			$flag_2 = 131;                                                     # Read 2 is on the + strand, score it for OT (1+2+128)
		}
	}
	elsif ($index == 3){    # OB
		unless ($old_flag){
			$flag_1 = 83;                                                      # Read 1 is on the - strand, mapped in proper pair and Read 1 is reversed  (1+2+16+64)
			$flag_2 = 163;                                                     # Read 2 is on the - strand, mapped in proper pair and Read 1 is reversed  (1+2+32+128)
		}
		else{
			$flag_1 = 115;                                                     # Read 1 is on the - strand, we score for OB  (1+2+16+32+64)
			$flag_2 = 179;                                                     # Read 2 is on the + strand  (1+2+16+32+128)
		}
	}

	#####

	my $mapq = $methylation_call_params->{$id}->{mapq};
	
	#####
	
	my $cigar_1 = $methylation_call_params->{$id}->{CIGAR_1};             # Actual CIGAR string reported by Bowtie 2
	my $cigar_2 = $methylation_call_params->{$id}->{CIGAR_2};
	
	#####

	my $rnext = '=';                                                     # Chromosome of mate; applies to both reads

	#####

	my $pnext_1 = $start_read_2;                                         # Leftmost position of mate
	my $pnext_2 = $start_read_1;

	#####

	my $tlen_1;                                                          # signed observed Template LENgth (or inferred fragment size)
	my $tlen_2;

	if ($start_read_1 <= $start_read_2){

		# Read 1 alignment is leftmost

		if ($end_read_2 >= $end_read_1){

			if ($flag_1 == 83 and $dovetail){   # R1 has a reverse orientation
				#         ----------------->     read 2   reads are dovetailing, that is one mate alignment extends past the beginning of the other
				#  <-------------------          read 1   such that the wrong mate begins upstream
				# warn "FLAG 1: $flag_1\nFLAG 2: $flag_2\n";
				# warn "Reads are dovetailing\n";
				$tlen_1 = $start_read_1 - $end_read_2 - 1;     # Read 1 still receives a - sign even though it is the leftmost one
				$tlen_2 = $end_read_2 - $start_read_1 + 1;     # Read 2 receives a + sign,
				# warn "TLEN 1: $tlen_1\nTLEN 2: $tlen_2\n";
			}
			else{
				# ------------->                 read 1   reads not overlapping
				#                 <----------    read 2
				#             or
				# 	------------------->           read 1   reads overlapping
				#        <-------------------    read 2
				#             or
				# ------------------------->     read 1
				#   <-----------------------     read 2   read 2 contained within read 1
				#             or
				# ------------------------->     read 1   reads 1 and 2 exactly overlapping
				# <-------------------------     read 2
				#

				$tlen_1 = $end_read_2 - $start_read_1 + 1;                         # Leftmost read has a + sign,
				$tlen_2 = $start_read_1 - $end_read_2 - 1;                         # Rightmost read has a - sign
				# warn "Reads are non/overlapping\nTLEN 1: $tlen_1\nTLEN 2: $tlen_2\n";
			}
		}
		elsif ($end_read_2 < $end_read_1){

			# ------------------------->     read 1
			#       <-----------             read 2   read 2 contained within read 1
			#
			# or
			#
			# ------------------------->     read 1
			# <------------------------      read 2   read 2 contained within read 1

			# start and end of read 2  are fully contained within read 1, using the length of read 1 for the TLEN variable
			$tlen_1 = $end_read_1 - $start_read_1 + 1;          # Set to length of read 1   Leftmost read has a + sign,
			$tlen_2 = ($end_read_1 - $start_read_1 + 1) * -1;   # Set to length of read 1   Rightmost read has a - sign. well this is debatable. Changed this
			### as a request by frozenlyse on SeqAnswers on 24 July 2013
		}
	}
	elsif ($start_read_2 < $start_read_1){

		# Read 2 alignment is leftmost

		if ($end_read_1 >= $end_read_2){

			# Read 2 alignment is leftmost
			if ($flag_1 == 99 and $dovetail){   # R1 has a forward orientation

				#         ----------------->     read 1   reads are dovetailing, that is one mate alignment extends past the beginning of the other
				#  <-------------------          read 2   such that the wrong mate begins upstream

				# warn "FLAG 1: $flag_1\nFLAG 2: $flag_2\n";
				# warn "Reads are dovetailing\n";
				$tlen_1 = $end_read_1 - $start_read_2 + 1;     # Read 1 still receives a + sign even though it is not leftmost
				$tlen_2 = $start_read_2 - $end_read_1 - 1;
				# warn "TLEN 1: $tlen_1\nTLEN 2: $tlen_2\n";
			}
			else{
				# ------------->                 read 2   reads not overlapping
				#                 <----------    read 1
				#             or
				# ------------------------->     read 2   reads overlapping
				#  <-------------------------    read 1
				#             or
				# ------------------------->     read 2
				#   <-----------------------     read 1   read 1 contained within read 2
				#             or
				# ------------------------->     read 2
				#   <-----------------------     read 1   read 1 contained within read 2
				# warn "FLAG 1: $flag_1\nFLAG 2: $flag_2\n";
				# warn "Read 2 has a forward orientation\n";
				$tlen_2 = $end_read_1 - $start_read_2 + 1;                         # Leftmost read has a + sign,
				$tlen_1 = $start_read_2 - $end_read_1 - 1;                         # Rightmost read has a - sign
			}
		}
		elsif ($end_read_1 < $end_read_2){

				# ------------------------->     read 2
				#       <-----------             read 1   read 1 contained within read 2
				#
				# or
				#
				# ------------------------->     read 2
				#  <------------------------      read 1   read 1 contained within read 2

				# start and end of read 1  are fully contained within read 2, using the length of read 2 for the TLEN variable
				$tlen_1 = ($end_read_2 - $start_read_2 + 1) * -1;          # Set to length of read 2   Shorter read receives a - sign,
				$tlen_2 = $end_read_2 - $start_read_2 + 1;                 # Set to length of read 2   Longer read receives a +. Well this is debatable. Changed this
				### as a request by frozenlyse on SeqAnswers on 24 July 2013
		}
	}


	#####

	# adjusting the strand of the sequence before we use them to generate mismatch strings
	if ($strand_1 eq '-'){
		$actual_seq_1 = revcomp($actual_seq_1);                            # Sequence represented on the forward genomic strand
		$ref_seq_1 = revcomp($ref_seq_1);                                  # Required for comparison with actual sequence
		if ($cigar_1 =~ /[D]/){ # deletion or spliced read
			$methylation_call_params->{$id}->{genomic_seq_for_MD_tag_1} = revcomp( $methylation_call_params->{$id}->{genomic_seq_for_MD_tag_1} );
		}
		$qual_1 = reverse $qual_1;                                         # we need to reverse the quality string as well
	}
	if ($strand_2 eq '-'){
		$actual_seq_2 = revcomp($actual_seq_2);                            # Mate sequence represented on the forward genomic strand
		$ref_seq_2 = revcomp($ref_seq_2);                                  # Required for comparison with actual sequence
		if ($cigar_2 =~ /[D]/){ # deletion or spliced read
			$methylation_call_params->{$id}->{genomic_seq_for_MD_tag_2} = revcomp( $methylation_call_params->{$id}->{genomic_seq_for_MD_tag_2} );
		}
		$qual_2 = reverse $qual_2;                                         # If the sequence gets reverse complemented we reverse the quality string as well
	}

	# print "$actual_seq_1\n$ref_seq_1\n\n";
	# print "$actual_seq_2\n$ref_seq_2\n\n";

	#####

	my $hemming_dist_1 = hemming_dist($actual_seq_1,$ref_seq_1);         # Minimal number of one-nucleotide edits needed to transform the read string into the reference sequence
	my $hemming_dist_2 = hemming_dist($actual_seq_2,$ref_seq_2);
	$hemming_dist_1 += $methylation_call_params->{$id}->{indels_1};    # Adding the number of inserted/deleted bases which we parsed while getting the genomic sequence
	$hemming_dist_2 += $methylation_call_params->{$id}->{indels_2};    # Adding the number of inserted/deleted bases which we parsed while getting the genomic sequence
	
	my $NM_tag_1 = "NM:i:$hemming_dist_1";                               # Optional tag NM: edit distance based on nucleotide differences
	my $NM_tag_2 = "NM:i:$hemming_dist_2";                               # Optional tag NM: edit distance based on nucleotide differences

	#####

	my $MD_tag_1 = make_mismatch_string($actual_seq_1,$ref_seq_1,$cigar_1,$methylation_call_params->{$id}->{genomic_seq_for_MD_tag_1}); # Optional tag MD: String providing mismatched reference bases in the alignment (including indel information)
	my $MD_tag_2 = make_mismatch_string($actual_seq_2,$ref_seq_2,$cigar_2,$methylation_call_params->{$id}->{genomic_seq_for_MD_tag_2});

	#  my $XX_tag_1 = make_mismatch_string($actual_seq_1,$ref_seq_1);       # Optional tag XX: String providing mismatched reference bases in the alignment (NO indel information!)
	#  my $XX_tag_2 = make_mismatch_string($actual_seq_2,$ref_seq_2);

	#####

	my $XM_tag_1;                                                        # Optional tag XM: Methylation call string
	my $XM_tag_2;

	if ($strand_1 eq '-'){
		$XM_tag_1 = 'XM:Z:'.reverse $methcall_1;                           # Needs to be reversed if the sequence was reverse complemented
	}
	else{
		$XM_tag_1 = "XM:Z:$methcall_1";
	}

	if ($strand_2 eq '-'){
		$XM_tag_2 = 'XM:Z:'.reverse $methcall_2;                           # Needs to be reversed if the sequence was reverse complemented
	}
	else{
		$XM_tag_2 = "XM:Z:$methcall_2";
	}

	#####

	my $XR_tag_1 = "XR:Z:$read_conversion_1";                            # Optional tag XR: Read 1 conversion state
	my $XR_tag_2 = "XR:Z:$read_conversion_2";                            # Optional tag XR: Read 2 conversion state

	#####

	my $XG_tag = "XG:Z:$genome_conversion";                              # Optional tag XG: Genome Conversion state; valid for both reads
	
	#####

	# Optionally calculating number of mismatches for Bowtie 2 alignments

	if ($non_bs_mm) {
		$number_of_mismatches_1 =~ s/-//; # removing the minus sign
		$number_of_mismatches_2 =~ s/-//;

		### We need to analyse the CIGAR strings whether the reads contained any indels to determine the number of mismatches

		### CIGAR 1
		if ($cigar_1 =~ /(D|I)/) {
			# warn "$cigar_1\n";

			# parsing CIGAR string
			my @len = split (/\D+/,$cigar_1); # storing the length per operation
			my @ops = split (/\d+/,$cigar_1); # storing the operation
			shift @ops;		# remove the empty first element
			die "CIGAR string '$cigar_1' contained a non-matching number of lengths and operations\n" unless (scalar @len == scalar @ops);

			foreach (0..$#len) {
				if ($ops[$_] eq 'M') {
					# warn "skipping\n";
					next;		# irrelevant
				}
				elsif ($ops[$_] eq 'I') {	# insertion in the read sequence
					$number_of_mismatches_1 -= $insertion_open;
					$number_of_mismatches_1 -= $len[$_] * $insertion_extend;
					# warn "Insertion: Subtracting $ops[$_], length $len[$_], open: $insertion_open, extend: $insertion_extend\n";
				}
				elsif ($ops[$_] eq 'D') {	# deletion in the read sequence
					$number_of_mismatches_1 -= $deletion_open;
					$number_of_mismatches_1 -= $len[$_] * $deletion_extend;
					# warn "Deletion: Subtracting $ops[$_], length $len[$_], open: $deletion_open, extend: $deletion_extend\n";
				}
				elsif ($ops[$_] eq 'N') { # skipped portion, spliced read
					# warn "skipping\n";
					next;		# irrelevant
				}
				elsif ($cigar_1 =~ tr/[HPSX=]//) {	# if these (for standard mapping) illegal characters exist we die
					die "The CIGAR string contained illegal CIGAR operations in addition to 'M', 'I', 'D' and 'N': $cigar_1";
				}
				else {
					die "The CIGAR string contained undefined CIGAR operations in addition to 'M', 'I', 'D' and 'N': $cigar_1";
				}
			}

			# warn "Alignment score $number_of_mismatches_1\n";
			# print "Mismatches $number_of_mismatches_1\n\n";
		}

		### CIGAR 2
		if ($cigar_2 =~ /(D|I)/) {
			# warn "$cigar_2\n";

			# parsing CIGAR string
			my @len = split (/\D+/,$cigar_2); # storing the length per operation
			my @ops = split (/\d+/,$cigar_2); # storing the operation
			shift @ops;		# remove the empty first element
			die "CIGAR string '$cigar_2' contained a non-matching number of lengths and operations\n" unless (scalar @len == scalar @ops);

			foreach (0..$#len) {
				if ($ops[$_] eq 'M') {
					# warn "skipping\n";
					next; #irrelevant
				}
				elsif ($ops[$_] eq 'I') {	# insertion in the read sequence
					$number_of_mismatches_2 -= $insertion_open;
					$number_of_mismatches_2 -= $len[$_] * $insertion_extend;
					# warn "Insertion: Subtracting $ops[$_], length $len[$_], open: $insertion_open, extend: $insertion_extend\n";
				}
				elsif ($ops[$_] eq 'D') {	# deletion in the read sequence
					$number_of_mismatches_2 -= $deletion_open;
					$number_of_mismatches_2 -= $len[$_] * $deletion_extend;
					# warn "Deletion: Subtracting $ops[$_], length $len[$_], open: $deletion_open, extend: $deletion_extend\n";
				}
				elsif ($ops[$_] eq 'N') { # skipped portion, spliced-read
					# warn "skipping\n";
					next; #irrelevant
				}
				elsif ($cigar_2 =~ tr/[SHPX=]//) {	# if these (for standard mapping) illegal characters exist we die
					die "The CIGAR string contained illegal CIGAR operations in addition to 'M', 'I','D' and 'N': $cigar_2";
				}
				else {
					die "The CIGAR string contained undefined CIGAR operations in addition to 'M', 'I', 'D' and 'N': $cigar_2";
				}
			}
		}


		### Now we have InDel corrected Alignment scores

		### if the actual sequence contained Ns we need to adjust the number of mismatches. Ns receive a penalty of -1, 
		### but normal mismatches receive -6. This might still break if the sequence contained more than 5 Ns, but this should occur close to never

		my $seq_1_N_count = $number_of_mismatches_1 % 6; # modulo 6 will return the integer rest after the division
		my $seq_2_N_count = $number_of_mismatches_2 % 6;
		#   warn "N count 1: $seq_1_N_count\n";
		#   warn "N count 2: $seq_2_N_count\n";

		$number_of_mismatches_1 = int ($number_of_mismatches_1 / 6) + $seq_1_N_count;
		$number_of_mismatches_2 = int ($number_of_mismatches_2 / 6) + $seq_2_N_count;

		# warn "MM1    $number_of_mismatches_1 \n";
		# warn "MM2    $number_of_mismatches_2 \n";
		
	}

	####

	my $XA_tag = "XA:Z:$number_of_mismatches_1";
	my $XB_tag = "XB:Z:$number_of_mismatches_2";

	####

	my $read_group; # optional
	if ($rg_tag){
		$read_group = "RG:Z:$rg_id";
	}

	####

	# SAM format: QNAME, FLAG, RNAME, 1-based POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, optional fields
	### optionally print number of non-bisulfite mismatches
	if ($non_bs_mm){
		if ($rg_tag){
			print OUT join("\t", ($id_1, $flag_1, $chr, $start_read_1, $mapq, $cigar_1, $rnext, $pnext_1, $tlen_1, $actual_seq_1, $qual_1, $NM_tag_1, $MD_tag_1, $XM_tag_1,$XR_tag_1,$XG_tag,$XA_tag,$read_group)), "\n";
			print OUT join("\t", ($id_2, $flag_2, $chr, $start_read_2, $mapq, $cigar_2, $rnext, $pnext_2, $tlen_2, $actual_seq_2, $qual_2, $NM_tag_2, $MD_tag_2, $XM_tag_2,$XR_tag_2,$XG_tag,$XB_tag,$read_group)), "\n";
		}
		else{
			print OUT join("\t", ($id_1, $flag_1, $chr, $start_read_1, $mapq, $cigar_1, $rnext, $pnext_1, $tlen_1, $actual_seq_1, $qual_1, $NM_tag_1, $MD_tag_1, $XM_tag_1,$XR_tag_1,$XG_tag,$XA_tag)), "\n";
			print OUT join("\t", ($id_2, $flag_2, $chr, $start_read_2, $mapq, $cigar_2, $rnext, $pnext_2, $tlen_2, $actual_seq_2, $qual_2, $NM_tag_2, $MD_tag_2, $XM_tag_2,$XR_tag_2,$XG_tag,$XB_tag)), "\n";
		}
	}
	else{ # default
		if ($rg_tag){
			print OUT join("\t", ($id_1, $flag_1, $chr, $start_read_1, $mapq, $cigar_1, $rnext, $pnext_1, $tlen_1, $actual_seq_1, $qual_1, $NM_tag_1, $MD_tag_1, $XM_tag_1,$XR_tag_1,$XG_tag,$read_group)), "\n";
			print OUT join("\t", ($id_2, $flag_2, $chr, $start_read_2, $mapq, $cigar_2, $rnext, $pnext_2, $tlen_2, $actual_seq_2, $qual_2, $NM_tag_2, $MD_tag_2, $XM_tag_2,$XR_tag_2,$XG_tag,$read_group)), "\n";
		}
		else{
			print OUT join("\t", ($id_1, $flag_1, $chr, $start_read_1, $mapq, $cigar_1, $rnext, $pnext_1, $tlen_1, $actual_seq_1, $qual_1, $NM_tag_1, $MD_tag_1, $XM_tag_1,$XR_tag_1,$XG_tag)), "\n";
			print OUT join("\t", ($id_2, $flag_2, $chr, $start_read_2, $mapq, $cigar_2, $rnext, $pnext_2, $tlen_2, $actual_seq_2, $qual_2, $NM_tag_2, $MD_tag_2, $XM_tag_2,$XR_tag_2,$XG_tag)), "\n";
		
			#print join("\t", ("READ 1:",$id_1, $flag_1, $chr, $start_read_1, $mapq, $cigar_1, $rnext, $pnext_1, $tlen_1, $actual_seq_1, $qual_1, $NM_tag_1, $MD_tag_1, $XM_tag_1,$XR_tag_1,$XG_tag)), "\n";
			#print join("\t", ("READ 2:",$id_2, $flag_2, $chr, $start_read_2, $mapq, $cigar_2, $rnext, $pnext_2, $tlen_2, $actual_seq_2, $qual_2, $NM_tag_2, $MD_tag_2, $XM_tag_2,$XR_tag_2,$XG_tag)), "\n";
			#print "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
		}
	}
}

sub revcomp{
  my $seq = shift or die "Missing seq to reverse complement\n";
  $seq = reverse $seq;
  $seq =~ tr/ACTGactg/TGACTGAC/;
  return $seq;
}

sub hemming_dist{
	my $matches = 0;
	my @actual_seq = split //,(shift @_);
	my @ref_seq = split //,(shift @_);

	foreach (0..$#actual_seq){
		++$matches if ($actual_seq[$_] eq $ref_seq[$_]);
	}
	return my $hd = scalar @actual_seq - $matches;
}

sub merge_individual_BAM_files{

  my ($tempbam,$original_filename,$single_end) = @_;
  my $merged_name = $original_filename;

  #warn "merged name is: $merged_name\n";
  $merged_name =~ s/.*\///; # deleting path information
  # warn "merged name is: $merged_name\n";
  $merged_name =~ s/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//; # attempting to remove fastq.gz etc to make filename a little shorter
  # warn "merged name is: $merged_name\n"; sleep(5);

  foreach my $temp_bam(@$tempbam){
      $temp_bam =~ s/.*\///; # deleting path information
  }

  if ($prefix){
    $merged_name = "$prefix.$merged_name";
  }

  if ($single_end){
    if ($bowtie2){ # BAM format is the default for Bowtie 2
      $merged_name .= '_bismark_bt2.bam';
    }
    else{          # BAM is the default output
      $merged_name .= '_bismark_hisat2.bam';
    }

    if ($basename){ # Output file basename is set using the -B argument
      $merged_name = "${basename}.bam";
    }
  }
  else{
    if ($bowtie2){ # BAM format is the default for Bowtie 2
      $merged_name .= '_bismark_bt2_pe.bam';
    }
    else{          # BAM is the default output
      $merged_name .= '_bismark_hisat2_pe.bam';
    }

    if ($basename){ # Output file basename is set using the -B argument
      $merged_name = "${basename}_pe.bam";
    }
  }


  if ($cram){
      $merged_name =~ s/bam$/cram/;
      warn "At this stage we write out a single CRAM file and delete all temporary BAM files\n";
      warn "Now merging BAM files @$tempbam into >>> $merged_name <<<\n";
      $final_output_filename = "${output_dir}${merged_name}";

      open (OUT,"| $samtools_path view -h -C -T $cram_ref 2>/dev/null - > ${output_dir}${merged_name}") or die "Failed to write to CRAM file $merged_name: $!\nPlease note that this option requires Samtools version 1.2 or higher!\n\n";
  }
  else{
      $final_output_filename = "${output_dir}${merged_name}";
      warn "Now merging BAM files @$tempbam into >>> $merged_name <<<\n";
      open (OUT,"| $samtools_path view -bSh 2>/dev/null - > ${output_dir}${merged_name}") or die "Failed to write to $merged_name: $!\n";
  }

  my $first = 0;

  foreach my $temp_bam(@$tempbam){
    # $temp_bam =~ s/.*\///; # deleting path information

    warn "Merging from file >> $temp_bam <<\n";

    if ($first > 0){
      open (IN,"$samtools_path view ${output_dir}${temp_bam} |") or die "Failed to read from BAM file ${output_dir}${temp_bam}\n";
    }
    else{ # only for the first file we print the header as well
      open (IN,"$samtools_path view -h ${output_dir}${temp_bam} |") or die "Failed to read from BAM file ${output_dir}${temp_bam}\n";
    }

    while (<IN>){
      print OUT;
    }
    close IN or warn "Failed to close filehandle\n";
    ++$first;
  }
  warn "\n";

  close OUT or warn "Failed to close output filehandle\n\n";

}



sub merge_individual_mapping_reports{

  my ($temp_reports,$original_filename_1,$single_end,$original_filename_2) = @_;
  my $report_file = $original_filename_1;
  $report_file =~ s/.*\///; # removing path information
  $report_file =~ s/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//; # attempting to remove fastq.gz etc to make filename a little shorter

  if ($prefix){
    $report_file = "${prefix}.${report_file}";
  }

  if ($basename){ # Output file basename is set using the -B argument
    $report_file = ${basename};
  }

	if ($single_end){
		if ($bowtie2){
			$report_file .= '_bismark_bt2_SE_report.txt';
		}
		else{
			$report_file .= '_bismark_hisat2_SE_report.txt';
		}
	}
	else{
		if ($bowtie2){
			$report_file .= '_bismark_bt2_PE_report.txt';
		}
		else{
			$report_file .= '_bismark_hisat2_PE_report.txt';
		}
	}
	warn "Writing report to ${output_dir}${report_file}\n";
	open (REPORT,'>',"$output_dir$report_file") or die "Failed to write to ${output_dir}${report_file}: $!\n";

	foreach my $temp(@$temp_reports){
		$temp =~ s/.*\///; # removing path information
	}

	warn "Now merging temporary reports @$temp_reports into >>> ${output_dir}${report_file} <<<\n";

	if ($single_end){
		print REPORT "Bismark report for: $original_filename_1 (version: $bismark_version)\n";
	}
	else{ # paired-end
		print REPORT "Bismark report for: $original_filename_1 and $original_filename_2 (version: $bismark_version)\n";
	}


  my $first = 0;

  foreach my $temp(@$temp_reports){
    # $temp =~ s/.*\///; # removing path information

    warn "Merging from file >> $temp <<\n";
    open (IN,"${output_dir}${temp}") or die "Failed to read from temporary mapping report '${output_dir}${temp}'\n";

    ### this is printing the first couple of lines
    while (<IN>){
      chomp;
      if ($_ =~ /^Bismark report/){
	next;
      }

      unless ($first){ # only happens for the first run we are processing
	if ($_ =~ /^Final Alignment/){
	  ++$first;
	  last;
	}
	else{
	  print REPORT "$_\n";
	}
      }
    }
    close IN or warn "Failed to close filehandle\n";

    ### Simon says: You are going to regret this in the future. Just for the record. He might be right...
    read_alignment_report($temp,$single_end);

  }
  warn "\n";

}

sub read_alignment_report{
  my ($report,$single_end) = @_;

  my $unique;
  my $no_aln;
  my $multiple;
  my $no_genomic;
  my $total_seqs;
  my $bismark_version;
  my $input_filename;

  my $unique_text;
  my $no_aln_text;
  my $multiple_text;
  my $total_seq_text;

  my $total_C_count;
  my ($meth_CpG,$meth_CHG,$meth_CHH,$meth_unknown);
  my ($unmeth_CpG,$unmeth_CHG,$unmeth_CHH,$unmeth_unknown);

  my $number_OT;
  my $number_CTOT;
  my $number_CTOB;
  my $number_OB;

  open (ALN,"${output_dir}${report}") or die "Failed to read from temporary mapping report '$output_dir$report'\n";

  while (<ALN>){
    chomp;

    ### General Alignment stats
    if ($_ =~ /^Sequence pairs analysed in total:/ ){ ## Paired-end
      (undef,$total_seqs) = split /\t/;
      # warn "Total paired seqs: >> $total_seqs <<\n";
    }
    elsif ($_ =~ /^Sequences analysed in total:/ ){   ## Single-end
      (undef,$total_seqs) = split /\t/;
      # warn "total single-end seqs >> $total_seqs <<\n";
    }

    elsif($_ =~ /^Number of paired-end alignments with a unique best hit:/){ ## Paired-end
      (undef,$unique) = split /\t/;
      # warn "Unique PE>> $unique <<\n";
    }
    elsif($_ =~ /^Number of alignments with a unique best hit from/){        ## Single-end
      (undef,$unique) = split /\t/;
      # warn "Unique SE>> $unique <<\n";
    }

    elsif($_ =~ /^Sequence pairs with no alignments under any condition:/){  ## Paired-end
      (undef,$no_aln) = split /\t/;
      # warn "No alignment PE >> $no_aln <<\n";
    }
    elsif($_ =~ /^Sequences with no alignments under any condition:/){  ## Single-end
      (undef,$no_aln) = split /\t/;
      # warn "No alignments SE>> $no_aln <<\n";
    }

    elsif($_ =~ /^Sequence pairs did not map uniquely:/){ ## Paired-end
      (undef,$multiple) = split /\t/;
      # warn "Multiple alignments PE >> $multiple <<\n";
    }
    elsif($_ =~ /^Sequences did not map uniquely:/){ ## Single-end
      (undef,$multiple) = split /\t/;
      # warn "Multiple alignments SE >> $multiple <<\n";
    }

    elsif($_ =~ /^Sequence pairs which were discarded because genomic sequence could not be extracted:/){ ## Paired-end
      (undef,$no_genomic) = split /\t/;
      # warn "No genomic sequence PE >> $no_genomic <<\n";
    }
    elsif($_ =~ /^Sequences which were discarded because genomic sequence could not be extracted:/){ ## Single-end
      (undef,$no_genomic) = split /\t/;
      # warn "No genomic sequence SE>> $no_genomic <<\n";
    }

    ### Context Methylation
    elsif($_ =~ /^Total number of C/ ){
      (undef,$total_C_count) = split /\t/;
      # warn "Total number C >> $total_C_count <<\n";
    }

    elsif($_ =~ /^Total methylated C\'s in CpG context:/ ){
      (undef,$meth_CpG) = split /\t/;
      # warn "meth CpG >> $meth_CpG <<\n" ;
    }
    elsif($_ =~ /^Total methylated C\'s in CHG context:/ ){
      (undef,$meth_CHG) = split /\t/;
      # warn "meth CHG >> $meth_CHG <<\n" ;
    }
    elsif($_ =~ /^Total methylated C\'s in CHH context:/ ){
      (undef,$meth_CHH) = split /\t/;
      # warn "meth CHH >> $meth_CHH <<\n" ;
    }
    elsif($_ =~ /^Total methylated C\'s in Unknown context:/ ){
      (undef,$meth_unknown) = split /\t/;
      # warn "meth Unknown >> $meth_unknown <<\n" ;
    }

    elsif($_ =~ /^Total unmethylated C\'s in CpG context:/ or $_ =~ /^Total C to T conversions in CpG context:/){
      (undef,$unmeth_CpG) = split /\t/;
      # warn "unmeth CpG >> $unmeth_CpG <<\n" ;
    }
    elsif($_ =~ /^Total unmethylated C\'s in CHG context:/ or $_ =~ /^Total C to T conversions in CHG context:/){
      (undef,$unmeth_CHG) = split /\t/;
      # warn "unmeth CHG >> $unmeth_CHG <<\n" ;
    }
    elsif($_ =~ /^Total unmethylated C\'s in CHH context:/ or $_ =~ /^Total C to T conversions in CHH context:/){
      (undef,$unmeth_CHH) = split /\t/;
      # warn "unmeth CHH >> $unmeth_CHH <<\n";
    }
    elsif($_ =~ /^Total unmethylated C\'s in Unknown context:/ or $_ =~ /^Total C to T conversions in Unknown context:/){
      (undef,$unmeth_unknown) = split /\t/;
      # warn "unmeth Unknown >> $unmeth_unknown <<\n" ;
    }

    ### Strand Origin

    elsif($_ =~ /^CT\/GA\/CT:/ ){             ## Paired-end
      (undef,$number_OT) = split /\t/;
      # warn "Number OT PE>> $number_OT <<\n" ;
    }
    elsif($_ =~ /^CT\/CT:/ ){                 ## Single-end
      (undef,$number_OT) = split /\t/;
      # warn "Number OT SE>> $number_OT <<\n" ;
    }

    elsif($_ =~ /^GA\/CT\/CT:/ ){             ## Paired-end
      (undef,$number_CTOT) = split /\t/;
      # warn "Number CTOT PE >> $number_CTOT <<\n" ;
    }
    elsif($_ =~ /^GA\/CT:/ ){                 ## Single-end
      (undef,$number_CTOT) = split /\t/;
      # warn "Number CTOT SE >> $number_CTOT <<\n" ;
    }

    elsif($_ =~ /^GA\/CT\/GA:/ ){             ## Paired-end
      (undef,$number_CTOB) = split /\t/;
      # warn "Number CTOB PE >> $number_CTOB <<\n" ;
    }
    elsif($_ =~ /^GA\/GA:/ ){                 ## Single-end
      (undef,$number_CTOB) = split /\t/;
      # warn "Number CTOB SE >> $number_CTOB <<\n";
    }

    elsif($_ =~ /^CT\/GA\/GA:/ ){             ## Paired-end
      (undef,$number_OB) = split /\t/;
      # warn "Number OB PE >> $number_OB <<\n";
    }
    elsif($_ =~ /^CT\/GA:/ ){                 ## Single-end
      (undef,$number_OB) = split /\t/;
      # warn "Number OB SE >> $number_OB <<\n";
    }
  }

  $counting{sequences_count}                               += $total_seqs;
  $counting{unique_best_alignment_count}                   += $unique;
  $counting{no_single_alignment_found}                     += $no_aln;
  $counting{unsuitable_sequence_count}                     += $multiple;
  $counting{genomic_sequence_could_not_be_extracted_count} += $no_genomic;

  $counting{total_meCHH_count}                             += $meth_CHH;
  $counting{total_meCHG_count}                             += $meth_CHG;
  $counting{total_meCpG_count}                             += $meth_CpG;
  if ($bowtie2){
    $counting{total_meC_unknown_count}                     += $meth_unknown;
  }

  $counting{total_unmethylated_CHH_count}                  += $unmeth_CHH;
  $counting{total_unmethylated_CHG_count}                  += $unmeth_CHG;
  $counting{total_unmethylated_CpG_count}                  += $unmeth_CpG;
  if ($bowtie2){
    $counting{total_unmethylated_C_unknown_count}          += $unmeth_unknown;
  }

  if ($single_end){
    $counting{CT_CT_count}    += $number_OT;
    $counting{CT_GA_count}    += $number_OB;
    $counting{GA_CT_count}    += $number_CTOT;
    $counting{GA_GA_count}    += $number_CTOB;
  }
  else{
    # paired-end
    $counting{GA_CT_CT_count} += $number_CTOT;
    $counting{CT_GA_CT_count} += $number_OT;
    $counting{GA_CT_GA_count} += $number_CTOB;
    $counting{CT_GA_GA_count} += $number_OB;
  }
}

sub print_final_analysis_report_paired_ends{
	my ($C_to_T_infile_1,$G_to_A_infile_1,$C_to_T_infile_2,$G_to_A_infile_2,$pid,$merge_multi) = @_;

	if ($merge_multi){
		warn "Printing a final merged alignment report for all individual sub-reports\n\n";
	}
	else{
		### All sequences from the original sequence file have been analysed now, therefore deleting temporary C->T or G->A infiles
		if ($directional){
			my $deletion_successful =  unlink "$temp_dir$C_to_T_infile_1","$temp_dir$G_to_A_infile_2";
			if ($deletion_successful == 2){
				warn "\nSuccessfully deleted the temporary files $temp_dir$C_to_T_infile_1 and $temp_dir$G_to_A_infile_2\n\n";
			}
			else{
				warn "Could not delete temporary files $temp_dir$C_to_T_infile_1 and $temp_dir$G_to_A_infile_2 properly: $!\n";
			}
		}
		elsif($pbat){
			# PBAT data should only have 2 files to delete, similar to directional files
			my $deletion_successful =  unlink "$temp_dir$G_to_A_infile_1","$temp_dir$C_to_T_infile_2";
			if ($deletion_successful == 2){
				warn "\nSuccessfully deleted the temporary files $temp_dir$G_to_A_infile_1 and $temp_dir$C_to_T_infile_2\n\n";
			}
			else{
				warn "Could not delete temporary files $temp_dir$G_to_A_infile_1 and $temp_dir$C_to_T_infile_2 properly: $!\n";
			}
		}
		else{ # non-directional
			my $deletion_successful =  unlink "$temp_dir$C_to_T_infile_1","$temp_dir$G_to_A_infile_1","$temp_dir$C_to_T_infile_2","$temp_dir$G_to_A_infile_2";
			if ($deletion_successful == 4){
				warn "\nSuccessfully deleted the temporary files $temp_dir$C_to_T_infile_1, $temp_dir$G_to_A_infile_1, $temp_dir$C_to_T_infile_2 and $temp_dir$G_to_A_infile_2\n\n";
			}
			else{
				warn "Could not delete temporary files properly: $!\n";
			}
		}
	}

	### printing a final report for the alignment procedure
	warn "Final Alignment report\n",'='x22,"\n";
	print REPORT "Final Alignment report\n",'='x22,"\n";
	#  foreach my $index (0..$#fhs){
	#    print "$fhs[$index]->{name}\n";
	#    print "$fhs[$index]->{seen}\talignments on the correct strand in total\n";
	#    print "$fhs[$index]->{wrong_strand}\talignments were discarded (nonsensical alignments)\n\n";
	#  }

	### printing a final report for the methylation call procedure
	warn "Sequence pairs analysed in total:\t$counting{sequences_count}\n";
	print REPORT "Sequence pairs analysed in total:\t$counting{sequences_count}\n";

	my $percent_alignable_sequence_pairs;
	if ($counting{sequences_count} == 0){
		$percent_alignable_sequence_pairs = 0;
	}
	else{
		$percent_alignable_sequence_pairs = sprintf ("%.1f",$counting{unique_best_alignment_count}*100/$counting{sequences_count});
	}	
	print "Number of paired-end alignments with a unique best hit:\t$counting{unique_best_alignment_count}\nMapping efficiency:\t${percent_alignable_sequence_pairs}%\n\n";
	print REPORT "Number of paired-end alignments with a unique best hit:\t$counting{unique_best_alignment_count}\nMapping efficiency:\t${percent_alignable_sequence_pairs}% \n";
	
	print "Sequence pairs with no alignments under any condition:\t$counting{no_single_alignment_found}\n";
	print "Sequence pairs did not map uniquely:\t$counting{unsuitable_sequence_count}\n";
	print "Sequence pairs which were discarded because genomic sequence could not be extracted:\t$counting{genomic_sequence_could_not_be_extracted_count}\n\n";
	print "Number of sequence pairs with unique best (first) alignment came from the bowtie output:\n";
	print join ("\n","CT/GA/CT:\t$counting{CT_GA_CT_count}\t((converted) top strand)","GA/CT/CT:\t$counting{GA_CT_CT_count}\t(complementary to (converted) top strand)","GA/CT/GA:\t$counting{GA_CT_GA_count}\t(complementary to (converted) bottom strand)","CT/GA/GA:\t$counting{CT_GA_GA_count}\t((converted) bottom strand)"),"\n\n";


	print REPORT "Sequence pairs with no alignments under any condition:\t$counting{no_single_alignment_found}\n";
	print REPORT "Sequence pairs did not map uniquely:\t$counting{unsuitable_sequence_count}\n";
	print REPORT "Sequence pairs which were discarded because genomic sequence could not be extracted:\t$counting{genomic_sequence_could_not_be_extracted_count}\n\n";
	print REPORT "Number of sequence pairs with unique best (first) alignment came from the bowtie output:\n";
	print REPORT join ("\n","CT/GA/CT:\t$counting{CT_GA_CT_count}\t((converted) top strand)","GA/CT/CT:\t$counting{GA_CT_CT_count}\t(complementary to (converted) top strand)","GA/CT/GA:\t$counting{GA_CT_GA_count}\t(complementary to (converted) bottom strand)","CT/GA/GA:\t$counting{CT_GA_GA_count}\t((converted) bottom strand)"),"\n\n";
	### detailed information about Cs analysed

	if ($directional){
		print "Number of alignments to (merely theoretical) complementary strands being rejected in total:\t$counting{alignments_rejected_count}\n\n";
		print REPORT "Number of alignments to (merely theoretical) complementary strands being rejected in total:\t$counting{alignments_rejected_count}\n\n";
	}

	warn "Final Cytosine Methylation Report\n",'='x33,"\n";
	print REPORT "Final Cytosine Methylation Report\n",'='x33,"\n";

	my $total_number_of_C = $counting{total_meCHG_count}+ $counting{total_meCHH_count}+$counting{total_meCpG_count}+$counting{total_unmethylated_CHG_count}+$counting{total_unmethylated_CHH_count}+$counting{total_unmethylated_CpG_count};
	warn "Total number of C's analysed:\t$total_number_of_C\n\n";
	warn "Total methylated C's in CpG context:\t$counting{total_meCpG_count}\n";
	warn "Total methylated C's in CHG context:\t$counting{total_meCHG_count}\n";
	warn "Total methylated C's in CHH context:\t$counting{total_meCHH_count}\n";
	warn "Total methylated C's in Unknown context:\t$counting{total_meC_unknown_count}\n\n";
	
	warn "Total unmethylated C's in CpG context:\t$counting{total_unmethylated_CpG_count}\n";
	warn "Total unmethylated C's in CHG context:\t$counting{total_unmethylated_CHG_count}\n";
	warn "Total unmethylated C's in CHH context:\t$counting{total_unmethylated_CHH_count}\n";
	warn "Total unmethylated C's in Unknown context:\t$counting{total_unmethylated_C_unknown_count}\n\n";

	print REPORT "Total number of C's analysed:\t$total_number_of_C\n\n";
	print REPORT "Total methylated C's in CpG context:\t$counting{total_meCpG_count}\n";
	print REPORT "Total methylated C's in CHG context:\t$counting{total_meCHG_count}\n";
	print REPORT "Total methylated C's in CHH context:\t$counting{total_meCHH_count}\n";
	print REPORT "Total methylated C's in Unknown context:\t$counting{total_meC_unknown_count}\n\n";
	  
	print REPORT "Total unmethylated C's in CpG context:\t$counting{total_unmethylated_CpG_count}\n";
	print REPORT "Total unmethylated C's in CHG context:\t$counting{total_unmethylated_CHG_count}\n";
	print REPORT "Total unmethylated C's in CHH context:\t$counting{total_unmethylated_CHH_count}\n";
	print REPORT "Total unmethylated C's in Unknown context:\t$counting{total_unmethylated_C_unknown_count}\n\n";
  
	my $percent_meCHG;
	if (($counting{total_meCHG_count}+$counting{total_unmethylated_CHG_count}) > 0){
		$percent_meCHG = sprintf("%.1f",100*$counting{total_meCHG_count}/($counting{total_meCHG_count}+$counting{total_unmethylated_CHG_count}));
	}

	my $percent_meCHH;
	if (($counting{total_meCHH_count}+$counting{total_unmethylated_CHH_count}) > 0){
		$percent_meCHH = sprintf("%.1f",100*$counting{total_meCHH_count}/($counting{total_meCHH_count}+$counting{total_unmethylated_CHH_count}));
	}

	my $percent_meCpG;
	if (($counting{total_meCpG_count}+$counting{total_unmethylated_CpG_count}) > 0){
		$percent_meCpG = sprintf("%.1f",100*$counting{total_meCpG_count}/($counting{total_meCpG_count}+$counting{total_unmethylated_CpG_count}));
	}

	my $percent_meC_unknown;
	if (($counting{total_meC_unknown_count}+$counting{total_unmethylated_C_unknown_count}) > 0){
		$percent_meC_unknown = sprintf("%.1f",100*$counting{total_meC_unknown_count}/($counting{total_meC_unknown_count}+$counting{total_unmethylated_C_unknown_count}));
	}

	### printing methylated CpG percentage if applicable
	if ($percent_meCpG){
		warn "C methylated in CpG context:\t${percent_meCpG}%\n";
		print REPORT "C methylated in CpG context:\t${percent_meCpG}%\n";
	}
	else{
		warn "Can't determine percentage of methylated Cs in CpG context if value was 0\n";
		print REPORT "Can't determine percentage of methylated Cs in CpG context if value was 0\n";
	}

	### printing methylated C percentage in CHG context if applicable
	if ($percent_meCHG){
		warn "C methylated in CHG context:\t${percent_meCHG}%\n";
		print REPORT "C methylated in CHG context:\t${percent_meCHG}%\n";
	}
	else{
		warn "Can't determine percentage of methylated Cs in CHG context if value was 0\n";
		print REPORT "Can't determine percentage of methylated Cs in CHG context if value was 0\n";
	}

	### printing methylated C percentage in CHH context if applicable
	if ($percent_meCHH){
		warn "C methylated in CHH context:\t${percent_meCHH}%\n";
		print REPORT "C methylated in CHH context:\t${percent_meCHH}%\n";
	}
	else{
		warn "Can't determine percentage of methylated Cs in CHH context if value was 0\n";
		print REPORT "Can't determine percentage of methylated Cs in CHH context if value was 0\n";
	}

	### printing methylated C percentage (Unknown C context) if applicable
	if ($percent_meC_unknown){	
		warn "C methylated in unknown context (CN or CHN):\t${percent_meC_unknown}%\n";
		print REPORT "C methylated in unknown context (CN or CHN):\t${percent_meC_unknown}%\n";
	}
	else{
		warn "Can't determine percentage of methylated Cs in unknown context (CN or CHN) if value was 0\n";
		print REPORT "Can't determine percentage of methylated Cs in unknown context (CN or CHN) if value was 0\n";
	}

	print REPORT "\n\n";
	warn "\n\n";

}
















