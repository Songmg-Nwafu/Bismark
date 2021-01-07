#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;
$|++;

#This script is supposed to convert a specified reference genome into two different bisulfite
#converted versions and index them for alignments with Bowtie 2 (default), or HISAT2. The first
#bisulfite genome will have all Cs converted to Ts (C->T), and the other one will have all Gs
#converted to As (G->A). Both bisulfite genomes will be stored in subfolders within the reference
#genome folder. 


my $verbose;
my $help;
my $version;
my $man;
my $parallel;
my $multi_fasta;
my $single_fasta;
my $slam;
my $large_index;

my $parent_dir = getcwd();

my $genomic_composition;
my %genomic_freqs; # storing the genomic sequence composition
my %freqs;

my $bismark_version = 'v0.1.1';
my $last_modified = "7 September 2020";

my $command_line = GetOptions ('verbose'         => \$verbose,
			'help'                => \$help,
			'man'                 => \$man,
			'version'             => \$version,
			'slam'                => \$slam,
			'parallel:i'          => \$parallel,
			'single_fasta'        => \$single_fasta,
			'genomic_composition' => \$genomic_composition,
			'large-index'         => \$large_index,
	);
 
die "Please respecify command line options\n\n" unless ($command_line);


if ($help or $man){
  print_helpfile();
  exit;
}

if ($version){
  print << "VERSION";

          Bismark - Bisulfite Mapper and Methylation Caller.

          Bismark Genome Preparation Version: $bismark_version
        Copyright 2010-19 Felix Krueger, Babraham Bioinformatics
              www.bioinformatics.babraham.ac.uk/projects/
               https://github.com/FelixKrueger/Bismark
			  

VERSION
  exit;
}

my $genome_folder = shift @ARGV; # mandatory
my %chromosomes; # checking if chromosome names are unique (required)

# Ensuring a genome folder has been specified
if ($genome_folder){
    unless ($genome_folder =~ /\/$/){
	$genome_folder =~ s/$/\//;
    }
    $verbose and print "Path to genome folder specified as: $genome_folder\n";
    chdir $genome_folder or die "Could't move to directory $genome_folder. Make sure the directory exists! $!";
    
    # making the genome folder path abolsolute so it won't break if the path was specified relative
    $genome_folder = getcwd();	
    unless ($genome_folder =~ /\/$/){
	$genome_folder =~ s/$/\//;
    }
}
else{
    die "Please specify a genome folder to be used for bisulfite conversion\n\n";
}

if (defined $parallel){
    die "--parallel should have a value of 2 or more. Please respecify\n\n" unless ($parallel > 1);

    warn "Using $parallel threads for the top and bottom strand indexing processes each, so using ",$parallel*2," cores in total\n"; sleep(1);
}
if ($slam){
	warn "This is a SLAM-seq test version. Proceed with care!\n\n"; sleep(3);
}	

my $CT_dir;
my $GA_dir;


if ($single_fasta){
    warn "Writing individual genomes out into single-entry fasta files (one per chromosome)\n\n";
    $multi_fasta = 0;
}
else{
    warn "Writing bisulfite genomes out into a single MFA (multi FastA) file\n\n";
    $single_fasta = 0;
    $multi_fasta = 1;
}

if ($slam){
	warn "Genome will be generated and indexed in with in-silico T->C transitions, and NOT in BISULFITE MODE\n"; sleep(3);
}

if($large_index){
	warn "Large-index specified. Forcing generated index to be 'large', even if reference has fewer than 4 billion nucleotides.\n\n"; sleep(1);
	$large_index = '--large-index';
}
else{	
	$large_index = '';
}

my @filenames = create_bisulfite_genome_folders();

if ($genomic_composition){
    get_genomic_frequencies();
    warn "Finished processing genomic nucleotide frequencies\n\n";
    %chromosomes = (); # resetting
}

process_sequence_files ();

sub process_sequence_files {
    
    
    my ($total_CT_conversions,$total_GA_conversions) = (0,0);
    $verbose and print "Bismark Genome Preparation - Step II: Bisulfite converting reference genome\n\n";
    sleep (1);
    
    chdir $genome_folder or die "Could't move to directory $genome_folder $!";
    
    $verbose and print "conversions performed:\n";
    if ($slam){
		$verbose and print join("\t",'chromosome','T->C','A->G'),"\n";
	}
	else{
		$verbose and print join("\t",'chromosome','C->T','G->A'),"\n";
	}	
    
  
    ### If someone wants to index a genome which consists of thousands of contig and scaffold files we need to write the genome conversions into an MFA file
    ### Otherwise the list of comma separated chromosomes we provide for bowtie-build will get too long for the kernel to handle
    ### This is now the default option
    
    if ($multi_fasta){
		### Here we just use one multi FastA file name, append .CT_conversion or .GA_conversion and print all sequence conversions into these files
		my $bisulfite_CT_conversion_filename = "$CT_dir/genome_mfa.CT_conversion.fa";
		open (CT_CONVERT,'>',$bisulfite_CT_conversion_filename) or die "Can't write to file $bisulfite_CT_conversion_filename: $!\n";
	
		my $bisulfite_GA_conversion_filename = "$GA_dir/genome_mfa.GA_conversion.fa";
		open (GA_CONVERT,'>',$bisulfite_GA_conversion_filename) or die "Can't write to file $bisulfite_GA_conversion_filename: $!\n";
    }
    
    foreach my $filename(@filenames){
	my ($chromosome_CT_conversions,$chromosome_GA_conversions) = (0,0);
	if ($filename =~ /\.gz$/){
	    open (IN,"gunzip -c $filename |") or die "Failed to read from sequence file $filename $!\n";
	}
	else{
	    open (IN,$filename) or die "Failed to read from sequence file $filename $!\n";
	}
	# warn "Reading chromosome information from $filename\n\n";
	
	### first line needs to be a fastA header
	my $first_line = <IN>;
	chomp $first_line;
	
	### Extracting chromosome name from the FastA header
	my $chromosome_name = extract_chromosome_name($first_line);
	
	### Exiting if a chromosome with the same name was present already
	if (exists $chromosomes{$chromosome_name}){
	    die "Exiting because chromosome name '$chromosome_name' already exists. Please make sure all chromosomes have a unique name!\n";
	}
	else{
	    $chromosomes{$chromosome_name}++;
	}
	
	### alternatively, chromosomes can be written out into single-entry FastA files. This will only work for genomes with up to a few hundred chromosomes.
	unless ($multi_fasta){
	    my $bisulfite_CT_conversion_filename = "$CT_dir/$chromosome_name";
	    $bisulfite_CT_conversion_filename =~ s/$/.CT_conversion.fa/;
	    open (CT_CONVERT,'>',$bisulfite_CT_conversion_filename) or die "Can't write to file $bisulfite_CT_conversion_filename: $!\n";
	    
	    my $bisulfite_GA_conversion_filename = "$GA_dir/$chromosome_name";
	    $bisulfite_GA_conversion_filename =~ s/$/.GA_conversion.fa/;
	    open (GA_CONVERT,'>',$bisulfite_GA_conversion_filename) or die "Can't write to file $bisulfite_GA_conversion_filename: $!\n";
	}
	
	print CT_CONVERT ">",$chromosome_name,"_CT_converted\n"; # first entry
	print GA_CONVERT ">",$chromosome_name,"_GA_converted\n"; # first entry
	### TODO: Change this for GrandSlam
	
	while (<IN>){
	    
      ### in case the line is a new fastA header
      if ($_ =~ /^>/){
		### printing out the stats for the previous chromosome
		$verbose and print join ("\t",$chromosome_name,$chromosome_CT_conversions,$chromosome_GA_conversions),"\n";
		### resetting the chromosome transliteration counters
		($chromosome_CT_conversions,$chromosome_GA_conversions) = (0,0);
		
		### Extracting chromosome name from the additional FastA header
		$chromosome_name = extract_chromosome_name($_);

		### alternatively, chromosomes can be written out into single-entry FastA files. This will only work for genomes with up to a few hundred chromosomes.
		unless ($multi_fasta){
		  my $bisulfite_CT_conversion_filename = "$CT_dir/$chromosome_name";
		  $bisulfite_CT_conversion_filename =~ s/$/.CT_conversion.fa/;
		  open (CT_CONVERT,'>',$bisulfite_CT_conversion_filename) or die "Can't write to file $bisulfite_CT_conversion_filename: $!\n";
		
		  my $bisulfite_GA_conversion_filename = "$GA_dir/$chromosome_name";
		  $bisulfite_GA_conversion_filename =~ s/$/.GA_conversion.fa/;
		  open (GA_CONVERT,'>',$bisulfite_GA_conversion_filename) or die "Can't write to file $bisulfite_GA_conversion_filename: $!\n";
		}

		print CT_CONVERT ">",$chromosome_name,"_CT_converted\n";
		print GA_CONVERT ">",$chromosome_name,"_GA_converted\n";
      }

		else{
			my $sequence = uc$_;

			### (I) First replacing all ambiguous sequence characters (such as M,S,R....) by N (G,A,T,C,N and the line endings \r and \n are added to a character group)
			
			$sequence =~ s/[^ATCGN\n\r]/N/g;
			
			### (II) Writing the chromosome out into a C->T converted version (equals forward strand conversion)
			my $CT_sequence = $sequence;
			my $CT_transliterations_performed;
			if ($slam){
				$CT_transliterations_performed = ($CT_sequence =~ tr/T/C/); # converts all Ts into Cs for SLAM-Seq
			}
			else{
				$CT_transliterations_performed = ($CT_sequence =~ tr/C/T/); # converts all Cs into Ts
			}
			$total_CT_conversions += $CT_transliterations_performed;
			$chromosome_CT_conversions += $CT_transliterations_performed;
			
			print CT_CONVERT $CT_sequence;
			
			### (III) Writing the chromosome out in a G->A converted version of the forward strand (this is equivalent to reverse-
			### complementing the forward strand and then C->T converting it)
			
			my $GA_sequence = $sequence;
			my $GA_transliterations_performed;
			if ($slam){
				$GA_transliterations_performed = ($GA_sequence =~ tr/A/G/); # converts all As to Gs on the forward strand, for SLAM-Seq
			}
			else{
				$GA_transliterations_performed = ($GA_sequence =~ tr/G/A/); # converts all Gs to As on the forward strand
			}	
			$total_GA_conversions += $GA_transliterations_performed;
			$chromosome_GA_conversions += $GA_transliterations_performed;
			
			print GA_CONVERT $GA_sequence;
			
		}
    }
    $verbose and print join ("\t",$chromosome_name,$chromosome_CT_conversions,$chromosome_GA_conversions),"\n";
  }
  
 close (CT_CONVERT) or die "Failed to close filehandle: $!\n";
 close (GA_CONVERT) or die "Failed to close filehandle: $!\n";


  print "\nTotal number of conversions performed:\n";
	if ($slam){
		print "T->C:\t$total_CT_conversions\n";
		print "A->G:\t$total_GA_conversions\n";
		warn "\nStep II - Genome SLAM conversions - completed\n\n\n";
	}	
	else{
		print "C->T:\t$total_CT_conversions\n";
		print "G->A:\t$total_GA_conversions\n";
		warn "\nStep II - Genome bisulfite conversions - completed\n\n\n";
	}
  
}

sub get_genomic_frequencies{
    
    warn "Calculating genomic frequencies (this may take several minutes depending on genome size) ...\n";
    warn "="x164,"\n";
    read_genome_into_memory($parent_dir);
    foreach my $chr (keys %chromosomes){
	warn "Processing chromosome >> $chr <<\n";
	process_sequence($chromosomes{$chr});
    }
    
    %genomic_freqs = %freqs;
    
    ### Attempting to write a genomic nucleotide frequency table out to the genome folder so we can re-use it next time without the need to re-calculate
    ### if this fails
    if  ( open (FREQS,'>',"${genome_folder}genomic_nucleotide_frequencies.txt") ){
	warn "Writing genomic nucleotide frequencies to the file >${genome_folder}genomic_nucleotide_frequencies.txt< for future re-use\n";
	foreach my $f(sort keys %genomic_freqs){
	    warn "Writing count of (di-)nucleotide: $f\t$genomic_freqs{$f}\n";
	    print FREQS "$f\t$genomic_freqs{$f}\n";
	}
	close FREQS or warn "Failed to close filehandle FREQS: $!\n\n";
    }
    else{
	warn "Failed to write out file ${genome_folder}genomic_nucleotide_frequencies.txt because of: $!. Skipping writing out genomic frequency table\n\n";
    }
}


sub process_sequence{
    
    my $seq = shift;
    my $mono;
    my $di;

    foreach my $index (0..(length$seq)-1){
		my $counted = 0;
		if ($index%10000000==0){
			# warn "Current index number is $index\n";
		}
		
		$mono = substr($seq,$index,1);
		unless ( $mono eq 'N'){
			$freqs{$mono}++;
		}
		
		unless ( ($index + 2) > length$seq ){
			$di = substr($seq,$index,2);
			if (index($di,'N') < 0) {
				$freqs{$di}++;
			}
		}
    }
}

sub extract_chromosome_name {
	## Bowtie extracts the first string after the inition > in the FASTA file, so we are doing this as well
	my $fasta_header = shift;
	if ($fasta_header =~ s/^>//){
		my ($chromosome_name) = split (/\s+/,$fasta_header);
			return $chromosome_name;
		}
	else{
		die "The specified chromosome ($fasta_header) file doesn't seem to be in FASTA format as required!\n";
	}
}


sub create_bisulfite_genome_folders{

	$verbose and print "Bismark Genome Preparation - Step I: Preparing folders\n\n";

	chdir $genome_folder or die "Could't move to directory $genome_folder. Make sure the directory exists! $!";


	# Exiting unless there are fastA files in the folder
	my @filenames = <*.fa>;
  
	### if there aren't any genomic files with the extension .fa we will look for files with the extension .fa.gz
	unless (@filenames){
		@filenames =  <*.fa.gz>;
	}
	### if there aren't any genomic files with the extension .fa or .fa.gz we will look for files with the extension .fasta
	unless (@filenames){
		@filenames =  <*.fasta>;
	}
	### if there aren't any genomic files with the extension .fa, .fa.gz or .fasta we will look for files with the extension .fasta.gz
	unless (@filenames){
		@filenames =  <*.fasta>;
	}
	unless (@filenames){
		die "The specified genome folder $genome_folder does not contain any sequence files in FastA format (with .fa, .fa.gz, .fasta or .fasta.gz file extensions)\n";
	}

	warn "Bisulfite Genome Indexer version $bismark_version (last modified: $last_modified)\n";
	sleep (1);

	# creating a directory inside the genome folder to store the bisfulfite genomes unless it already exists
	my $bisulfite_dir = "${genome_folder}Bisulfite_Genome/";
	unless (-d $bisulfite_dir){
		mkdir $bisulfite_dir or die "Unable to create directory $bisulfite_dir $!\n";
		$verbose and print "Created Bisulfite Genome folder $bisulfite_dir\n";
	}
	else{
		print "\nA directory called $bisulfite_dir already exists. Already existing converted sequences and/or already existing Bowtie 2 or HISAT2) indices will be overwritten!\n\n";
	}

	chdir $bisulfite_dir or die "Unable to move to $bisulfite_dir\n";
	$CT_dir = "${bisulfite_dir}CT_conversion/";
	$GA_dir = "${bisulfite_dir}GA_conversion/";

	# creating 2 subdirectories to store a C->T (forward strand conversion) and a G->A (reverse strand conversion)
	# converted version of the genome
	unless (-d $CT_dir){
		mkdir $CT_dir or die "Unable to create directory $CT_dir $!\n";
		$verbose and print "Created Bisulfite Genome folder $CT_dir\n";
	}
	unless (-d $GA_dir){
		mkdir $GA_dir or die "Unable to create directory $GA_dir $!\n";
		$verbose and print "Created Bisulfite Genome folder $GA_dir\n";
	}

	# moving back to the original genome folder
	chdir $genome_folder or die "Could't move to directory $genome_folder $!";
	# $verbose and print "Moved back to genome folder folder $genome_folder\n";
	warn "\nStep I - Prepare genome folders - completed\n\n\n";
	
	return @filenames;
	
}

sub read_genome_into_memory{

  ## reading in and storing the specified genome in the %chromosomes hash
  chdir ($genome_folder) or die "Can't move to $genome_folder: $!";
  warn "Now reading in and storing sequence information of the genome specified in: $genome_folder\n\n";

  my @chromosome_filenames =  <*.fa>;
  
  ### if there aren't any genomic files with the extension .fa we will look for files with the extension .fa.gz
  unless (@chromosome_filenames){
      @chromosome_filenames =  <*.fa.gz>;
  }
  
  ### if there aren't any genomic files with the extension .fa or .fq.gz we will look for files with the extension .fasta
  unless (@chromosome_filenames){
      @chromosome_filenames =  <*.fasta>;
  }
  ### if there aren't any genomic files with the extension .fa or .fa.gz or .fasta we will look for files with the extension .fasta.gz
  unless (@chromosome_filenames){
      @chromosome_filenames =  <*.fasta.gz>;
  }

  unless (@chromosome_filenames){
      die "The specified genome folder $genome_folder does not contain any sequence files in FastA format (with .fa, .fa.gz, .fasta or .fasta.gz file extensions)\n";
  }

  foreach my $chromosome_filename (@chromosome_filenames){
      
      # skipping the tophat entire mouse genome fasta file
      next if ($chromosome_filename eq 'Mus_musculus.NCBIM37.fa');
      if( $chromosome_filename =~ /\.gz$/){
	  open (CHR_IN,"gunzip -c $chromosome_filename |") or die "Failed to read from sequence file $chromosome_filename $!\n";	
      }
      else{
	  open (CHR_IN,$chromosome_filename) or die "Failed to read from sequence file $chromosome_filename $!\n";
      }
      
      ### first line needs to be a fastA header
      my $first_line = <CHR_IN>;
      chomp $first_line;
      $first_line =~ s/\r//; # removing /r carriage returns
      
      ### Extracting chromosome name from the FastA header
      my $chromosome_name = extract_chromosome_name($first_line);
      
      my $sequence;
      while (<CHR_IN>){
	  chomp;
	  $_ =~ s/\r//; # removing /r carriage returns
	  if ($_ =~ /^>/){
	      ### storing the previous chromosome in the %chromosomes hash, only relevant for Multi-Fasta-Files (MFA)
	      if (exists $chromosomes{$chromosome_name}){
		  warn "chr $chromosome_name (",length $sequence ," bp)\n";
		  die "Exiting because chromosome name '$chromosome_name' already exists. Please make sure all chromosomes have a unique name!\n";
	      }
	      else {
		  if (length($sequence) == 0){
		      warn "Chromosome $chromosome_name in the multi-fasta file $chromosome_filename did not contain any sequence information!\n";
		  }
		  warn "chr $chromosome_name (",length $sequence ," bp)\n";
	  $chromosomes{$chromosome_name} = $sequence;
	      }
	      ### resetting the sequence variable
	      $sequence = '';
	      ### setting new chromosome name
	      $chromosome_name = extract_chromosome_name($_);
	  }
	  else{
	      $sequence .= uc$_;
	  }
      }

    if (exists $chromosomes{$chromosome_name}){
      warn "chr $chromosome_name (",length $sequence ," bp)\t";
      die "Exiting because chromosome name '$chromosome_name' already exists. Please make sure all chromosomes have a unique name.\n";
    }
    else{
      if (length($sequence) == 0){
	warn "Chromosome $chromosome_name in the file $chromosome_filename did not contain any sequence information!\n";
      }
      warn "chr $chromosome_name (",length $sequence ," bp)\n";
      $chromosomes{$chromosome_name} = $sequence;
    }
  }
  warn "\n";
  chdir $parent_dir or die "Failed to move to directory $parent_dir\n";
}




sub print_helpfile{
  print << 'HOW_TO';


DESCRIPTION

This script is supposed to convert a specified reference genome into two different bisulfite
converted versions and index them for alignments with Bowtie 2 (default), or HISAT2. The first
bisulfite genome will have all Cs converted to Ts (C->T), and the other one will have all Gs
converted to As (G->A). Both bisulfite genomes will be stored in subfolders within the reference
genome folder. 

The new, and still experimental, --slam mode will produce T->C and A->G converted genomes instead.
The structure of the genome folder will remain the same as for BS-Seq data. This might, or might not,
change in a future release.

The following is a brief description of command line options and arguments to control the
Bismark Genome Preparation:


USAGE: bismark_genome_preparation [options] <arguments>


OPTIONS:

--help                   Displays this help file and exits.

--version                Displays version information and exits.

--verbose                Print verbose output for more details or debugging.

--parallel INT           Use several threads for each indexing process to speed up the genome preparation step.
                         Remember that the indexing is run twice in parallel already (for the top and bottom strand
                         separately), so e.g. '--parallel 4' will use 8 threads in total. Please also see --large-index
                         for parallel processing of VERY LARGE genomes (e.g. the axolotl)

--single_fasta           Instruct the Bismark Indexer to write the converted genomes into
                         single-entry FastA files instead of making one multi-FastA file (MFA)
                         per chromosome. This might be useful if individual bisulfite converted
                         chromosomes are needed (e.g. for debugging), however it can cause a
                         problem with indexing if the number of chromosomes is vast (this is likely
                         to be in the range of several thousand files; the operating system can
                         only handle lists up to a certain length, and some newly assembled
                         genomes may contain 20000-500000 contigs of scaffold files which do exceed
                         this list length limit).

--genomic_composition    Calculate and extract the genomic sequence composition for mono and di-nucleotides
                         and write the genomic composition table 'genomic_nucleotide_frequencies.txt' to the
                         genome folder. This may be useful later on when using bam2nuc or the Bismark option
                         --nucleotide_coverage.
						 
						 
--slam                   Instead of performing an in-silico bisulfite conversion, this mode transforms T to C (forward strand),
                         or A to G (reverse strand). The folder structure and rest of the indexing process is currently
                         exactly the same as for bisulfite sequences, but this might change at some point. This means
                         that a genome prepared in --slam mode is currently indistinguishable from a true Bisulfite Genome,
                         so please make sure you name the genome folder appropriately to avoid confusion.

--large-index            Force generated index to be 'large', even if reference has fewer than 4 billion nucleotides. At the time
                         of writing this is required for parallel processing of VERY LARGE genomes (e.g. the axolotl)

ARGUMENTS:

<path_to_genome_folder>  The path to the folder containing the genome to be bisulfite converted. The Bismark Genome
                         Preparation expects one or more fastA files in the folder (with the file extension: .fa or
                         .fasta (also ending in .gz)). Specifying this path is mandatory.


This script was last modified on 7 September 2020.

HOW_TO
}
