#!/usr/bin/perl

use Getopt::Long;
use FileHandle;
use DBI();

# Autoflush STDERR

STDERR->autoflush(1);

# Setup universal paths, database details

do "path_setup.pl";
do "database_setup.pl";

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Fetch COSMIC fasta files for protein and cDNA sequences\n";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_cosmic_genes = $dbh_titanic->prepare(q{SELECT cosmic_gene_name FROM cosmic_hgnc});

# Retrieve gene names and fetch files for each name

my $cosmic_gene_name;
my $protein_data_root = $mokcatanic_data_root . "/fasta/cosmic/protein/";
my $cdna_data_root = $mokcatanic_data_root . "/fasta/cosmic/cdna/";


$dbq_cosmic_genes->execute;
while (( $cosmic_gene_name) = $dbq_cosmic_genes->fetchrow_array) {
	
	my $protein_url, $protein_fn;
	my $cdna_url, $cdna_fn;
	
	if (substr($cosmic_gene_name,0,1) =~ /[0-9]/) {
		$protein_url = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/fasta_files/0-9/" . $cosmic_gene_name . "_protein.txt";
		$cdna_url = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/fasta_files/0-9/" . $cosmic_gene_name . "_cdna.txt";
	} else {
		$protein_url = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/fasta_files/" . uc(substr($cosmic_gene_name, 0, 1)) . "/" . $cosmic_gene_name . "_protein.txt";
		$cdna_url = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/fasta_files/" . uc(substr($cosmic_gene_name, 0, 1)) . "/" . $cosmic_gene_name . "_cdna.txt";
	}
		
	$protein_fn = $protein_data_root . $cosmic_gene_name . "_protein.txt";
	$cdna_fn = $cdna_data_root . $cosmic_gene_name . "_cdna.txt";

	if (! -e $protein_fn) {
		system("wget", $protein_url, "-q", "-N", "-P", $protein_data_root) == 0 || die "\n+++Cannot fetch $protein_url:$?\n";
	}
	if (! -e $cdna_fn) {
		system("wget", $cdna_url, "-q", "-N", "-P", $cdna_data_root) == 0 || die "\n+++Cannot fetch $cdna_url:$?\n";
	}
	
	# Pointless ticker
	
	if ($tick_c % $tick_step == 0) {
		print STDERR $tick_char;
		if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
			print STDERR "\n";
		}
	}
	$tick_c++;
}
$dbq_cosmic_genes->finish;

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
