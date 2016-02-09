#!/usr/bin/perl

use Getopt::Long;
use FileHandle;
use DBI();

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Import CosmicCompleteExport\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_drop = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS cosmic_complete_export});
my $dbq_create = $dbh_titanic->prepare(q{CREATE TABLE cosmic_complete_export (
											complete_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
											gene_name VARCHAR(48) DEFAULT '' NOT NULL,
											accession_number VARCHAR(32) DEFAULT '' NOT NULL,
                                            gene_cds_length INT UNSIGNED,
											hgnc_id INT UNSIGNED,
											sample_name VARCHAR(32) DEFAULT '',
											id_sample INT UNSIGNED,
											id_tumour INT UNSIGNED,
											primary_site VARCHAR(64) DEFAULT '',
											site_subtype VARCHAR(64) DEFAULT '',
											primary_histology VARCHAR(128) DEFAULT '',
											histology_subtype VARCHAR(128) DEFAULT '',
											genome_wide_screen CHAR DEFAULT '',
											mutation_id INT UNSIGNED,
											mutation_cds VARCHAR(256) DEFAULT '',
											mutation_aa VARCHAR(256) DEFAULT '',
											mutation_description ENUM('Complex - compound substitution','Complex - deletion inframe',
												'Complex - frameshift','Complex - insertion inframe','Deletion - Frameshift','Deletion - In frame',
												'Insertion - Frameshift','Insertion - In frame','No detectable mRNA/protein',
												'Nonstop extension','Substitution - Missense','Substitution - Nonsense','Substitution - coding silent',
												'Unknown','Whole gene deletion'),
											mutation_zygosity ENUM('het','hom'),
											mutation_ncbi36_genome_position VARCHAR(32) DEFAULT '',
											mutation_ncbi36_strand CHAR DEFAULT '',
											mutation_grch37_genome_position VARCHAR(32) DEFAULT '',
											mutation_grch37_strand CHAR DEFAULT '',
											mutation_somatic_status ENUM('Confirmed germline variant','Confirmed somatic variant','Not specified',
														'Reported in another cancer sample as somatic','Reported in another sample as germline',
														'Variant of unknown origin'),
											pubmed_pmid INT UNSIGNED,
											sample_source VARCHAR(32) DEFAULT '',
											tumour_origin ENUM('NS','adenoma adjacent to primary tumour','hyperplasia adjacent to primary tumour',
														'metastasis','primary','recurrent','secondary','surgery fresh/frozen'),
											comments TEXT,
											PRIMARY KEY (complete_id),
                                            KEY (gene_name),
                                            KEY (mutation_aa),
                                            KEY (mutation_description))});
my $dbq_insert = $dbh_titanic->prepare(q{INSERT INTO cosmic_complete_export VALUES(NULL,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)});

# Drop and create table

$dbq_drop->execute;
$dbq_create->execute;

# Read CosmicCompleteExport file

my $cosmic_fn = $mokcatanic_data_root . "/CosmicCompleteExport.tsv.gz";

open(COSM, "gzip -dc < $cosmic_fn |") || die "\n+++ Cannot open $cosmic_fn to read complete export: $!\n";

<COSM>; # Skip first line
while (<COSM>) {
	my $gene_name, $accession_number, $gene_cds_length, $hgnc_id, $sample_name, $id_sample, $id_tumour, $primary_site, $site_subtype, $primary_histology, $histology_subtype, $genome_wide_screen,
		$mutation_id, $mutation_cds, $mutation_aa, $mutation_description, $mutation_zygosity, $mutation_ncbi36_genome_position, $mutation_ncbi36_strand, $mutation_grch37_genome_position, $mutation_grch37_strand,
		$mutation_somatic_status, $pubmed_pmid, $sample_source, $tumour_origin, $comments;
	
	chomp;
	( $gene_name, $accession_number, $gene_cds_length, $hgnc_id, $sample_name, $id_sample, $id_tumour, $primary_site, $site_subtype, $primary_histology, $histology_subtype, $genome_wide_screen,
	$mutation_id, $mutation_cds, $mutation_aa, $mutation_description, $mutation_zygosity, $mutation_ncbi36_genome_position, $mutation_ncbi36_strand, $mutation_grch37_genome_position, $mutation_grch37_strand,
	$mutation_somatic_status, $pubmed_pmid, $sample_source, $tumour_origin, $comments ) = split(/\t/);
	
	$dbq_insert->execute($gene_name, $accession_number, $gene_cds_length, $hgnc_id, $sample_name, $id_sample, $id_tumour, $primary_site, $site_subtype, $primary_histology, $histology_subtype, $genome_wide_screen,
	$mutation_id, $mutation_cds, $mutation_aa, $mutation_description, $mutation_zygosity, $mutation_ncbi36_genome_position, $mutation_ncbi36_strand, $mutation_grch37_genome_position, $mutation_grch37_strand,
	$mutation_somatic_status, $pubmed_pmid, $sample_source, $tumour_origin, $comments);
    $dbq_insert->finish;
    
	# Pointless ticker
	
	if ($tick_c % $tick_step == 0) {
		print STDERR $tick_char;
		if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
			print STDERR "\n";
		}
	}
	$tick_c++;
}

close(COSM) || die "\n+++ Cannot close $cosmic_fn to stop reading complete export: $!\n";

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
