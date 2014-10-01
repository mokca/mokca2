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

print STDERR "::: Import CosmicHGNC - list of genes and cross references\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_drop = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS cosmic_hgnc});
my $dbq_create = $dbh_titanic->prepare(q{CREATE TABLE cosmic_hgnc (
											cosmic_id INT UNSIGNED NOT NULL,
											cosmic_gene_name VARCHAR(64) DEFAULT '' NOT NULL,
											entrez_id INT UNSIGNED,
											hgnc_id INT UNSIGNED,
											mutated BOOL DEFAULT FALSE,
											cancer_census BOOL DEFAULT FALSE,
											PRIMARY KEY (cosmic_id)
	)});
my $dbq_insert = $dbh_titanic->prepare(q{INSERT INTO cosmic_hgnc VALUES(?,?,?,?,?,?)});

# Drop and create table cosmic_hgnc

$dbq_drop->execute;
$dbq_create->execute;

# Open CosmicHGNC file

my $cosmic_fn = $mokcatanic_data_root . "/CosmicHGNC_" . $cosmic_release . ".tsv.gz";

open(HGNC, "gzip -dc < $cosmic_fn |") || die "\n+++ Cannot open $cosmic_fn to read gene data: $!\n";

<HGNC>;	# Skip first line
while (<HGNC>) {
	my $cosmic_id, $cosmic_gene_name, $entrez_id, $hgnc_id, $mutated, $cancer_census, $mutated_bool, $cancer_census_bool;
	
	chomp;
	( $cosmic_id, $cosmic_gene_name, $entrez_id, $hgnc_id, $mutated, $cancer_census) = split(/\t/);
	
	# Translate y/n to bool for Mutated? and Cancer Census?
	
	if ($cosmic_id && $mutated && $cancer_census) {
		if ($mutated eq 'y') {
			$mutated_bool = 1;
		} elsif ($mutated eq 'n') {
			$mutated_bool = 0;
		} else {
			die "+++ Unexpected value for 'mutated': $mutated\n";
		}
		
		if ($cancer_census eq 'y') {
			$cancer_census_bool = 1;
		} elsif ($cancer_census eq 'n') {
			$cancer_census_bool = 0;
		} else {
			die "+++ Unexpected value for 'cancer census': $cancer_census\n";
		}
				
		$dbq_insert->execute($cosmic_id, $cosmic_gene_name, $entrez_id, $hgnc_id, $mutated_bool, $cancer_census_bool);
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

close(HGNC) || die "\n+++ Cannot close $cosmic_fn to stop reading gene data: $!\n";

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
