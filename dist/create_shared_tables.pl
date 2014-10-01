#!/usr/bin/perl

use Getopt::Long;
use FileHandle;
use DBI();

# Autoflush STDERR

STDERR->autoflush(1);

# Explain purpose

print STDERR "::: Create shared database tables\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_drop_uniprot = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot});
my $dbq_create_uniprot = $dbh_titanic->prepare(q{CREATE TABLE uniprot(
	uniprot_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
	dataset ENUM('Swiss-Prot','TrEMBL'),
    accession VARCHAR(10),
	name VARCHAR(12),
    recommended_name VARCHAR(255),
    alternative_name_count INT UNSIGNED,
    submitted_name_count INT UNSIGNED,
    gene_name VARCHAR(16),
    gene_synonym_count INT UNSIGNED,
    refseq_count INT UNSIGNED,
    pdb_count INT UNSIGNED,
    mim_count INT UNSIGNED,
    pid_count INT UNSIGNED,
    reactome_count INT UNSIGNED,
    drugbank_count INT UNSIGNED,
    ensembl_count INT UNSIGNED,
    go_count INT UNSIGNED,
	hgnc_numeric INT UNSIGNED,
	hgnc_name VARCHAR(64),
	isoform_count INT UNSIGNED,
    kegg VARCHAR(16),
    geneid VARCHAR(16),
    chembl VARCHAR(16),
	PRIMARY KEY (uniprot_id),
    KEY (accession))});

my $dbq_drop_uniprot_alternative_name = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_alternative_name});
my $dbq_create_uniprot_alternative_name = $dbh_titanic->prepare(q{CREATE TABLE uniprot_alternative_name(
    uniprot_id INT UNSIGNED NOT NULL,
    alternative_name VARCHAR(255),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_submitted_name = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_submitted_name});
my $dbq_create_uniprot_submitted_name = $dbh_titanic->prepare(q{CREATE TABLE uniprot_submitted_name(
    uniprot_id INT UNSIGNED NOT NULL,
    submitted_name VARCHAR(255),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_gene_synonym = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_gene_synonym});
my $dbq_create_uniprot_gene_synonym = $dbh_titanic->prepare(q{CREATE TABLE uniprot_gene_synonym(
    uniprot_id INT UNSIGNED NOT NULL,
    gene_synonym VARCHAR(80),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_refseq = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_refseq });
my $dbq_create_uniprot_refseq = $dbh_titanic->prepare(q{CREATE TABLE uniprot_refseq(
    uniprot_id INT UNSIGNED NOT NULL,
    refseq_id VARCHAR(14),
    nucleotide_id VARCHAR(14),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_pdb = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_pdb });
my $dbq_create_uniprot_pdb = $dbh_titanic->prepare(q{CREATE TABLE uniprot_pdb(
    uniprot_id INT UNSIGNED NOT NULL,
    pdb_id VARCHAR(4),
    method VARCHAR(5),
    resolution VARCHAR(6),
    chains VARCHAR(64),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_mim = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_mim });
my $dbq_create_uniprot_mim = $dbh_titanic->prepare(q{CREATE TABLE uniprot_mim(
    uniprot_id INT UNSIGNED NOT NULL,
    mim_id INT UNSIGNED,
    type ENUM('phenotype','gene'),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_pid = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_pid });
my $dbq_create_uniprot_pid = $dbh_titanic->prepare(q{CREATE TABLE uniprot_pid(
    uniprot_id INT UNSIGNED NOT NULL,
    pid_id VARCHAR(32),
    name VARCHAR(128),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_reactome = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_reactome });
my $dbq_create_uniprot_reactome = $dbh_titanic->prepare(q{CREATE TABLE uniprot_reactome(
    uniprot_id INT UNSIGNED NOT NULL,
    reactome_id VARCHAR(12),
    name VARCHAR(64),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_drugbank = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_drugbank });
my $dbq_create_uniprot_drugbank = $dbh_titanic->prepare(q{CREATE TABLE uniprot_drugbank(
    uniprot_id INT UNSIGNED NOT NULL,
    drugbank_id VARCHAR(7),
    name VARCHAR(32),
    KEY (uniprot_id))});

my $dbq_drop_uniprot_go = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_go});
my $dbq_create_uniprot_go = $dbh_titanic->prepare(q{CREATE TABLE uniprot_go(
    uniprot_id INT UNSIGNED NOT NULL,
    go_id VARCHAR(10),
    term VARCHAR(256),
    evidence VARCHAR(11),
    project VARCHAR(64),
    KEY (uniprot_id))});

my $dbq_drop_ensembl = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS ensembl_ref});
my $dbq_create_ensembl = $dbh_titanic->prepare(q{CREATE TABLE ensembl_ref (
	ensembl_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
	reference_id INT UNSIGNED NOT NULL,
	reference_type ENUM('uniprot', 'cosmic'),
	ensg VARCHAR(16),
	enst VARCHAR(16),
	ensp VARCHAR(16),
	PRIMARY KEY (ensembl_id),
    KEY (reference_id))});

my $dbq_drop_sequence = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS sequence});
my $dbq_create_sequence = $dbh_titanic->prepare(q{CREATE TABLE sequence (
	sequence_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
	source ENUM('uniprot','cosmic', 'pdb'),
	source_id INT UNSIGNED NOT NULL,
	isoform INT UNSIGNED,
	sequence TEXT,
	sequence_md5 VARCHAR(32),
	PRIMARY KEY (sequence_id))});

my $dbq_drop_domain_map = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS domain_map});
my $dbq_create_domain_map = $dbh_titanic->prepare(q{CREATE TABLE domain_map (
    aggregate_id INT UNSIGNED NOT NULL,
	cosmic_id INT UNSIGNED NOT NULL,
	domain_type ENUM ('pfam_A', 'pfam_B', 'smart', 'prosite', 'other'),																 
	domain_acc VARCHAR(16),
	domain_name VARCHAR(16),
	domain_id INT,
	pfam_start INT,
	pfam_end INT,
	seq_start INT,
	seq_end INT,
	KEY (aggregate_id)
	)});

# Drop/create tables

$dbq_drop_uniprot->execute;
$dbq_create_uniprot->execute;

$dbq_drop_uniprot_alternative_name->execute;
$dbq_create_uniprot_alternative_name->execute;

$dbq_drop_uniprot_submitted_name->execute;
$dbq_create_uniprot_submitted_name->execute;

$dbq_drop_uniprot_gene_synonym->execute;
$dbq_create_uniprot_gene_synonym->execute;

$dbq_drop_uniprot_refseq->execute;
$dbq_create_uniprot_refseq->execute;

$dbq_drop_uniprot_pdb->execute;
$dbq_create_uniprot_pdb->execute;

$dbq_drop_uniprot_mim->execute;
$dbq_create_uniprot_mim->execute;

$dbq_drop_uniprot_pid->execute;
$dbq_create_uniprot_pid->execute;

$dbq_drop_uniprot_reactome->execute;
$dbq_create_uniprot_reactome->execute;

$dbq_drop_uniprot_drugbank->execute;
$dbq_create_uniprot_drugbank->execute;

$dbq_drop_uniprot_go->execute;
$dbq_create_uniprot_go->execute;

$dbq_drop_ensembl->execute;
$dbq_create_ensembl->execute;

$dbq_drop_sequence->execute;
$dbq_create_sequence->execute;

$dbq_drop_domain_map->execute;
$dbq_create_domain_map->execute;

print STDERR "\n";

$dbh_titanic->disconnect;
