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

my $dbq_drop_uniprot_accession = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_accession});
my $dbq_create_uniprot_accession = $dbh_titanic->prepare(q{CREATE TABLE uniprot_accession (
    uniprot_id INT UNSIGNED NOT NULL,
    accession VARCHAR(10),
    KEY(uniprot_id))});

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
    codes_for ENUM('gene','protein'),
	sequence TEXT,
	sequence_md5 VARCHAR(32),
	PRIMARY KEY (sequence_id))});

my $dbq_drop_domain_map = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS domain_map});
my $dbq_create_domain_map = $dbh_titanic->prepare(q{CREATE TABLE domain_map (
    aggregate_id INT UNSIGNED NOT NULL,
    uniprot_domain_id INT UNSIGNED NOT NULL,
	KEY (aggregate_id),
    KEY (uniprot_domain_id)
	)});

my $dbq_drop_cvu = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS cosmic_vs_uniprot});
my $dbq_create_cvu = $dbh_titanic->prepare(q{CREATE TABLE cosmic_vs_uniprot (
    cosmic_id INT UNSIGNED NOT NULL,
    uniprot_id INT UNSIGNED,
    identity double,
    PRIMARY KEY (cosmic_id))});

my $dbq_drop_alignment_mapping = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS alignment_mapping});
my $dbq_create_alignment_mapping = $dbh_titanic->prepare(q{CREATE TABLE alignment_mapping(
    query_db enum('cosmic','uniprot','domain','pdb') default NULL,
    hit_db enum('cosmic','uniprot','domain','pdb') default NULL,
    query_id int(10) unsigned default NULL,
    hit_id int(10) unsigned default NULL,
    query_pos int(10) unsigned default NULL,
    hit_pos int(10) unsigned default NULL,
    query_res char(1) default NULL,
    hit_res char(1) default NULL,
    gap enum('none','query','hit') default NULL,
    KEY (query_db),
    KEY (query_id)
    )});

my $dbq_drop_counts = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS agg_sites});
my $dbq_create_counts = $dbh_titanic->prepare(q{CREATE TABLE agg_sites(
    aggregate_id INT UNSIGNED,
    Agl INT UNSIGNED,
    CNS INT UNSIGNED,
    Eye INT UNSIGNED,
    Mng INT UNSIGNED,
    Pty INT UNSIGNED,
    Adr INT UNSIGNED,
    Pth INT UNSIGNED,
    Sal INT UNSIGNED,
    Thy INT UNSIGNED,
    Bil INT UNSIGNED,
    GIT INT UNSIGNED,
    LIn INT UNSIGNED,
    Oes INT UNSIGNED,
    Pnc INT UNSIGNED,
    SIn INT UNSIGNED,
    Sto INT UNSIGNED,
    UAT INT UNSIGNED,
    Bon INT UNSIGNED,
    Ski INT UNSIGNED,
    SoT INT UNSIGNED,
    Bre INT UNSIGNED,
    Cvx INT UNSIGNED,
    End INT UNSIGNED,
    FTu INT UNSIGNED,
    GTr INT UNSIGNED,
    Ova INT UNSIGNED,
    Pen INT UNSIGNED,
    Pla INT UNSIGNED,
    Pro INT UNSIGNED,
    Tes INT UNSIGNED,
    Vag INT UNSIGNED,
    Vul INT UNSIGNED,
    HLT INT UNSIGNED,
    Thm INT UNSIGNED,
    Kid INT UNSIGNED,
    Lvr INT UNSIGNED,
    Lng INT UNSIGNED,
    Plr INT UNSIGNED,
    UTr INT UNSIGNED,
    oth INT UNSIGNED,
    NS INT UNSIGNED,
    PRIMARY KEY (aggregate_id))});

my $dbq_drop_aggmem = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS aggregate_membership});
my $dbq_create_aggmem = $dbh_titanic->prepare(q{CREATE TABLE aggregate_membership (
    complete_id INT UNSIGNED,
    aggregate_id INT UNSIGNED,
    KEY(aggregate_id))});

my $dbq_drop_uniprot_domain = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS uniprot_domain});
my $dbq_create_uniprot_domain = $dbh_titanic->prepare(q{CREATE TABLE uniprot_domain (
    uniprot_domain_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    domain_type ENUM('pfam_A','pfam_B','inter'),
    domain_acc VARCHAR(16),
    domain_id VARCHAR(16),
    uniprot_id INT UNSIGNED,
    uniprot_start INT UNSIGNED,
    uniprot_end INT UNSIGNED,
    map_start INT UNSIGNED,
    map_end INT UNSIGNED,
    mutated INT UNSIGNED,
    PRIMARY KEY (uniprot_domain_id),
    KEY (uniprot_id))});

my $dbq_drop_pdb = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS pdb});
my $dbq_create_pdb = $dbh_titanic->prepare(q{CREATE TABLE pdb (
    pdb_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    pdb_code VARCHAR(4),
    pdb_chain CHAR,
    PRIMARY KEY (pdb_id))});

# Drop/create tables

$dbq_drop_uniprot->execute;
$dbq_create_uniprot->execute;

$dbq_drop_uniprot_accession->execute;
$dbq_create_uniprot_accession->execute;

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

$dbq_drop_cvu->execute;
$dbq_create_cvu->execute;

$dbq_drop_alignment_mapping->execute;
$dbq_create_alignment_mapping->execute;

$dbq_drop_counts->execute;
$dbq_create_counts->execute;

$dbq_drop_aggmem->execute;
$dbq_create_aggmem->execute;

$dbq_drop_uniprot_domain->execute;
$dbq_create_uniprot_domain->execute;

$dbq_drop_pdb->execute;
$dbq_create_pdb->execute;

print STDERR "\n";

$dbh_titanic->disconnect;
