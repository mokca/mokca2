#!/usr/bin/perl

use FileHandle;

# Autoflush STDERR

STDERR->autoflush(1);

# Explain purpose

print STDERR "::: Fetch UNIPROT Human fasta lists\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

my $uniprot_fasta_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/HUMAN.fasta.gz";
my $uniprot_trembl_fasta_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz";
my $uniprot_sprot_fasta_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
my $uniprot_xml_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz";
my $uniprot_trembl_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.xml.gz";
my $uniprot_isoforms_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz";

system("wget", $uniprot_fasta_url, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++ Cannot fetch $uniprot_fasta_url: $?\n";
system("wget", $uniprot_trembl_fasta_url, "-N", "-P", $mokcatanic_local_root) == 0 || die "+++ Cannot fetch $$uniprot_trembl_fasta_url: $?\n";
system("wget", $uniprot_sprot_fasta_url, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++ Cannot fetch $$uniprot_sprot_fasta_url: $?\n";
system("wget", $uniprot_xml_url, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++ Cannot fetch $uniprot_xml_url: $?\n";
system("wget", $uniprot_trembl_url, "-N", "-P", $mokcatanic_local_root) == 0 || die "+++ Cannot fetch $uniprot_trembl_url: $?\n";
system("wget", $uniprot_isoforms_url, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++ Cannot fetch $uniprot_isoforms_url: $?\n";
