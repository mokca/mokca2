#!/usr/bin/perl -w

use Getopt::Long;
use FileHandle;
use DBI();
use XML::Twig;

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Import Uniprot XML dump - accessions, descriptions and cross-references for Human genes\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Local uniprot files

my $split_base_fn = $mokcatanic_local_root . "/uniprot_split_";

# Check we got chunk number from command line

if ( $xml_chunk eq "unset" ) {
    die "\n+++ Need chunk number\n";
}

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_create_entry = $dbh_titanic->prepare(q{INSERT INTO uniprot VALUES(NULL, ?, ?, ?, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)});
my $dbq_fetch_id = $dbh_titanic->prepare(q{SELECT uniprot_id FROM uniprot WHERE accession = ?});
my $dbq_add_altname = $dbh_titanic->prepare(q{INSERT INTO uniprot_alternative_name VALUES(?,?)});
my $dbq_add_subname = $dbh_titanic->prepare(q{INSERT INTO uniprot_submitted_name VALUES(?,?)});
my $dbq_add_gene_synonym = $dbh_titanic->prepare(q{INSERT INTO uniprot_gene_synonym VALUES(?,?)});
my $dbq_ensembl_insert = $dbh_titanic->prepare(q{INSERT INTO ensembl_ref VALUES(NULL, ?, 'uniprot', ?, ?, ?)});
my $dbq_PDB_insert = $dbh_titanic->prepare(q{INSERT INTO uniprot_pdb VALUES(?, ?, ?, ?, ?)});
my $dbq_MIM_insert = $dbh_titanic->prepare(q{INSERT INTO uniprot_mim VALUES(?, ?, ?)});
my $dbq_pid_insert = $dbh_titanic->prepare(q{INSERT INTO uniprot_pid VALUES(?,?,?)});
my $dbq_reactome_insert = $dbh_titanic->prepare(q{INSERT INTO uniprot_reactome VALUES(?,?,?)});
my $dbq_go_insert = $dbh_titanic->prepare(q{INSERT INTO uniprot_go VALUES(?,?,?,?,?)});
my $dbq_drugbank_insert = $dbh_titanic->prepare(q{INSERT INTO uniprot_drugbank VALUES(?,?,?)});
my $dbq_refseq_insert = $dbh_titanic->prepare(q{INSERT INTO uniprot_refseq VALUES(?,?,?)});
my $dbq_uniprot_update = $dbh_titanic->prepare(q{UPDATE uniprot SET recommended_name = ?, alternative_name_count = ?, submitted_name_count = ?, gene_name = ?, gene_synonym_count = ?, refseq_count = ?, pdb_count = ?,mim_count = ?, pid_count=?, reactome_count=?, drugbank_count=?, ensembl_count=?, go_count = ?, hgnc_numeric = ?, hgnc_name = ?, isoform_count = ?,
    kegg = ?, geneid = ?, chembl = ? WHERE uniprot_id = ?});

# Open chunk

my $chunk_fn = $split_base_fn . $xml_chunk . ".xml";
open my $uniprot_fh, '<', $chunk_fn || die "\n+++ Cannot open $chunk_fn to read data: $?\n";

# Parse XML from chunk

my $name_lmax = 0;

my $twig = XML::Twig->new( twig_handlers => { entry => \&uniprot_entry } );
$twig->parse($uniprot_fh);
$twig->purge;

#print "LONGESTESTEST STRING: $name_lmax\n";

# Close chunk

close $uniprot_fh || die "\n+++ Cannot close $chunk_fn\n";

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;

# XML Twig handler for UNIPROT entry entities

sub uniprot_entry {
	my ( $l_twig, $entry) = @_;
    
	my $isoform_c = 0;
	my ( $hgnc_id, $hgnc_numeric, $hgnc_gene );
    my ( $kegg_id, $gene_id_id, $chembl_id );
    $kegg_id = "";
    $gene_id_id = "";
    $chembl_id = "";
	my ( $sequence_l, $sequence_aa );
    my ( $refseq_c, $pdb_c, $mim_c, $pid_c, $reactome_c, $drugbank_c, $ensembl_c, $go_c );
    $refseq_c = 0;
    $pdb_c = 0;
    $mim_c = 0;
    $pid_c = 0;
    $reactome_c = 0;
    $drugbank_c = 0;
    $ensembl_c = 0;
    $go_c = 0;
    
    # Get basic details: accession, name, recommendedName (if any)
    
    my $dataset = $entry->{'att'}->{'dataset'};
    my $accession = $entry->first_child('accession')->text;
    my $name = $entry->first_child('name')->text;
    my $protein = $entry->first_child('protein');
    my $protein_recname = $protein->first_child('recommendedName');
    my $recommended_name;
    if ($protein_recname) {
        $recommended_name = $protein_recname->first_child('fullName')->text;
    }
    
    #print "$accession $name\n";
    #print "$recommended_name\n";
    
    # Create entry, then fetch ID
    
    $dbq_create_entry->execute($dataset, $accession, $name);
    $dbq_fetch_id->execute($accession);
    ( my $uniprot_id ) = $dbq_fetch_id->fetchrow_array;
    $dbq_fetch_id->finish;
    
    # Get alternative names, submitted names, add to tables linked to uniprot_id
    # Make recommended name first submitted name if missing
    
    my @protein_subnames = $protein->children('submittedName');
    my $subname_c = @protein_subnames;
    foreach my $protein_subname (@protein_subnames) {
        $submitted_name = $protein_subname->first_child('fullName')->text;
        if (!$recommended_name) {
            $recommended_name = $submitted_name;
            #print ">>>$recommended_name\n";
        }
        $dbq_add_subname->execute($uniprot_id, $submitted_name);
        #print "[$submitted_name] ";
    }
    
    my @protein_altnames = $protein->children('alternativeName');
    my $altname_c = @protein_altnames;
    foreach my $protein_altname (@protein_altnames) {
        $alternative_name = $protein_altname->first_child('fullName')->text;
        #print "($alternative_name) ";
        $dbq_add_altname->execute($uniprot_id, $alternative_name);
    }
    #print "\n";
    
    # Get gene names, first one in list is used as gene name, rest as
    # synonyms (even if they're ORF names)
    
    my $gene = $entry->first_child('gene');
    my $gene_name;
    my $synonym_c = 0;
    if ($gene) {
        my @entry_gene_names = $gene->children('name');
        my $first_gene_name = shift(@entry_gene_names);
        $gene_name = $first_gene_name->text;
        foreach my $entry_gene_name (@entry_gene_names) {
            my $gene_synonym = $entry_gene_name->text;
            $dbq_add_gene_synonym->execute($uniprot_id, $gene_synonym);
            $synonym_c++;
        }
    }
    
    # Fetch organism, currently only used as a sanity test
    
	my $organism = $entry->first_child('organism');
	my @names= $organism->children('name');
	my $common_name;
	foreach my $name (@names) {
		my $name_type = $name->{'att'}->{'type'};
		my $name_text = $name->text;
        
		if ($name_type eq 'common') {
			$common_name = $name_text;
		}
	}
    
	if ($common_name ne "Human") {
        print "\n+++ Organism not human ($common_name)\n";
    }
    
    # Count isoforms
    
    my @comments = $entry->children('comment');
    foreach my $comment (@comments) {
        my $comment_type = $comment->{'att'}->{'type'};
        if ($comment_type eq "alternative products") {
            #my @events = $comment->children('event');
            #foreach my $event (@events) {
            #	my $event_type = $event->{'att'}->{'type'};
            #	if ($event_type ne "alternative splicing") {
            #		print "$accession : $comment_type : $event_type\n";
            #	}
            #}
            
            my @isoforms = $comment->children('isoform');
            foreach $isoform (@isoforms) {
                
                # We don't currently do anything with the data describing isoforms
                
                my $isoform_name_text;
                
                $isoform_c++;
                my $isoform_id = $isoform->first_child('id')->text;
                my @isoform_names = $isoform->children('name');
                foreach $isoform_name(@isoform_names) {
                    $isoform_name_text = $isoform_name->text;
                }
                my $sequence_tag = $isoform->first_child('sequence');
                my $sequence_type = $sequence_tag->{'att'}->{'type'};
                my $sequence_ref;
                if ($sequence_type eq "described") {
                    $sequence_ref = $sequence_tag->{'att'}->{'ref'};
                }
            }
            
        }
    }
    if ($isoform_c == 0) { # If we have one isoform, the alternative products section will not be present
        $isoform_c = 1;
    }
    
    # Process database references
    
    my $hgnc_c = 0;
    $ensembl_c = 0;
    my @dbRefs = $entry->children('dbReference');
    foreach $dbRef (@dbRefs) {
        my $db_ref_type = $dbRef->{'att'}->{'type'};
        if ($db_ref_type eq "HGNC") {                           # HGNC
            
            $hgnc_c++;
            
            $hgnc_id = $dbRef->{'att'}->{'id'};
            if ($hgnc_id =~ /HGNC:(\d+)/) {
                $hgnc_numeric = $1;
            } else {
                print "$accession: Unusual HGNC ID: $hgnc_id\n";
            }
            
            my @properties = $dbRef->children('property');
            my $property;
            if ((scalar @properties) > 1) {
                print "$accession: Unusual number of properties: " . scalar(@properties) . "\n";
            }
            foreach $property (@properties) {
                if ($property->{'att'}->{'type'} eq 'gene designation') {
                    $hgnc_gene = $property->{'att'}->{'value'};
                }
            }
            
        } elsif ($db_ref_type eq "Ensembl") {                   # Ensembl
            
            my ( $ensg, $enst, $ensp );
            $enst = $dbRef->{'att'}->{'id'};
            my @ensembl_props = $dbRef->children('property');
            my $ensembl_prop;
            foreach $ensembl_prop (@ensembl_props) {
                if ($ensembl_prop->{'att'}->{'type'} eq "protein sequence ID") {
                    $ensp = $ensembl_prop->{'att'}->{'value'};
                } elsif ($ensembl_prop->{'att'}->{'type'} eq "gene ID") {
                    $ensg = $ensembl_prop->{'att'}->{'value'};
                }
            }
            
            $dbq_ensembl_insert->execute($uniprot_id, $ensg, $enst, $ensp);
            $ensembl_c++;
        } elsif ($db_ref_type eq "RefSeq") {                    # RefSeq
            
            my $refseq_id = $dbRef->{'att'}->{'id'};
            my $nucleotide_id = $dbRef->first_child('property')->{'att'}->{'value'};
            $dbq_refseq_insert->execute($uniprot_id, $refseq_id, $nucleotide_id);
            $refseq_c++;
            
        } elsif ($db_ref_type eq "PDB") {                       # PDB
            
            my $pdb_id = $dbRef->{'att'}->{'id'};
            my @PDB_props = $dbRef->children('property');
            my ( $pdb_method, $pdb_resolution, $pdb_chains );
            foreach my $PDB_prop (@PDB_props) {
                my $PDB_prop_type = $PDB_prop->{'att'}->{'type'};
                if ($PDB_prop_type eq "method") {
                    $pdb_method = $PDB_prop->{'att'}->{'value'};
                } elsif ($PDB_prop_type eq "resolution") {
                    $pdb_resolution = $PDB_prop->{'att'}->{'value'};
                } elsif ($PDB_prop_type eq "chains") {
                    $pdb_chains = $PDB_prop->{'att'}->{'value'};
                }
            }
            $dbq_PDB_insert->execute($uniprot_id, $pdb_id, $pdb_method, $pdb_resolution, $pdb_chains);
            $pdb_c++;
            
        } elsif ($db_ref_type eq "MIM") {                       # MIM
            
            my $mim_id = $dbRef->{'att'}->{'id'};
            my $mim_type;
            if ($dbRef->first_child('property')) {
                $mim_type = $dbRef->first_child('property')->{'att'}->{'value'};
            } else {
                $mim_type = "";
            }
            $dbq_MIM_insert->execute($uniprot_id, $mim_id, $mim_type);
            $mim_c++;
            
        } elsif ($db_ref_type eq "Pathway_Interaction_DB") {    # Pathway_IDB
            
            my $pathway_id = $dbRef->{'att'}->{'id'};
            my $pathway_name = $dbRef->first_child('property')->{'att'}->{'value'};
            $dbq_pid_insert->execute($uniprot_id, $pathway_id, $pathway_name);
            $pid_c++;

        } elsif ($db_ref_type eq "Reactome") {                  # Reactome
            
            my $reactome_id = $dbRef->{'att'}->{'id'};
            my $reactome_name = $dbRef->first_child('property')->{'att'}->{'value'};
            $dbq_reactome_insert->execute($uniprot_id, $reactome_id, $reactome_name);
            $reactome_c++;
            
        } elsif ($db_ref_type eq "GO") {                        # GO
            
            my ( $go_term, $go_evidence, $go_project );
            my $go_id = $dbRef->{'att'}->{'id'};
            my @go_props = $dbRef->children('property');
            foreach my $go_prop (@go_props) {
                my $go_prop_type = $go_prop->{'att'}->{'type'};
                if ($go_prop_type eq 'term') {
                    $go_term = $go_prop->{'att'}->{'value'};
                } elsif ($go_prop_type eq 'evidence') {
                    $go_evidence = $go_prop->{'att'}->{'value'};
                } elsif ($go_prop_type eq 'project') {
                    $go_project = $go_prop->{'att'}->{'value'};
                }
                $dbq_go_insert->execute($uniprot_id, $go_id, $go_term, $go_evidence, $go_project);
            }
            $go_c++
            
        } elsif ($db_ref_type eq "DrugBank") {                  # DrugBank
            
            my $drugbank_id = $dbRef->{'att'}->{'id'};
            my $drugbank_name = $dbRef->first_child('property')->{'att'}->{'value'};
            $dbq_drugbank_insert->execute($uniprot_id, $drugbank_id, $drugbank_name);
            $drugbank_c++;
            
        } elsif ($db_ref_type eq "KEGG") {                      # KEGG
            
            $kegg_id = $dbRef->{'att'}->{'id'};
            
        } elsif ($db_ref_type eq "GeneId") {                    # GeneID
            
            $gene_id_id = $dbRef->{'att'}->{'id'};
            
        } elsif ($db_ref_type eq "ChEMBL") {                    # ChEMBL
            
            $chembl_id = $dbRef->{'att'}->{'id'};
            
        }
    }
    
    $dbq_uniprot_update->execute($recommended_name, $altname_c, $subname_c, $gene_name, $synonym_c, $refseq_c, $pdb_c, $mim_c, $pid_c, $reactome_c, $drugbank_c, $ensembl_c, $go_c, $hgnc_numeric, $hgnc_gene, $isoform_c, $kegg_id, $gene_id_id, $chembl_id, $uniprot_id);
    
    # Pointless ticker
	
	if ($tick_c % $tick_step == 0) {
		print STDERR $tick_char;
		if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
			print STDERR "\n";
		}
	}
	$tick_c++;
    
    $entry->purge;
}
