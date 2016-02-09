#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();
use LWP::UserAgent;

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Fetch Biomart - retrieve ensembl-Uniprot cross-references via Biomart web/xml interface\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

my $xml_query_fn = $mokcatanic_data_root . "biomart_query.xml";
my $biomart_fn = $mokcatanic_data_root . "biomart.tsv";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

# Go

# Read XML query

my $xml_query;
open my $QUERY, "< $xml_query_fn" || die "+++ Cannot open $xml_query_fn to read query: $?\n";
while (<$QUERY>) {
    $xml_query .= $_;
}
close $QUERY || die "+++ Cannot close $xml_query_fn: $?\n";

# Query BIOMART

my $path="http://www.biomart.org/biomart/martservice?";
my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml_query."\n");
my $ua = LWP::UserAgent->new;

my $response;

open my $RESULT, "> $biomart_fn" || die "+++ Cannot open $biomart_fn to write result data: $?\n";
$ua->request($request,
sub{
    my($data, $response) = @_;
    if ($response->is_success) {
        print $RESULT "$data";
    }
    else {
        warn ("Problems with the web server: ".$response->status_line);
    }
},1000);
close $RESULT || die "+++ Cannot close $biomart_fn: $?\n";

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
