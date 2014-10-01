package database_setup;
use warnings;
use parent 'Exporter';
# Package to setup database names, usernames, hosts

our $mokca_db = "mokcatanic_rebuild";
our $mokca_host = "umma";
our $mokca_user = "mokca";
our $mokca_pw = "iwt4lp";

our @EXPORT = qw($mokca_db $mokca_host $mokca_user $mokca_pw);

1;