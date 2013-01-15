#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'metagenome_assembly' ) || print "Bail out!\n";
}

diag( "Testing metagenome_assembly $metagenome_assembly::VERSION, Perl $], $^X" );
