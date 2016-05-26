#!/usr/local/gnu/bin/perl -w
####  _______              __
#### / ___/ /  ___  __ _  / /  ___
####/ /__/ _ \/ _ \/  ' \/ _ \/ _ \
####\___/_//_/\___/_/_/_/_.__/\___/ 
####
####
#### This software is copyright (C) by the Lawrence Berkeley
#### National Laboratory.  Permission is granted to reproduce
#### this software for non-commercial purposes provided that
#### this notice is left intact.
#### 
#### It is acknowledged that the U.S. Government has rights to
#### this software under Contract DE-AC03-765F00098 between
#### the U.S.  Department of Energy and the University of
#### California.
####
#### This software is provided as a professional and academic
#### contribution for joint exchange. Thus it is experimental,
#### is provided ``as is'', with no warranties of any kind
#### whatsoever, no support, no promise of updates, or printed
#### documentation. By using this software, you acknowledge
#### that the Lawrence Berkeley National Laboratory and
#### Regents of the University of California shall have no
#### liability with respect to the infringement of other
#### copyrights by any part of this software.
####
#################################################
###  This is the MultiDo processer.
### Interface is
### sub MultiDoProc::procMultiDoMacros(inputfile, outputfile,
###                                SpaceDim, debug)
###
###  reads in input file 
#################################################

package MakeAllCapsProc;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.23;
@ISA = qw(Exporter);
@EXPORT = qw(&procMakeAllCaps);
@EXPORT_OK = qw(&procMakeAllCaps);
@EXPORT_TAGS= ();

sub MakeAllCapsProc::makeAllCaps
{
    use strict;
    my ($FINFile, $FOUTFile, $SpaceDim, $debug) = @_;
    if($debug)
    {
        print "MakeAllCapsProc: \n";
        print "input file  = $FINFile \n";
        print "output file = $FOUTFile \n";
        print "SpaceDim    = $SpaceDim \n";
    }
    
    open(FOUT,">" . $FOUTFile) 
        or die "Error: cannot open output file " . $FOUTFile . "\n";
    open(FIN,"<" . $FINFile) 
        or die "Error: cannot open input file " . $FINFile . "\n";

    my $ibuf = "";
    my $obuf = "";
    ###put the entire input file into a string buffer
    while (defined( $ibuf = <FIN> )) 
    {
        ##only uppercase lines that do not start in # 
        if($ibuf =~ m/^\#/i)
        {
            $obuf .= $ibuf;
        }
        else
        {
            $obuf .= uc($ibuf);
        }
    }
    
    print FOUT $obuf;

    ###need to close files to make this modular.
    close(FIN);
    close(FOUT);
    return 1;   
}
###i have no idea why this is here.
###the perl cookbook book told me to put it there.
###really.
1;
