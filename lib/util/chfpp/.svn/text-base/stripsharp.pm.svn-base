#!/usr/local/gnu/bin/perl -w
# stdin.pl 
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
###  This is the DTerm processer.
### Interface is
### sub DTermProc::procDTermMacros(inputfile, outputfile,
###                                SpaceDim, debug)
###
###  reads in input file 
#################################################

package StripSharpProc;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.23;
@ISA = qw(Exporter);
@EXPORT = qw(&stripSharp);
@EXPORT_OK = qw(&stripSharp);
@EXPORT_TAGS= ();

sub StripSharpProc::StripSharp
{

    use strict;
    my ($debug) = @_;
    if($debug)
    {
        print "StripSharpProc: \n";
    }
    
    while (defined(my $ibuf = <STDIN> )) 
    {

###     skip lines that start with # 
###     (or any number of spaces before it) -JNJ
        if($ibuf =~ m/^. \#/i)
        {
            next;
        }
        else
        {
            print $ibuf;
        }
    }
    
   return 1;   
}
###i have no idea why this is here.
###the perl cookbook book told me to put it there.
###really.
1;
