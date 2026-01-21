#!/bin/perl

my @rea = ();
my $rc = 0;

# TWOSCALING SINGLE-ARG
$twoscaling_s = qr{ PROD\(
    ( #body
        (?:
            [^(),]+? # no parens or commas (non-greedy)
            |
            (\(  # inner-paren-group
                (?:
                    (?> [^()]+ ) # Non-parens without backtracking
                    |
                    (?-1) # Recurse to inner-paren-group
                )*
            \)) #end of inner-paren-group
        )*?
    )
    \s*\bTWOSCALING\s*\)  # terminal 'G)'
}x;

# TWOSCALING MULTI-ARG
$twoscaling_m = qr{ PROD\(
    ( #body
        (?:
            [^()]+? # non-greedy non-parens
            |
            (\(  # inner-paren-group
                (?:
                    (?> [^()]+ ) # Non-parens without backtracking
                    |
                    (?-1) # Recurse to inner-paren-group
                )*
            \)) #end of inner-paren-group
        )*?
    )
    \s*\bTWOSCALING\s*\)  # terminal 'G)'
}x;

# TWOSCALING SINGLE-ARG
$scaling_s = qr{ PROD\(
    ( #body
        (?:
            [^(),]+? # no parens or commas (non-greedy)
            |
            (\(  # inner-paren-group
                (?:
                    (?> [^()]+ ) # Non-parens without backtracking
                    |
                    (?-1) # Recurse to inner-paren-group
                )*
            \)) #end of inner-paren-group
        )*?
    )
    \s*\bSCALING\s*\)  # terminal 'G)'
}x;

# TWOSCALING MULTI-ARG
$scaling_m = qr{ PROD\(
    ( #body
        (?:
            [^()]+? # non-greedy non-parens
            |
            (\(  # inner-paren-group
                (?:
                    (?> [^()]+ ) # Non-parens without backtracking
                    |
                    (?-1) # Recurse to inner-paren-group
                )*
            \)) #end of inner-paren-group
        )*?
    )
    \s*\bSCALING\s*\)  # terminal 'G)'
}x;


my $str='PROD(X,PROD(A,B(X,Y),C(D(E,F)) Q(98798 - * q234890kljasdlk)),C) PROD(A TWOSCALING) SUM(PROD(A((B(C)))) SUM(PROD((X,A(X,(B(C))),2),2 TWOSCALING)';
#my $str='PROD(X G)';
#my @matches = $str =~ /PROD\(([^()]++)\s*\bG\s*\)/xg;
#my @matches = $str =~ /PROD\(([^()]+?)( G)\)/xg;
#my @matches = $str =~ s/$re/PFSCALE($1)/g;
#my $ns = $str =~ s/$re/PFSCALE($1)/xg and print $ns;
#foreach my $m (@matches) {
#  print "Match: $m\n";
#}
$str =~ s/$twoscaling_s/PFSCALE2($1)/xg and print "After 1: $str\n";
$str =~ s/$twoscaling_m/PFSCALE2(PROD($1))/xg and print "After 2: $str\n";