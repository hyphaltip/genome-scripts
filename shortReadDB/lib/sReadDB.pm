=head1 NAME

sReadDB - short read DB interface

=head1 SYNOPSIS

=head1 USAGE

=head1 AUTHORS Jason Stajich

=cut

package sReadDB;
use strict;
use Env qw(HOME);
use DBI;
our ($InnoDB, $DBTYPE);

# connection magic

sub connect {
    my ($user,$pass,$host,$dbname,$dbtype) = @_;
    my $dsn;
    if( $dbtype =~ /mysql/i ) {
	($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;    
	$dsn = sprintf('dbi:mysql:database=%s;host=%s:mysql_local_infile=1:mysql_server_prepare=1',$dbname,$host);
	$DBTYPE = 'mysql';	
    } elsif( $dbtype =~ /Pg|postgres/i) {
	$dsn = sprintf("dbi:Pg:dbname=%s",$dbname);
	$DBTYPE = 'pg';
    }
    return DBI->connect($dsn,$user,$pass,
			{RaiseError => 1,
			 AutoCommit => 0});
}

# loading routines

# this routine loads the project info into the project table
sub get_project_id {
    my ($dbh,$project,$desc,$code) = @_;
    my $sth = $dbh->prepare_cached(qq{SELECT project_id,name,description,code FROM project WHERE name = ?});
    $sth->execute($project);
    if( $sth->rows ) {
	my ($id,$r_name,$r_desc,$r_code) = @{$sth->fetchrow_arrayref};
	$sth->finish;
	if( (defined $desc && $r_desc ne $desc) || 
	    (defined $code && $r_code ne $code) ) {
	    $sth = $dbh->prepare_cached(qq{UPDATE project 
					       SET description= ? AND 
					       code = ? where project_id = ?
					   });
	    $sth->execute($desc,$code,$id);
	    if( $@ ) {
		$sth->finish;
		warn("Cannot update project $r_name with $code, $desc:\n",
		     $dbh->errstr,"\n");
		$dbh->rollback;
		return -1;
	    }
	}	
	$sth->finish;
	return $id;
    }
    $sth = $dbh->prepare_cached(qq{INSERT INTO project (name,description,code) 
				       VALUES (?,?,?)
				   });    
    $sth->execute($project,$desc,$code);
    
    my $id;
    if( $DBTYPE eq 'mysql' ) {
	$id = $dbh->{mysql_insertid};
    } else {
	$id = $dbh->last_insert_id(undef,undef,"project",undef);
    }
    if( $@ || ! $id) {
	warn("error ",$dbh->errstr,"\n");
	$dbh->rollback;
	return -1;
    }
    $dbh->commit if $InnoDB || $DBTYPE eq 'pg';
    $sth->finish;
    return $id;
}

# this routine loads the chromosome info into the chromosome table
# from a sequence file
sub load_chromosomes_seq {
    my ($dbh,$org,$file,$format) = @_;
    my $rc = 0;
    my $isth = $dbh->prepare_cached(qq{INSERT INTO chromosome
					   (organism,name,length,description)
					   VALUES (?,?,?,?)});
    my $in = Bio::SeqIO->new(-file   => $file,
			     -format => $format);
    while( my $seq = $in->next_seq ) {
	$isth->execute($org,$seq->display_id,$seq->length,$seq->description);
	if( $@ ) {
	    warn($dbh->errstr,"\n");
	    $dbh->rollback;
	    $rc = -1;
	    last;
	}	
    }
    if( $rc ) {
	$dbh->commit if $InnoDB || $DBTYPE eq 'pg';
    }
    $isth->finish;
    $rc;
}

# get the bulk chromsome info as a hash reference based on organism name
sub get_chromosome_info {
    my ($dbh,$organism) = @_;
    my $d = {};
    my ($id,$name,$length);
    my $sth = $dbh->prepare_cached(qq{SELECT chromosome_id, name, length 
					  FROM chromosome WHERE organism = ?});
    $sth->execute($organism);
    $sth->bind_columns(\$id,\$name,\$length);
    while( $sth->fetch ) {
	$d->{$name} = [ $id,$length];
    }
    $sth->finish;
    return $d;
}

# get the method id for a prog/param combo
#   and will load and create if necessary
sub get_method_id {
    my ($dbh,$prog,$param) = @_;
    my $qsth = $dbh->prepare_cached(qq{SELECT method_id FROM analysis_method 
				       WHERE program = ? AND program_params = ?
				   });
    $qsth->execute($prog,$param);
    my $id;
    if( $qsth->rows ) {
	($id) = @{$qsth->fetchrow_arrayref};
	$qsth->finish;
	return $id if $id;
    }
    $qsth->finish;

    my $isth = $dbh->prepare_cached(qq{INSERT INTO analysis_method
					   (program,program_params)
					   VALUES (?,?)});    
	warn("sending in $prog $param\n");
    $isth->execute($prog,$param);    
    if( $DBTYPE eq 'mysql' ) {
	$id = $dbh->{mysql_insertid};
    } else {
	$id = $dbh->last_insert_id(undef,undef,"analysis_method",undef);
    }
    if( $@ || ! $id) {
	warn("error ",$dbh->errstr,"\n");
	$dbh->rollback;
	return -1;
    }
    $dbh->commit if $InnoDB || $DBTYPE eq 'pg';
    $isth->finish;
    return $id;
}

# load my local .my.cnf for user/pass combo
sub read_cnf {
    my ($user,$pass) = @_;
    if( -f "$HOME/.my.cnf") {
        open(IN,"$HOME/.my.cnf");
        while(<IN>) {
            if(/user(name)?\s*=\s*(\S+)/ ) {
                $user = $2;
            } elsif( /pass(word)\s*=\s*(\S+)/ ) {
                $pass = $2;
            }
        }
        close(IN);
    }
    return ($user,$pass);
}

1;
