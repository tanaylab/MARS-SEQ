#!/usr/bin/env python

import sys, os, shutil
import subprocess
import time
import datetime

#
#   This script should be an analoug to the jobs runners script but not for jobs, but for scripts / commands.
#   So again basically there are 2 input files, one includes the commands that will be executed, and one is
#   the same csv as in the jobs_runner. The main motivation is to really replace the jobs_runners on the lab
#   cluster, which means that the command that we will run is the bsub command with the right arguments.
#   16-Mar-16
#   I am extending the script to work also by putting the command and / or the list in the running line.
#   I am also change the {NAME to {c for column.

def get_nof_jobs( queue = None ):

    cmd = [ 'bjobs' ]
    if ( queue ):
        cmd += [ '-q', queue ]

    jobs = subprocess.run( cmd, stdout = subprocess.PIPE )
    return( str( jobs.stdout ).count( '\\n' ) )


max_nof_jobs = -1
if ( len( sys.argv ) > 3 ):
    try:
        max_nof_jobs = int( sys.argv[ 3 ] )
    except:
        pass

queue_name = ''
if ( len( sys.argv ) > 4 and not sys.argv[ 4 ].startswith( "-" ) ):
    queue_name = sys.argv[ 4 ]

if ( max_nof_jobs > 0 and queue_name != '' ):
    keep_sending = 0
else:
    keep_sending = -1

debug = False
if ( "-d" in sys.argv ):
    debug = True

verbose = False
if ( "-v" in sys.argv ):
    verbose = True

arg1 = sys.argv[ 1 ]
if ( os.path.exists( arg1 ) ):
    cmd_file = open( arg1, 'r' )
    cmd_list = cmd_file.readlines()
    cmd_file.close()
    cmd_str = ''.join( cmd_list )
else:
    cmd_str = arg1

arg2 = sys.argv[ 2 ]
if ( arg2 == '-' ):
    parms_list = sys.stdin.readlines()
elif ( os.path.exists( arg2 ) and os.path.isfile( arg2 ) ):
    parms_file = open( arg2, 'r' )
    parms_list = parms_file.readlines()
    parms_file.close()
else:
    parms_list = []
    arr = arg2.split( ':' )
    if ( len( arr ) == 3 ):
        start = int( arr[ 0 ] )
        end = int( arr[ 1 ] )
        format_str = '\'{:0' + str( arr[ 2 ] ) + '}\''
        for n in range( start, end + 1 ):
            parms_list.append( format_str.format( n ) )
    elif ( len( arr ) == 2 ):
        start = int( arr[ 0 ] )
        end = int( arr[ 1 ] )
        for n in range( start, end + 1 ):
            parms_list.append( str( n ) )
    else:
        arr = arg2.split( ',' )
        for entry in arr:
            parms_list.append( entry.replace( '|', '\t' ) )



l = 1
for line in parms_list:

    if ( line.startswith( '#' ) ):
        continue

    arr = line.strip( '\n' ).split( '\t' )

    cmd = cmd_str
    n = 1
    for parm in arr:

        pattern = '{c' + str(n) + '}'
        cmd = cmd.replace( pattern, parm )
        n += 1

    cmd = cmd.replace( '{n}', str( l ) )
    cmd = cmd.replace( '{cc}', '_'.join( arr ) )
    if ( debug ):
        print( "#" + str( l ) )
        print( cmd )
    else:
        if ( max_nof_jobs < 0 or keep_sending < 0 ):
            if ( verbose ):
                print( "#" + str( l ) )
                print( cmd )
            os.system( cmd )
            keep_sending += 1
        else:
            nof_jobs = get_nof_jobs( queue_name )
            if ( nof_jobs >= max_nof_jobs ):
                while ( nof_jobs > max_nof_jobs ):
                    time.sleep( 10 )
                    nof_jobs = get_nof_jobs( queue_name )
            keep_sending = int( ( max_nof_jobs - nof_jobs ) / 2 )
            if ( verbose ):
                print( "#" + str( l ) )
                print( cmd )
            os.system( cmd )
            keep_sending += 1
    l += 1
