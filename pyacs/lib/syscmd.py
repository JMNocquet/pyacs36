"""
Runs a system commands
"""

class Error(Exception):
    """Base class for exceptions in module syscmd."""
    pass

class ReturnedCodeError(Exception):
    """Exception raised for returned code not 0."""
    pass


def getstatusoutput( cmd , verbose=False ):

#     if verbose:
#         return ( run_command_verbose(cmd) )
# 
#     else:
        
    import subprocess


# did not understand whether cmd should be a string or a list
#    if isinstance(cmd,str):
#        cmd = cmd.split()
#
#    print("-- Running system command:|%s|" % " ".join( cmd ) )

    if verbose:
        print("-- Running system command:|%s|" % cmd )

    try:
        completed = subprocess.run( \
            cmd,                    \
            check=True,             \
            shell=True,             \
            stdout=subprocess.PIPE, \
        )
    except subprocess.CalledProcessError as err:
        print('!!! ERROR in subprocess.run:', err)

    else:
        
        if completed.returncode != 0:

            raise ReturnedCodeError( "!!! ERROR returned code is not 0 %d:" , completed.returncode )
            print('!!! ERROR:', err)
            print('!!! stdout: {!r}'.format(
                len(completed.stdout),
                completed.stdout.decode('utf-8')))
            
        else:
            if verbose:
                print('stdout: {!r}'.format(
                    len(completed.stdout),
                    completed.stdout.decode('utf-8')))
############################################################################

def run_command_verbose(command):
    
    import subprocess
    import shlex
    
    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline().decode()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc


