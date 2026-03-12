"""Subprocess helpers."""


def run_cmd(cmd, shell=True, check=True, capture_output=True, text=True):
    """Run a command using subprocess and check for errors.

    Parameters
    ----------
    cmd : str or list
        Command to run. If str, will be run with shell=True
    shell : bool, optional
        Whether to use shell. Default is True
    check : bool, optional
        Whether to raise CalledProcessError if return code is non-zero. Default is True
    capture_output : bool, optional
        Whether to capture stdout and stderr. Default is True
    text : bool, optional
        Whether to return output as string. Default is True

    Returns
    -------
    subprocess.CompletedProcess
        Object containing returncode, stdout, stderr

    Raises
    ------
    subprocess.CalledProcessError
        If check=True and command returns non-zero exit code
    """
    import subprocess

    try:
        result = subprocess.run(cmd, shell=shell, check=check, capture_output=capture_output, text=text)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Command: {e.cmd}")
        print(f"Output: {e.output}")
        print(f"Error: {e.stderr}")
        raise

