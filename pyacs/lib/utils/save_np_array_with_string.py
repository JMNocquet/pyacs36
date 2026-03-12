"""Save numpy arrays with an additional string column."""


def save_np_array_with_string(A, S, fmt, outfile, comment=''):
    """Save a numpy array as a text file with a string column per row.

    Equivalent to np.savetxt but S is a 1D string array appended to each row.

    Parameters
    ----------
    A : numpy.ndarray
        2D array to save.
    S : array_like
        1D string array; len(S) must equal A.shape[0].
    fmt : str
        Format string for each row (e.g. "%03d %03d %s").
    outfile : str
        Output file path.
    comment : str, optional
        Comment line written at top. Default is ''.

    Examples
    --------
    >>> A = np.arange(4).reshape(-1, 2)
    >>> S = np.array(['lima', 'quito'])
    >>> save_np_array_with_string(A, S, "%03d %03d %s", 'test.dat')
    """
    import numpy as np

    if S.size != A.shape[0]:
        raise TypeError

    out = open(outfile, 'w+')

    if comment.strip() == '':
        comment = '# '

    if comment[0] != '#':
        comment = '# ' + comment

    out.write("%s\n" % comment)

    for i in np.arange(A.shape[0]):
        fmtstr = (fmt % tuple(list(A[i, :]) + [S[i]]))
        out.write("%s\n" % fmtstr)

    out.close()

