###################################################################
## smooth a Gts time series
###################################################################

def smooth(self, window_len=11, window='hanning', in_place=False, verbose=False, component='NEU'):
    """
    smooth a time series
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print(error)
        return (self)

    ###########################################################################

    ###################################################################
    ## smoothing routines from http://wiki.scipy.org/Cookbook/SignalSmooth
    # changes numpy to np # JMN 18/07/2014
    ###################################################################

    def smooth_scipy(x, window_len=11, window='hanning'):
        """smooth the data using a window with requested size.

        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.

        input:
            x: the input signal
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.

        output:
            the smoothed signal

        example:

        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)

        see also:

        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter

        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")

        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than window size.")

        if window_len < 3:
            return x

        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

        s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
        # print(len(s))
        if window == 'flat':  # moving average
            w = np.ones(window_len, 'd')
        else:
            w = eval('np.' + window + '(window_len)')

        y = np.convolve(w / w.sum(), s, mode='valid')
        return y

    new_east = smooth_scipy(self.data[:, 1], window_len=window_len, window=window)
    new_north = smooth_scipy(self.data[:, 2], window_len=window_len, window=window)
    new_up = smooth_scipy(self.data[:, 3], window_len=window_len, window=window)

    new_Gts = self.copy()

    if in_place:
        return (self)
        del new_Gts
    else:
        new_Gts.data[:, 1] = new_east[window_len // 2 - 1:new_Gts.data[:, 1].shape[0] + window_len // 2 - 1]
        new_Gts.data[:, 2] = new_north[window_len // 2 - 1:new_Gts.data[:, 1].shape[0] + window_len // 2 - 1]
        new_Gts.data[:, 3] = new_up[window_len // 2 - 1:new_Gts.data[:, 1].shape[0] + window_len // 2 - 1]
        return (new_Gts)
