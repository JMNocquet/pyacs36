def info(self):
    """
    Basic informations about time series included in the Gts instance.
    """
    import pyacs.message.message as MESSAGE
    import pyacs.lib.astrotime as at
    import numpy as np

    # number of ts
    MESSAGE("Number of time series: %d" % self.n())
    # date
    min_decyear = 3000.
    max_decyear = -1
    l1000days = []
    l25days = []
    lduration = []
    llon = []
    llat = []
    for code in self.lcode():
        lduration.append(self.__dict__[code].data[-1, 0] - self.__dict__[code].data[0, 0])
        llon.append(self.__dict__[code].lon)
        llat.append(self.__dict__[code].lat)
        if self.__dict__[code].data[0, 0] < min_decyear:
            min_decyear = self.__dict__[code].data[0, 0]
        if self.__dict__[code].data[-1, 0] > max_decyear:
            max_decyear = self.__dict__[code].data[-1, 0]
        if self.__dict__[code].data.shape[0] > 1000: l1000days.append(code)
        if self.__dict__[code].data[-1, 0] - self.__dict__[code].data[0, 0] > 2.5: l25days.append(code)
    MESSAGE("Start/End Duration in decimal year: %10.4lf %10.4lf %8.4lf" % (
    min_decyear, max_decyear, max_decyear - min_decyear))
    ndays = int(at.decyear2mjd(max_decyear) - at.decyear2mjd(min_decyear))
    MESSAGE("Start/End Duration doy: %04d-%03d %04d-%03d %d days" % (
    int(min_decyear), at.decyear2dayno(min_decyear)[0], int(max_decyear), at.decyear2dayno(max_decyear)[0], ndays))
    str_sdate = at.decyear2datetime(min_decyear).isoformat(" ", "minutes")
    str_edate = at.decyear2datetime(max_decyear).isoformat(" ", "minutes")
    MESSAGE("Start/End : %s   %s " % (str_sdate, str_edate))
    MESSAGE("Number of sites with more than 1000 observations: %d" % len(l1000days))
    MESSAGE("Number of sites with more than 2.5 years of observations: %d" % len(l25days))
    MESSAGE("Mean time series duration %.1lf years. Median time series duration %.1lf years" % (
    np.mean(np.array(lduration)), np.median(np.array(lduration))))
    MESSAGE("Network mean location: (%.1lf,%.1lf). Network median location: (%.1lf,%.1lf)" % (
    np.mean(np.array(llon)), np.mean(np.array(llat)), np.median(np.array(llon)), np.median(np.array(llat))))
    return self
