def find_offsets_ivel(self,
                      ivelts=None,
                      lcomponent='EN',
                      threshold_offset=100,
                      offsets_file=None):
    """
    Find offsets using l1trend filtering and instantaneous velocity

    parameters:
    -----------
    ivelts: Sgts
        ivel time series
    lcomponent: string
        component to use for the analysis
    threshold_offset: float
        threshold for the offset in mm/yr
    offsets_file: string
        name of the file to append the offsets. If None, no file is written

    notes:
    ------
    The function uses the l1trend of the time series. It then detects large ivel values and check whether they last for more than one day.
    If so, an offset is estimated by integrating ivel over two successive dates and offset is applied to the time series.
    Format of offsets_file is:
    code date (datetime.isoformat()) lon lat offset_e offset_n offset_u std_offset_e std_offset_n std_offset_u comment
    """
    # import
    import numpy as np
    import pyacs.lib.astrotime as at

    import threading

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug
    import pyacs.verbose
    from icecream import ic

    # offset file format
    ("%s %s %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %s")

    # new_ts
    new_ts = self.copy()

    # detrend
    try:
        dts = self.detrend_median()
    except:
        WARNING('Using detrend instead of detrend_median' % self.code)
        try:
            dts = self.detrend()
        except:
            ERROR('Error in detrend for site %s' % self.code)
            return None

    # compute l1trend
    pyacs.verbose('ERROR')
    try:
        l1ts = dts.l1trend(lcomponent=lcomponent)
    except:
        ERROR('Error in l1trend for site %s' % self.code)
        return None

    pyacs.verbose()

    # compute residual time series
    # rests =

    # compute ivel
    ivelts = l1ts.ivel()

    # compute max ivel

    norm_ivel = np.sqrt(ivelts.data[:, 1] ** 2 + ivelts.data[:, 2] ** 2) * 1.E3

    lidx = np.argwhere(norm_ivel > threshold_offset).flatten()

    print('dates to be checked ')
    for idx in lidx:
        print(at.decyear2datetime(ivelts.data[idx, 0]).isoformat())

    l3dates = []

    for idx in lidx:

        # check whether max_vel is only one date
        single_date = False
        if idx == ivelts.data.shape[0]:
            single_date = True
        if idx == 0:
            single_date = True

        component = np.argmax(np.fabs(ivelts.data[idx, 1:3])) + 1

        if np.fabs(ivelts.data[idx-1, component]) < 0.99 * np.fabs(ivelts.data[idx, component]) and \
                np.fabs(ivelts.data[idx+1, component]) < 0.99 * np.fabs(ivelts.data[idx, component]):
            single_date = True

        if single_date:
            # apply offset
            delta_date_year = ivelts.data[idx, 0] - ivelts.data[idx - 1, 0]
            offset_date = np.array([[ivelts.data[idx, 0] - 0.25 / 365.25, -ivelts.data[idx, 1] * delta_date_year,
                                     -ivelts.data[idx, 2] * delta_date_year, 0]])
            str_date = at.decyear2datetime(ivelts.data[idx, 0]).isoformat()[:10]
            VERBOSE('Correcting offset at %s %s %.1lf %.1lf mm' % (
            self.code, str_date, offset_date[0, 1] * 1.E3, offset_date[0, 2] * 1.E3))
            new_ts = new_ts.apply_offsets(offset_date)
            # write offset to file
            if offsets_file is not None:
                str_date = at.decyear2datetime(ivelts.data[idx, 0]).isoformat()
                lock = threading.Lock()

                with lock:
                    with open(offsets_file, 'a+') as f:
                        f.write("%s %s %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %s\n" %
                                (self.code, str_date, self.lon, self.lat, offset_date[0, 1] * 1.E3,
                                 offset_date[0, 2] * 1.E3, 0., 0., 0., 0., 'automated from PYACS'))

        else:
            # check duration of high velocity
            # forward in time
            VERBOSE('Max ivel for site %s is several days indicating a transient at %s' % (self.code,at.decyear2datetime(ivelts.data[idx,0]).isoformat()[:10]) )
            idx_ref = idx
            ldix_large_ivel = [idx_ref]
            print(np.fabs(ivelts.data[idx+1, component])*1E3 , threshold_offset / 4)
            while np.fabs(ivelts.data[idx+1, component])*1E3 > threshold_offset / 4:
                print(idx)
                print(ldix_large_ivel)
                idx += 1
                if idx == ivelts.data.shape[0]:
                    break
                ldix_large_ivel.append(idx)

            idx = idx_ref

            print(np.fabs(ivelts.data[idx-1, component])*1E3 , threshold_offset / 4)
            while np.fabs(ivelts.data[idx-1, component])*1E3 > threshold_offset / 4 :
                print(idx)
                print(ldix_large_ivel)
                idx -= 1
                if idx == 0:
                    break
                ldix_large_ivel.append(idx)
            VERBOSE('Duration of high ivel is %d days' % len(ldix_large_ivel))
            print(ldix_large_ivel)

            if len(ldix_large_ivel) == 3:
                l3dates.append(sorted(ldix_large_ivel))

    # handle the case of 3 consecutive days of large velocity
    # makes unique
    ul3dates = [list(x) for x in set(tuple(x) for x in l3dates)]
    for period in ul3dates:
        # check if an offset is hidden in 3 days high velocity
        print(period)
        print([x.isoformat() for x in at.decyear2datetime(ivelts.data[period, 0])])

        idx_ref = period[0]
        idx = idx_ref - 1
        print(dts.data[idx:idx + 6, 2] * 1.E3)
        print(l1ts.data[idx:idx + 6, 2] * 1.E3)
        new_segmentation = np.zeros(6)
        new_segmentation[0] = l1ts.data[idx, 2] * 1.E3
        new_segmentation[1] = l1ts.data[idx + 1, 2] * 1.E3
        ivb = (new_segmentation[1] - new_segmentation[0]) / (l1ts.data[idx + 1, 0] - l1ts.data[idx, 0])
        new_segmentation[2] = new_segmentation[1] + ivb * (l1ts.data[idx + 2, 0] - l1ts.data[idx + 1, 0])
        new_segmentation[4] = l1ts.data[idx + 4, 2] * 1.E3
        new_segmentation[5] = l1ts.data[idx + 5, 2] * 1.E3
        iva = (new_segmentation[5] - new_segmentation[4]) / (l1ts.data[idx + 5, 0] - l1ts.data[idx + 4, 0])
        new_segmentation[3] = new_segmentation[4] - iva * (l1ts.data[idx + 5, 0] - l1ts.data[idx + 4, 0])
        print(new_segmentation)

        chi2_orig = np.sum( (dts.data[idx:idx + 6, 2] * 1.E3 - l1ts.data[idx:idx + 6, 2] * 1.E3) ** 2)
        print('chi2 original ', chi2_orig)
        chi2_new = np.sum((dts.data[idx:idx + 6, 2] * 1.E3 - new_segmentation) ** 2)
        print('chi2 new      ', chi2_new)



    return new_ts
