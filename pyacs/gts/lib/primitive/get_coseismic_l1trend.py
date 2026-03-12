def get_coseismic_l1trend(self, date):
    """
    Get coseismic offsets using pyacs customized l1 trend filter

    param:
    date: date in decimal year, or datetime.datetime, or (mday,month,year)
    """

    # import
    from datetime import datetime
    import pyacs.lib.astrotime as at
    import numpy as np

    # decipher date and convert it into integer mjd
    if type(date) is tuple:
        mjd_eq = int(at.cal2mjd(date[0], date[1], date[2]))
    if type(date) is float and date > 1980:
        mjd_eq = int(at.decyear2mjd(date))

    if type(date) is datetime:
        mjd_eq = int(at.datetime2mjd(date))

    # extract one year of data surrounding the eq_date
    date_mjd = np.array(at.decyear2mjd(self.data[:, 0]), dtype=int)
    np_idx = np.where((date_mjd > mjd_eq - 183) & (date_mjd < mjd_eq + 183) & (date_mjd != mjd_eq))[0]
    np_mjd = date_mjd[np_idx]
    wts = self.copy()
    wts.data_xyz = None
    wts.data = self.data[np_idx, :]
    # apply l1 trend
    l1wts = wts.l1trend(min_samples_per_segment=2)
    # get coseismic offset from l1-trend filtered data
    idx_before_eq = np.where(np_mjd < mjd_eq)[0][-1]
    idx_after_eq = np.where(np_mjd > mjd_eq)[0][0]
    coseismic = l1wts.data[idx_after_eq, 1:4] - l1wts.data[idx_before_eq, 1:4]

    # get uncertainty
    # uses 10 days before and after
    idx_before_eq = np.where(np_mjd < mjd_eq)[0][-10:]
    idx_after_eq = np.where(np_mjd > mjd_eq)[0][:10]
    diff_before = np.std(wts.data[idx_before_eq, 1:4] - l1wts.data[idx_before_eq, 1:4], axis=0)
    diff_after = np.std(wts.data[idx_after_eq, 1:4] - l1wts.data[idx_after_eq, 1:4], axis=0)
    uncertainty = np.sqrt(diff_before ** 2 + diff_after ** 2)

    # print(coseismic*1.E3 , uncertainty)
    # wts.plot(superimposed=l1wts,center=False)

    return np.append(coseismic, uncertainty)
