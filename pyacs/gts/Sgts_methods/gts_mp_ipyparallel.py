# def gts_mp(self, method, *args, **kwargs):
#     """
#     Sgts.gts_mp is the parallelized version of Sgts.gts, using ipyparallel

#     :param method: the Gts method to be applied to all time series (gts) of the current Sgts
#     :param *args, **kwargs: arguments of method to be used when calling method
#     :param ncpu: ncpu is an additional kwarg to define the number of CPU used. default is 4

#     :return: the new Sgts including the gts processed by method
#     """
#     from pyacs.gts.Sgts import Sgts
#     from tqdm import tqdm

#     import logging
#     import pyacs.message.message as MESSAGE
#     import pyacs.message.verbose_message as VERBOSE
#     import pyacs.message.error as ERROR
#     import pyacs.message.warning as WARNING
#     import pyacs.message.debug_message as DEBUG
#     import pyacs.debug

#     import ipyparallel as ipp
#     from icecream import ic

#     import inspect
#     VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

#     sgts = self
#     n = self.n()

#     # get func

#     import importlib
#     import inspect
#     import os

#     try:
#         func = getattr(self.__dict__[self.lcode()[0]], method)
#     except:
#         ERROR(("There is no Gts method called %s. Exit" % method), exit=True)

#     # get func path
#     ldir = []
#     OK = False
#     for cdir in inspect.getfile(func).split('/'):
#         if OK:
#             ldir.append(cdir)
#         if cdir == 'gts':
#             OK = True
#     ldir.pop()

#     str_4_import = 'pyacs.gts.' + '.'.join(ldir) + '.' + method
#     mod = importlib.import_module(str_4_import)
#     func = getattr(mod, method)

#     arg = [args, kwargs]

#     verbose = kwargs.get('verbose', False)
#     ncpu = kwargs.get('ncpu', 4)

#     if 'ncpu' in kwargs.keys():
#         kwargs.pop('ncpu')

#     MESSAGE("Will run %s method on %d Gts using %d CPU" % (method, sgts.n(), ncpu))

#     def worker_wrapper(gts, arg):
#         [args, kwargs] = arg
#         #return func(gts, *args, **kwargs)

#         try:
#             return func(gts, *args, **kwargs)
#         except:
#             gts.data = None
#             return gts

#     # request a cluster
#     with ipp.Cluster(n=ncpu, monitor_engines=True, log_level=40) as rc:
#     #with ipp.Cluster(n=ncpu) as rc:
#         # get a view on the cluster
#         view = rc.load_balanced_view()
#         # submit the tasks
#         asyncresult = view.map_async(worker_wrapper, sgts.lGts(), [arg] * n)
#         # wait interactively for results
#         asyncresult.wait_interactive()
#         # retrieve actual results
#         result = asyncresult.get()
#     # at this point, the cluster processes have been shutdown

#     # get results
#     new_ts = Sgts(read=False)
#     for gts in result:
#         if gts is not None:
#             new_ts.append(gts)

#     # print missing
#     lgood_code = new_ts.lcode(force_4_digit=True)
#     for code in sorted(sgts.lcode()):
#         if code not in lgood_code:
#             ERROR(("Problem while applying %s to %s. Skipped from output Sgts." % (method, code)))

#     return (new_ts)
def gts_mp_ipyparallel(self, method, *args, **kwargs):
    """Parallel version of Sgts.gts using ipyparallel with fault-tolerant recovery. Deprecated.

    Parameters
    ----------
    method : str
        Gts method to apply.
    *args : tuple
        Positional arguments for the method.
    **kwargs : dict
        Keyword arguments; ncpu (int) sets number of CPUs (default 4).

    Returns
    -------
    Sgts
        New Sgts with processed Gts.
    """
    from pyacs.gts.Sgts import Sgts
    import pyacs.message.message as MESSAGE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.verbose_message as VERBOSE
    from ipyparallel import Cluster
    from ipyparallel.error import EngineError
    import traceback
    import importlib
    import inspect

    VERBOSE("Running Sgts.gts_mp")

    sgts = self
    n = self.n()
    ncpu = kwargs.pop('ncpu', 4)
    verbose = kwargs.get('verbose', False)

    # Get the method
    try:
        func = getattr(self.__dict__[self.lcode()[0]], method)
    except:
        ERROR(f"No such method: {method}", exit=True)

    # Resolve correct import path for func
    ldir = []
    OK = False
    for cdir in inspect.getfile(func).split('/'):
        if OK:
            ldir.append(cdir)
        if cdir == 'gts':
            OK = True
    ldir.pop()

    str_4_import = 'pyacs.gts.' + '.'.join(ldir) + '.' + method
    mod = importlib.import_module(str_4_import)
    func = getattr(mod, method)

    arg = [args, kwargs]
    MESSAGE(f"Will run {method} on {n} Gts using {ncpu} CPU(s)")

    def worker_wrapper(gts, arg):
        [args, kwargs] = arg
        #return func(gts, *args, **kwargs)

        try:
            return func(gts, *args, **kwargs)
        except:
            gts.data = None
            return gts

    new_ts = Sgts(read=False)

    with Cluster(n=ncpu, monitor_engines=True) as rc:
        view = rc.load_balanced_view()
        asyncresult = view.map_async(worker_wrapper, sgts.lGts(), [arg] * n)

        print("⏳ Waiting for tasks...")
        asyncresult.wait(timeout=600)

        partial_results = []

        try:
            partial_results = asyncresult.get()
        except EngineError as ee:
            ERROR(f"💥 Engine crash: {ee}")
            # Attempt to retrieve what we can
            for i in range(len(sgts.lGts())):
                try:
                    partial_results.append(asyncresult[i])
                except Exception as e:
                    ERROR(f"⚠️ Task {i} retrieval failed: {type(e).__name__}: {e}")
                    partial_results.append(None)
        except Exception as e:
            ERROR(f"❌ Unexpected exception: {e}")
#            partial_results = [None] * len(sgts.lGts())
            # Attempt to retrieve what we can
            for i in range(len(sgts.lGts())):
                try:
                    partial_results.append(asyncresult[i])
                except Exception as e:
                    ERROR(f"⚠️ Task {i} retrieval failed: {type(e).__name__}: {e}")
                    partial_results.append(None)


        for res in partial_results:
            if isinstance(res, tuple) and len(res) == 3:
                code, gts_result, error = res
                if gts_result is None:
                    ERROR(f"⚠️ {code} failed:\n{error}")
                else:
                    new_ts.append(gts_result)
            elif res is not None:
                new_ts.append(res)
            else:
                ERROR(f"⚠️ Missing result (engine crash or task failure)")


    # Collect successful results
    for gts_result in partial_results:
        if gts_result is not None:
            new_ts.append(gts_result)
        else:
            ERROR(f"⚠️ Failed to process a Gts")

    # Check missing outputs

    lpb = []

    lgood_code = new_ts.lcode(force_4_digit=True)
    for code in sorted(sgts.lcode()):
        if code not in lgood_code:
            ERROR(f"❌ Missing output for {code} — skipped in final Sgts.")
            lpb.append(code)

    # Check for inconsistencies between input and output Sgts
    if new_ts.n() != self.n():
        ERROR(f"Number of Gts in output Sgts ({new_ts.n()}) does not match input Sgts ({self.n()}).")
    lcode_output = sorted([x[:4] for x in new_ts.lcode()])
    for code in self.lcode():
        if code not in lcode_output:
            WARNING(f"{code} is missing in output Sgts. Check why.")
            lpb.append(code)
        else:
            # Check for None values in l1trend_sts_aicc
            for code in self.lcode():
                if new_ts.__dict__ is None:
                    WARNING(f"{code} was not properly processed and is None.")
                    lpb.append(code)
                else:
                    if new_ts.__dict__[code].data is None:
                        WARNING(f"{code} has None data in l1trend_sts_aicc")
                        lpb.append(code)

    if len(lpb) > 0:
        WARNING(f"Number of missing output time series: {lpb}")
    else:
        MESSAGE("All time series were properly processed.")

    return new_ts.delnone()