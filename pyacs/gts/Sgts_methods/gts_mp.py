def _worker_wrapper(gts, args, kwargs, method_path, method_name):
    """Worker for parallel processing of Gts methods.

    Parameters
    ----------
    gts : Gts
        Gts instance to process.
    args : tuple
        Positional arguments for the method.
    kwargs : dict
        Keyword arguments for the method.
    method_path : str
        Module path (e.g. 'pyacs.gts.lib.primitive.copy').
    method_name : str
        Method name to call.

    Returns
    -------
    tuple
        (code, result, error); error is None on success.
    """
    import traceback
    import importlib
    import pyacs.message.error as ERROR
    
    try:
        # Validate Gts object
        if not hasattr(gts, 'code'):
            return (None, None, "Invalid Gts object: missing 'code' attribute")
            
        # Import and get method
        try:
            # Import the module
            mod = importlib.import_module(method_path)
            
            # Get the method directly from the module
            func = getattr(mod, method_name)
            
            if not callable(func):
                return (gts.code, None, f"Method {method_name} is not callable in {method_path}")
                
        except ImportError as e:
            return (gts.code, None, f"Failed to import module {method_path}: {str(e)}")
        except AttributeError as e:
            return (gts.code, None, f"Method {method_name} not found in {method_path}: {str(e)}")
            
        # Execute method
        result = func(gts, *args, **kwargs)
        return (gts.code, result, None)
        
    except Exception as e:
        return (gts.code if hasattr(gts, 'code') else None, 
                None, 
                f"{str(e)}\n{traceback.format_exc()}")

def gts_mp(self, method, *args, **kwargs):
    """Apply a Gts method to each series in parallel (ProcessPoolExecutor).

    Parameters
    ----------
    method : str
        Gts method name (e.g. 'detrend').
    *args : tuple
        Positional arguments for the method.
    **kwargs : dict
        Keyword arguments; ncpu (int) sets number of CPUs (default 4).

    Returns
    -------
    Sgts
        New Sgts with processed Gts instances.
    """
    from pyacs.gts.Sgts import Sgts
    import inspect
    import importlib
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import pyacs.message.message as MESSAGE
    import pyacs.message.error as ERROR
    import pyacs.message.verbose_message as VERBOSE
    import multiprocessing
    import time
    from tqdm import tqdm

    start_time = time.time()
    VERBOSE("Running Sgts.gts_mp2")

    sgts = self
    # Get available CPUs, but leave one free for system
    available_cpus = multiprocessing.cpu_count()
    ncpu = min(kwargs.pop('ncpu', 4), max(1, available_cpus - 1))
    verbose = kwargs.get('verbose', False)

    if sgts.n() == 0:
        ERROR("No Gts objects to process", exit=True)

    # Resolve method location for re-importing in subprocess
    try:
        if not self.lcode():
            ERROR("No Gts objects available in the instance", exit=True)
        example_gts = self.__dict__[self.lcode()[0]]
        if not hasattr(example_gts, method):
            ERROR(f"Method {method} not found in Gts instance", exit=True)
        func = getattr(example_gts, method)
    except Exception as e:
        ERROR(f"Failed to resolve method {method}: {str(e)}", exit=True)

    # Get the module path
    try:
        func_file = inspect.getfile(func)
        if not func_file:
            ERROR(f"Could not determine file location for method {method}", exit=True)
            
        # Get the module name from the function
        module = inspect.getmodule(func)
        if module is None:
            ERROR(f"Could not determine module for method {method}", exit=True)
            
        method_path = module.__name__
        method_name = method
        
        VERBOSE(f"Method path: {method_path}, Method name: {method_name}")
        
    except Exception as e:
        ERROR(f"Failed to resolve method path: {str(e)}", exit=True)

    MESSAGE(f"Will run {method} on {sgts.n()} Gts using {ncpu} CPU(s)")

    # Prepare job arguments
    arglist = [(gts, args, kwargs, method_path, method_name) for gts in sgts.lGts()]
    new_ts = Sgts(read=False)
    processed = 0
    total = len(arglist)
    failed = 0

    # Parallel execution
    try:
        with ProcessPoolExecutor(max_workers=ncpu) as executor:
            futures = [executor.submit(_worker_wrapper, *job_args) for job_args in arglist]
            
            # Create progress bar if verbose
            pbar = tqdm(total=total, desc=f"Processing {method}", 
                       unit="site",
                       bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')

            for fut in as_completed(futures):
                processed += 1
                code, gts_result, error = fut.result()
                if gts_result is not None:
                    new_ts.append(gts_result)
                    if verbose:
                        elapsed = time.time() - start_time
                        rate = processed / elapsed if elapsed > 0 else 0
                        pbar.set_postfix({'rate': f'{rate:.1f} sites/sec'})
                else:
                    failed += 1
                    ERROR(f"⚠️ {code or 'Unknown'} failed:\n{error}")
                pbar.update(1)

            pbar.close()

    except Exception as e:
        ERROR(f"Parallel processing failed: {str(e)}", exit=True)

    # Check for missing outputs
    lgood_code = new_ts.lcode(force_4_digit=True)
    missing = [code for code in sorted(sgts.lcode()) if code not in lgood_code]
    if missing:
        ERROR(f"❌ Missing outputs for {len(missing)} sites: {', '.join(missing)}")

    elapsed = time.time() - start_time
    MESSAGE(f"Successfully processed {new_ts.n()} out of {sgts.n()} Gts in {elapsed:.1f} seconds "
            f"({new_ts.n()/elapsed:.1f} sites/sec)")
    if failed > 0:
        MESSAGE(f"⚠️ {failed} sites failed processing")
        
    return new_ts