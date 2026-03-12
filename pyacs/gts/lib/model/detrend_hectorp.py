from __future__ import annotations

def detrend_hectorp(self, component='NEU', **kwargs):
    """Detrend a time series using the Hector model (HECTOR from Machiel Bos).

    References: Bos et al. (2013), J. Geodesy 87(4), 351-360.
    See https://gitlab.com/machielsimonbos/hectorp and estimatetrend wiki.

    Parameters
    ----------
    component : str, optional
        Component(s) to detrend ('N', 'E', 'U' or 'NEU'). Default is 'NEU'.
    **kwargs : dict, optional
        Additional arguments for the Hector control file (estimatetrend).

    Returns
    -------
    Gts
        Residual time series with velocity (and sigmas), offsets attributes.
    """
    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect
    import os
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # check if hectorp is installed by running os command estimatetrend
    import subprocess
    result = subprocess.run(['estimatetrend', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        ERROR("hectorp is not installed. Please install it from https://gitlab.com/machielsimonbos/hectorp", exit=True)

    # function to convert a Gts object to a HectorP .mom format 
    def gts2mom(gts):
        import pyacs.lib.astrotime as at
        import numpy as np

        date_mom = at.decyear2mjd(gts.data[:, 0])
        # force mid-day (integer MJD day)
        date_mom = np.floor(date_mom).astype(int) + 0.5
        pos_mom = gts.data[:, 1:4]*1.E3

        # header
        header = "sampling period 1.000000"
        # offsets
        if gts.offsets_dates is not None :
            for date in gts.offsets_dates:
                header += "\noffset %f " % (at.decyear2mjd(date))

        # create output directory if it does not exist
        import os
        if not os.path.exists('tmp_hectorp'):
            os.makedirs('tmp_hectorp')

        # north
        out_mom_north = 'tmp_hectorp/'+gts.code + "_north.mom"
        np.savetxt(out_mom_north, np.c_[date_mom, pos_mom[:, 0]], fmt="%10.5f %15.6f", header=header)
        # east
        out_mom_east = 'tmp_hectorp/'+gts.code + "_east.mom"
        np.savetxt(out_mom_east, np.c_[date_mom, pos_mom[:, 1]], fmt="%10.5f %15.6f", header=header)
        # up
        out_mom_up = 'tmp_hectorp/'+gts.code + "_up.mom"
        np.savetxt(out_mom_up, np.c_[date_mom, pos_mom[:, 2]], fmt="%10.5f %15.6f", header=header)  

    # parse estimatetrend report
    def parse_estimatetrend_report(text: str) -> Dict[str, Any]:
        """
        Parse 'estimatetrend' text output and extract:
        - noise_models: list of dicts
        - bias: dict(value, sigma, epoch)
        - trend: dict(value, sigma, unit)
        - seasonal: dict(period -> {cos/sin/amp/pha: {value, sigma, unit/deg}})
        - offsets: list of dicts(epoch, value, sigma, unit)
        """


        import re
        from dataclasses import dataclass, asdict
        from pathlib import Path
        from typing import Any, Dict, List, Optional, Union

        @dataclass
        class NoiseModel:
            name: str
            fraction: Optional[float] = None
            sigma: Optional[float] = None
            sigma_unit: Optional[str] = None
            d: Optional[float] = None
            kappa: Optional[float] = None
            one_minus_phi: Optional[float] = None


        def _to_float(x: str) -> float:
            # Handles "  3.0183", "-1.0000", "0.0000", etc.
            return float(x.strip())

        out: Dict[str, Any] = {
            "noise_models": [],
            "bias": None,
            "trend": None,
            "seasonal": {},  # keyed by period (float)
            "offsets": [],
        }

        # -------- Noise Models block --------
        # We'll locate "Noise Models" section and parse model sub-blocks like:
        # FlickerGGM:
        # fraction  = 0.36643
        # sigma     =  3.0183 mm/yr^0.25
        # d         =  0.5000 (fixed)
        # kappa     = -1.0000 (fixed)
        # 1-phi     =  0.0000 (fixed)
        #
        # White:
        # fraction  = 0.63357
        # sigma     =  0.9079 mm
        # ...
        noise_section_match = re.search(
            r"(?ms)^Noise Models\s*\n-+\s*\n(?P<body>.*?)(?:\n\s*\n|^bias\s*:|^trend\s*:|^offset\b|$)",
            text,
        )
        if noise_section_match:
            body = noise_section_match.group("body")

            # Split into model blocks by headers "Name:"
            # Keep header lines using a regex finditer.
            model_header_re = re.compile(r"(?m)^(?P<name>[A-Za-z0-9_\-]+)\s*:\s*$")
            headers = list(model_header_re.finditer(body))

            for i, h in enumerate(headers):
                name = h.group("name")
                start = h.end()
                end = headers[i + 1].start() if i + 1 < len(headers) else len(body)
                block = body[start:end]

                nm = NoiseModel(name=name)

                m = re.search(r"(?m)^\s*fraction\s*=\s*(?P<v>[+-]?\d+(?:\.\d+)?)\s*$", block)
                if m:
                    nm.fraction = _to_float(m.group("v"))

                m = re.search(
                    r"(?m)^\s*sigma\s*=\s*(?P<v>[+-]?\d+(?:\.\d+)?)\s*(?P<u>.+?)\s*$",
                    block,
                )
                if m:
                    nm.sigma = _to_float(m.group("v"))
                    nm.sigma_unit = m.group("u").strip()

                m = re.search(r"(?m)^\s*d\s*=\s*(?P<v>[+-]?\d+(?:\.\d+)?)\b", block)
                if m:
                    nm.d = _to_float(m.group("v"))

                m = re.search(r"(?m)^\s*kappa\s*=\s*(?P<v>[+-]?\d+(?:\.\d+)?)\b", block)
                if m:
                    nm.kappa = _to_float(m.group("v"))

                m = re.search(r"(?m)^\s*1-phi\s*=\s*(?P<v>[+-]?\d+(?:\.\d+)?)\b", block)
                if m:
                    nm.one_minus_phi = _to_float(m.group("v"))

                out["noise_models"].append(asdict(nm))

        # -------- bias --------
        # bias : 4.320 +/- 0.946 (at 55965.50)
        m = re.search(
            r"(?m)^\s*bias\s*:\s*(?P<v>[+-]?\d+(?:\.\d+)?)\s*\+/-\s*(?P<s>[+-]?\d+(?:\.\d+)?)\s*\(at\s*(?P<t>[+-]?\d+(?:\.\d+)?)\)\s*$",
            text,
        )
        if m:
            out["bias"] = {"value": _to_float(m.group("v")), "sigma": _to_float(m.group("s")), "epoch": _to_float(m.group("t"))}

        # -------- trend --------
        # trend: 0.955 +/- 0.164 mm/days
        m = re.search(
            r"(?m)^\s*trend\s*:\s*(?P<v>[+-]?\d+(?:\.\d+)?)\s*\+/-\s*(?P<s>[+-]?\d+(?:\.\d+)?)\s*(?P<u>.+?)\s*$",
            text,
        )
        if m:
            out["trend"] = {"value": _to_float(m.group("v")), "sigma": _to_float(m.group("s")), "unit": m.group("u").strip()}

        # -------- seasonal terms --------
        # cos  365.250 : -0.387 +/- 0.139 mm
        # sin  365.250 : -0.158 +/- 0.141 mm
        # amp  365.250 : 0.442 +/- 0.135 mm
        # pha  365.250 : -111.046 +/- 108.700 degrees
        seasonal_re = re.compile(
            r"(?m)^\s*(?P<kind>cos|sin|amp|pha)\s+(?P<per>\d+(?:\.\d+)?)\s*:\s*(?P<v>[+-]?\d+(?:\.\d+)?)\s*\+/-\s*(?P<s>[+-]?\d+(?:\.\d+)?)\s*(?P<u>mm|degrees)\s*$"
        )
        for m in seasonal_re.finditer(text):
            kind = m.group("kind")
            per = float(m.group("per"))
            v = _to_float(m.group("v"))
            s = _to_float(m.group("s"))
            u = m.group("u").strip()

            if per not in out["seasonal"]:
                out["seasonal"][per] = {}
            out["seasonal"][per][kind] = {"value": v, "sigma": s, "unit": u}

        # -------- offsets --------
        # offset at 56129.5000 :   -4.77 +/-  0.61 mm
        offset_re = re.compile(
            r"(?m)^\s*offset\s+at\s+(?P<t>[+-]?\d+(?:\.\d+)?)\s*:\s*(?P<v>[+-]?\d+(?:\.\d+)?)\s*\+/-\s*(?P<s>[+-]?\d+(?:\.\d+)?)\s*(?P<u>mm)\s*$"
        )
        for m in offset_re.finditer(text):
            out["offsets"].append(
                {"epoch": _to_float(m.group("t")), "value": _to_float(m.group("v")), "sigma": _to_float(m.group("s")), "unit": m.group("u").strip()}
            )

        return out    


    def run_hectorp(gts, component='ENU', **kwargs):
        import os
        import subprocess
        from icecream import ic
        import numpy as np
        import pyacs.lib.astrotime as at

        import logging
        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        # check tmp_hector exists
        if not os.path.exists('tmp_hectorp'):
            os.makedirs('tmp_hectorp')

        # new data
        new_data = np.zeros((gts.data.shape[0],3))

        # create control file
        base = gts.code
        default_ctl = {
            "DataFile": f"{base}_north.mom",
            "DataDirectory": os.getcwd()+'/tmp_hectorp',
            "OutputFile": os.getcwd()+'/tmp_hectorp/'+f"{base}_north.out",
            "TimeUnit": "days",
            "interpolate": "no",
            "PhysicalUnit": "mm",
            "ScaleFactor": 1.0,
            "periodicsignals": "365.25 182.625",
            "DegreePolynomial": 1,
            "NoiseModels": "GGM White",
            "GGM_1mphi": "6.9e-06",
            "useRMLE" :  "no",
            "Verbose": "yes"
        }

        # update with user-provided parameters
        default_ctl.update(kwargs)
        link_path = 'estimatetrend.ctl'

        # North component
        if 'N' in component:
            default_ctl["DataFile"] = f"{base}_north.mom"
            default_ctl["OutputFile"] = os.getcwd()+'/tmp_hectorp/'+f"{base}_north.out"
            # write control file for noth component
            ctl_file = 'tmp_hectorp/'+gts.code + "_north.ctl"
            if os.path.exists(ctl_file):
                os.remove(ctl_file)
            with open(ctl_file, "w") as f:
                for k, v in default_ctl.items():
                    f.write(f"{k} {v}\n")
            # run HectorP command
            # Remove existing file/symlink (including broken symlink)
            if os.path.lexists(link_path):
                os.unlink(link_path)
            os.symlink(ctl_file, link_path)
            # run HectorP command
            cmd = ['estimatetrend']
            VERBOSE("Running HectorP for north component (%d data) with Noise Models = %s" % (gts.data.shape[0], default_ctl["NoiseModels"]))
            completed = subprocess.run(cmd, check=True, text=True, capture_output=True)
            # get HectorP output
            results_north = parse_estimatetrend_report(completed.stdout)
            # check processing is OK
            if results_north["trend"] is None:
                ERROR("HectorP processing failed for north component. ", exit=True)
            # populate new_data
            new_data[:, 0] = np.genfromtxt(default_ctl["OutputFile"], usecols=2)*1.E-3
            # write completed.stdout to file
            with open(os.path.join('tmp_hectorp', f"{base}_north.log"), "w") as f:
                f.write(completed.stdout)


        # East component
        if 'E' in component:
            default_ctl["DataFile"] = f"{base}_east.mom"
            default_ctl["OutputFile"] = os.getcwd()+'/tmp_hectorp/'+f"{base}_east.out"
            # write control file for east component
            ctl_file = 'tmp_hectorp/'+gts.code + "_east.ctl"
            if os.path.exists(ctl_file):
                os.remove(ctl_file)
            with open(ctl_file, "w") as f:
                for k, v in default_ctl.items():
                    f.write(f"{k} {v}\n")
            # run HectorP command
            # Remove existing file/symlink (including broken symlink)
            if os.path.lexists(link_path):
                os.unlink(link_path)
            os.symlink(ctl_file, link_path)
            cmd = ['estimatetrend']
            VERBOSE("Running HectorP for  east component (%d data) with Noise Models = %s" % (gts.data.shape[0], default_ctl["NoiseModels"]))
            completed = subprocess.run(cmd, check=True, text=True, capture_output=True)
            # get HectorP output
            results_east = parse_estimatetrend_report(completed.stdout)
            if results_east["trend"] is None:
                ERROR("HectorP processing failed for east component. ", exit=True)
            # populate new_data
            new_data[:, 1] = np.genfromtxt(default_ctl["OutputFile"], usecols=2)*1.E-3
            # write completed.stdout to file
            with open(os.path.join('tmp_hectorp', f"{base}_east.log"), "w") as f:
                f.write(completed.stdout)

        # Up component
        if 'U' in component:
            default_ctl["DataFile"] = f"{base}_up.mom"
            default_ctl["OutputFile"] = os.getcwd()+'/tmp_hectorp/'+f"{base}_up.out"
            # write control file for up component
            ctl_file = 'tmp_hectorp/'+gts.code + "_up.ctl"
            if os.path.exists(ctl_file):
                os.remove(ctl_file)
            with open(ctl_file, "w") as f:
                for k, v in default_ctl.items():
                    f.write(f"{k} {v}\n")
            # run HectorP command
            # Remove existing file/symlink (including broken symlink)
            if os.path.lexists(link_path):
                os.unlink(link_path)
            os.symlink(ctl_file, link_path)
            cmd = ['estimatetrend']
            VERBOSE("Running HectorP for    up component (%d data) with Noise Models = %s" % (gts.data.shape[0], default_ctl["NoiseModels"]))
            completed = subprocess.run(cmd, check=True, text=True, capture_output=True)
            # get HectorP output
            results_up = parse_estimatetrend_report(completed.stdout)
            if results_up["trend"] is None:
                ERROR("HectorP processing failed for up component. ", exit=True)
            # populate new_data
            new_data[:, 2] = np.genfromtxt(default_ctl["OutputFile"], usecols=2)*1.E-3
            # write completed.stdout to file
            with open(os.path.join('tmp_hectorp', f"{base}_up.log"), "w") as f:
                f.write(completed.stdout)

        # Populate and return residual gts
        new_gts = gts.copy()
        new_gts.data[:, 1:4] = new_gts.data[:, 1:4] - new_data
        # robust velocity assignment: keep existing values when available, fill only processed components
        if getattr(new_gts, "velocity", None) is not None and len(new_gts.velocity) >= 6:
            velocity = np.array(new_gts.velocity[:6], dtype=float)
        else:
            velocity = np.full(6, np.nan, dtype=float)

        if "results_north" in locals() and results_north.get("trend") is not None:
            velocity[0] = results_north["trend"].get("value", np.nan)
            velocity[3] = results_north["trend"].get("sigma", np.nan)

        if "results_east" in locals() and results_east.get("trend") is not None:
            velocity[1] = results_east["trend"].get("value", np.nan)
            velocity[4] = results_east["trend"].get("sigma", np.nan)

        if "results_up" in locals() and results_up.get("trend") is not None:
            velocity[2] = results_up["trend"].get("value", np.nan)
            velocity[5] = results_up["trend"].get("sigma", np.nan)

        new_gts.velocity = velocity * 1.E-3 # convert from mm/yr to m/yr

        # offsets
        if gts.offsets_dates is not None:
            np_offsets_values = np.zeros((len(gts.offsets_dates), 7)) # epoch + 3 components + 3 sigmas
            if "results_north" in locals() and results_north.get("offsets") != []:
                for i, offset in enumerate(results_north["offsets"]):
                    np_offsets_values[i, 0] = at.mjd2decyear(offset.get("epoch", np.nan))
                    np_offsets_values[i, 1] = offset.get("value", np.nan) * 1.E-3 # convert from mm to m
                    np_offsets_values[i, 4] = offset.get("sigma", np.nan) * 1.E-3 # convert from mm to m
            if "results_east" in locals() and results_east.get("offsets") != []:
                for i, offset in enumerate(results_east["offsets"]):
                    np_offsets_values[i, 0] = at.mjd2decyear(offset.get("epoch", np.nan))
                    np_offsets_values[i, 2] = offset.get("value", np.nan) * 1.E-3 # convert from mm to m
                    np_offsets_values[i, 5] = offset.get("sigma", np.nan) * 1.E-3 # convert from mm to m
            if "results_up" in locals() and results_up.get("offsets") != []:
                for i, offset in enumerate(results_up["offsets"]):
                    np_offsets_values[i, 0] = at.mjd2decyear(offset.get("epoch", np.nan))
                    np_offsets_values[i, 3] = offset.get("value", np.nan) * 1.E-3 # convert from mm to m
                    np_offsets_values[i, 6] = offset.get("sigma", np.nan) * 1.E-3 # convert from mm to m
            new_gts.offsets_values = np_offsets_values

        # annual and semi annual terms
        np_annual = np.zeros(6)
        np_semi_annual = np.zeros(6)
        if "results_north" in locals() and results_north.get("seasonal") != []:
            amp = results_north["seasonal"].get(365.25, {}).get("amp", {}).get("value", np.nan) 
            pha = results_north["seasonal"].get(365.25, {}).get("pha", {}).get("value", np.nan)
            np_annual[0] = amp * np.cos(pha/180.*np.pi) * 1.E-3 # convert from mm to m
            np_annual[1] = amp * np.sin(pha/180.*np.pi) * 1.E-3 # convert from mm to m

        if "results_east" in locals() and results_east.get("seasonal") != []:
            amp = results_east["seasonal"].get(365.25, {}).get("amp", {}).get("value", np.nan) 
            pha = results_east["seasonal"].get(365.25, {}).get("pha", {}).get("value", np.nan)
            np_annual[2] = amp * np.cos(pha/180.*np.pi) * 1.E-3 # convert from mm to m
            np_annual[3] = amp * np.sin(pha/180.*np.pi) * 1.E-3 # convert from mm to m
        if "results_up" in locals() and results_up.get("seasonal") != []:
            amp = results_up["seasonal"].get(365.25, {}).get("amp", {}).get("value", np.nan) 
            pha = results_up["seasonal"].get(365.25, {}).get("pha", {}).get("value", np.nan)
            np_annual[4] = amp * np.cos(pha/180.*np.pi) * 1.E-3 # convert from mm to m
            np_annual[5] = amp * np.sin(pha/180.*np.pi) * 1.E-3 # convert from mm to m

        if np_annual.any() != 0:
            new_gts.annual = np_annual

        # semi_annual

        if "results_north" in locals() and results_north.get("seasonal") != []:
            amp = results_north["seasonal"].get(182.625, {}).get("amp", {}).get("value", np.nan) 
            pha = results_north["seasonal"].get(182.625, {}).get("pha", {}).get("value", np.nan)
            np_semi_annual[0] = amp * np.cos(pha/180.*np.pi) * 1.E-3 # convert from mm to m
            np_semi_annual[1] = amp * np.sin(pha/180.*np.pi) * 1.E-3 # convert from mm to m

        if "results_east" in locals() and results_east.get("seasonal") != []:
            amp = results_east["seasonal"].get(182.625, {}).get("amp", {}).get("value", np.nan) 
            pha = results_east["seasonal"].get(182.625, {}).get("pha", {}).get("value", np.nan)
            np_semi_annual[2] = amp * np.cos(pha/180.*np.pi) * 1.E-3 # convert from mm to m
            np_semi_annual[3] = amp * np.sin(pha/180.*np.pi) * 1.E-3 # convert from mm to m
        if "results_up" in locals() and results_up.get("seasonal") != []:
            amp = results_up["seasonal"].get(182.625, {}).get("amp", {}).get("value", np.nan) 
            pha = results_up["seasonal"].get(182.625, {}).get("pha", {}).get("value", np.nan)
            np_semi_annual[4] = amp * np.cos(pha/180.*np.pi) * 1.E-3 # convert from mm to m
            np_semi_annual[5] = amp * np.sin(pha/180.*np.pi) * 1.E-3 # convert from mm to m

        if np_semi_annual.any() != 0:
            new_gts.semi_annual = np_semi_annual

        # populate the noise attribute - todo



        return new_gts

    # MAIN
    gts2mom(self)
    new_gts = run_hectorp(self, component=component, **kwargs)
    return new_gts