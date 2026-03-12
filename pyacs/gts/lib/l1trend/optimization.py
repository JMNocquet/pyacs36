"""
Optimization functions for L1-trend analysis.
"""

import numpy as np

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def best_l1trend_golden(x, y, criterion_idx, bounds=[-2, 1], tol=0.01, logger=None, component_mask=None):
    """
    Find the optimal hyperparameter alpha in l1trend using golden section search algorithm.
    
    Parameters
    ----------
    x : numpy.ndarray
        Input time array
    y : numpy.ndarray
        Input data array
    criterion_idx : int
        Index of the criterion to use (-1 for BIC, -2 for AICc, -3 for Cp)
    bounds : list
        Bounds for the search [lower, upper]
    tol : float
        Tolerance for convergence
    logger : logging.Logger, optional
        Logger instance for logging messages
        
    Returns
    -------
    numpy.ndarray
        Optimally filtered data
    """
    from time import time
    from .statistics import get_stats_l1_model

    try:
        from trendfilter import trend_filter
    except:
        from pathlib import Path
        theme_path = 'bokeh_theme.yaml'
        if not theme_path.exists():
            theme_path.write_text(
                "attrs:\n"
                "    Axis:\n"
                "        major_label_text_font_size: '12pt'\n"
                "        axis_label_text_font_size: '14pt'\n"
                "    Line:\n"
                "        line_width: 3\n"
                "    Title:\n"
                "        text_font_size: '12pt'\n"
                "    Legend:\n"
                "        background_fill_alpha: 0.8\n",
                encoding='utf-8'
            )
        from trendfilter import trend_filter


    # Start timing
    start_time = time()
    msg = f"Starting golden section search optimization ({y.shape[0]} epochs)"
    VERBOSE(msg)
    if logger:
        logger.info(msg)

    # if y.shape[0] <=5, return the best fitting line
    if y.shape[0] <= 5:
        # Compute best fitting line using least squares
        coeffs = np.polyfit(x, y, 1)  # 1st degree polynomial (line)
        fy = coeffs[0] * x + coeffs[1]  # y = mx + b
        mesg = (f"Less than 5 epochs in best_l1trend_golden. Using best fitting line.")
        if logger:
            logger.warning(mesg)
        else:
            WARNING(mesg)
        return fy

    def convergence_reached(score1, score2, score_ref, cp1, cp2, cpref, ny):
        """Check if convergence criteria are met."""
        if np.fabs(score1 - score_ref) < 0.01 and np.fabs(score2 - score_ref) < 0.01:
            return True, 'Convergence reached: delta_score < tol'
        if np.fabs(cp1-cp2) < 2:
            return True, 'Convergence reached: same number of parameters'
        if (cp1 == cpref) and (cp2 == cpref):
            return True, 'Convergence reached: same number of parameters'
        if (cp1 >= ny) and (cp2 >= ny):
            return True, 'Stop: too many breakpoints'
        return False, None

    def cost(alpha, x, y, criterion_idx, component_mask):
        """Compute cost function for given alpha."""
        try:
            fy = trend_filter(x, y, l_norm=1, alpha_2=np.float_power(10, alpha))['y_fit']
            H_alpha = get_stats_l1_model(x, y, fy, np.float_power(10, alpha), component_mask)
            return H_alpha[criterion_idx], H_alpha[-5], fy
        except Exception as e:
            warning_msg = f"Error in cost function for alpha={alpha}: {str(e)}"

            if logger:
                logger.warning(warning_msg)
            else:
                WARNING(warning_msg)

            return np.inf, 0, np.zeros_like(y)

    def golden_section_search(a, b, x, y, criterion_idx, component_mask, maxiter=40, tol=1e-5):
        """Perform golden section search optimization."""
        ny = y.shape[0]
        golden_ratio = (np.sqrt(5) - 1) / 2
        a_start, b_start = a, b

        # Validate bounds
        try:
            score_a, cpa, fa = cost(a, x, y, criterion_idx, component_mask)
            score_b, cpb, fb = cost(b, x, y, criterion_idx, component_mask)
        except Exception as e:
            warning_msg = f"Failed to validate bounds: {str(e)}"
            WARNING(warning_msg)
            if logger:
                logger.warning(warning_msg)
            return np.zeros_like(y)

        # Main optimization loop
        niter = 0
        while niter < maxiter:
            niter += 1
            c = b - golden_ratio * (b - a)
            d = a + golden_ratio * (b - a)

            try:
                score_c, cpc, fc = cost(c, x, y, criterion_idx, component_mask)
                score_d, cpd, fd = cost(d, x, y, criterion_idx, component_mask)

                # Check convergence
                END, message = convergence_reached(score_c, score_d, min(score_c, score_d), cpc, cpd, min(cpc, cpd), ny)
                if END:
                    VERBOSE(message)
                    if logger:
                        logger.info(message)
                    best_idx = np.argmin([score_c, score_d])
                    return [fc, fd][best_idx]

                # Update search interval
                if score_c < score_d:
                    b = d
                    score_d = score_c
                    cpd = cpc
                else:
                    a = c
                    score_c = score_d
                    cpc = cpd

            except Exception as e:
                warning_msg = f"Error in iteration {niter}: {str(e)}"
                WARNING(warning_msg)
                if logger:
                    logger.warning(warning_msg)
                break

        warning_msg = f"No convergence after {maxiter} iterations"
        WARNING(warning_msg)
        if logger:
            logger.warning(warning_msg)
        return fc if score_c < score_d else fd

    try:
        # Run golden section search
        result = golden_section_search(bounds[0], bounds[1], x, y, criterion_idx, component_mask, maxiter=40, tol=tol)
        
        # Log completion
        elapsed = time() - start_time
        VERBOSE(f"Golden section search completed in {elapsed:.2f} seconds")
        if logger:
            logger.info(f"Golden section search completed in {elapsed:.2f} seconds")
        
        return result

    except Exception as e:
        warning_msg = f"Golden section search failed: {str(e)}"
        WARNING(warning_msg)
        if logger:
            logger.warning(warning_msg)
        return np.zeros_like(y)


def best_l1trend_custom(x, y, criterion_idx, logger=None, component_mask=None):
    """
    Find the optimal hyperparameter alpha in l1trend using a custom search algorithm.
    
    Parameters
    ----------
    x : numpy.ndarray
        Input time array
    y : numpy.ndarray
        Input data array
    criterion_idx : int
        Index of the criterion to use (-1 for BIC, -2 for AICc, -3 for Cp)
    logger : logging.Logger, optional
        Logger instance for logging messages
        
    Returns
    -------
    tuple
        (optimal filtered data, history dictionary, optimal alpha)
    """
    # Start timing
    start_time = time()
    VERBOSE("Starting custom l1trend optimization")
    if logger:
        logger.info("Starting custom l1trend optimization")

    # Initialize variables
    H_alpha = {}
    step = 1
    maxiter = 40
    niter = 0
    alpha = 0

    try:
        # Start with alpha = 0
        fy = trend_filter(x, y, l_norm=1, alpha_2=np.float_power(10, alpha))['y_fit']
        H_alpha[alpha] = get_stats_l1_model(x, y, fy, np.float_power(10, alpha), component_mask)
        alpha_ref = alpha
        score_ref = H_alpha[alpha][criterion_idx]
        cpref = H_alpha[alpha][-5]
        fy_ref = fy

        # Choose initial search direction based on number of change points
        if cpref < 0.1 * fy.shape[0]:
            step = -1

        alpha = alpha + step

        # Main optimization loop
        while niter < maxiter:
            niter += 1
            try:
                fy = trend_filter(x, y, l_norm=1, alpha_2=np.float_power(10, alpha))['y_fit']
                H_alpha[alpha] = get_stats_l1_model(x, y, fy, np.float_power(10, alpha), component_mask)
                score = H_alpha[alpha][criterion_idx]
                cp = H_alpha[alpha][-5]

                if score < score_ref:
                    if np.fabs(score - score_ref) < 0.01:
                        VERBOSE('Convergence reached: delta_score < tol')
                        if logger:
                            logger.info('Convergence reached: delta_score < tol')
                        break
                    if cp == cpref:
                        VERBOSE('Convergence reached: number of parameters is stable')
                        if logger:
                            logger.info('Convergence reached: number of parameters is stable')
                        break
                    alpha_ref = alpha
                    score_ref = score
                    cpref = cp
                    fy_ref = fy
                else:
                    step = -step
                    if alpha_ref + step in H_alpha:
                        step = step / 2

                alpha = alpha_ref + step

            except Exception as e:
                warning_msg = f"Error in iteration {niter}: {str(e)}"
                WARNING(warning_msg)
                if logger:
                    logger.warning(warning_msg)
                break

        # Log results
        elapsed = time() - start_time
        VERBOSE(f"Custom optimization completed in {elapsed:.2f} seconds")
        VERBOSE(f"Results: niter={niter}, alpha={alpha_ref:.3f}, cp={cpref}, score={score_ref:.2f}, ndata={x.shape[0]}")
        if logger:
            logger.info(f"Custom optimization completed in {elapsed:.2f} seconds")
            logger.info(f"Results: niter={niter}, alpha={alpha_ref:.3f}, cp={cpref}, score={score_ref:.2f}, ndata={x.shape[0]}")

        return fy_ref, H_alpha, alpha_ref

    except Exception as e:
        warning_msg = f"Custom optimization failed: {str(e)}"
        WARNING(warning_msg)
        if logger:
            logger.warning(warning_msg)
        return np.zeros_like(y), {}, 0
