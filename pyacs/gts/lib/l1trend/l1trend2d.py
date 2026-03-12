"""
2D L1-trend filtering using CVXPY.

This module provides 2D L1-trend filtering functionality with support for
modern CVXPY solvers including the new default Clarabel solver.
"""

import numpy as np
import cvxpy as cp


def second_diff_matrix(n):
    """
    Generate the second difference matrix for L1-trend filtering.
    
    Parameters
    ----------
    n : int
        Size of the matrix
        
    Returns
    -------
    numpy.ndarray
        Second difference matrix of shape (n-2, n)
    """
    D = np.zeros((n-2, n))
    for i in range(n-2):
        D[i, i] = 1
        D[i, i+1] = -2
        D[i, i+2] = 1
    return D


def l1_trendfilter_2d(Y, lam=1.0, solver="CLARABEL", verbose=False, scs_kwargs=None, clarabel_kwargs=None, warm_start=None):
    Y = np.asarray(Y, dtype=float)
    T = Y.shape[0]
    if Y.ndim != 2 or Y.shape[1] != 2:
        raise ValueError("Y must have shape (T, 2).")

    D2 = second_diff_matrix(T)
    X  = cp.Variable((T, 2))
    data_fit = 0.5 * cp.sum_squares(X - Y)
    D2X = D2 @ X
    group_pen = cp.sum(cp.norm(D2X, axis=1))
    obj = cp.Minimize(data_fit + lam * group_pen)
    prob = cp.Problem(obj)

    # Warm-start (helps SCS a lot)
    if warm_start is not None:
        try:
            X.value = np.asarray(warm_start, float)
        except Exception:
            pass

    if solver.upper() == "SCS":
        opts = dict(max_iters=5000, eps=1e-4, scale=1.0)
        if scs_kwargs: opts.update(scs_kwargs)
        prob.solve(solver="SCS", verbose=verbose, **opts)
    elif solver.upper() == "CLARABEL":
        opts = dict(max_iter=10000)
        if clarabel_kwargs: opts.update(clarabel_kwargs)
        prob.solve(solver="CLARABEL", verbose=verbose, **opts)
    elif solver.upper() == "ECOS":
        # ECOS is deprecated, use CLARABEL instead
        import warnings
        warnings.warn("ECOS solver is deprecated. Using CLARABEL instead.", FutureWarning)
        opts = dict(max_iter=10000)
        if clarabel_kwargs: opts.update(clarabel_kwargs)
        prob.solve(solver="CLARABEL", verbose=verbose, **opts)
    else:
        prob.solve(solver=solver, verbose=verbose)

    return prob, X, D2


def find_optimal_lambda(Y, lambda_range=None, method='cross_validation', cv_folds=5, solver="CLARABEL", verbose=False):
    """
    Find the optimal lambda value for L1-trend filtering using various criteria.
    
    Parameters
    ----------
    Y : array_like
        Input data of shape (T, 2)
    lambda_range : array_like, optional
        Range of lambda values to test. If None, uses a logarithmic range.
    method : str, optional
        Method to use for optimization. Options: 'cross_validation', 'aic', 'bic', 'gcv'
    cv_folds : int, optional
        Number of folds for cross-validation
    solver : str, optional
        Solver to use for optimization
    verbose : bool, optional
        Whether to print progress information
        
    Returns
    -------
    float
        Optimal lambda value
    dict
        Dictionary containing optimization results
    """
    Y = np.asarray(Y, dtype=float)
    T = Y.shape[0]
    
    if lambda_range is None:
        # Create a logarithmic range of lambda values
        lambda_range = np.logspace(-2, 3, 20)  # From 0.001 to 100
    
    if method == 'cross_validation':
        return _cross_validation_lambda(Y, lambda_range, cv_folds, solver, verbose)
    elif method == 'aic':
        return _aic_lambda(Y, lambda_range, solver, verbose)
    elif method == 'bic':
        return _bic_lambda(Y, lambda_range, solver, verbose)
    elif method == 'gcv':
        return _gcv_lambda(Y, lambda_range, solver, verbose)
    else:
        raise ValueError("Method must be one of: 'cross_validation', 'aic', 'bic', 'gcv'")


def _cross_validation_lambda(Y, lambda_range, cv_folds, solver, verbose):
    """Cross-validation approach for lambda selection."""
    T = Y.shape[0]
    fold_size = T // cv_folds
    cv_scores = []
    
    for lam in lambda_range:
        if verbose:
            print(f"Testing lambda = {lam:.4f}")
        
        fold_scores = []
        for fold in range(cv_folds):
            # Create train/test split
            test_start = fold * fold_size
            test_end = (fold + 1) * fold_size if fold < cv_folds - 1 else T
            
            # Create mask for training data
            train_mask = np.ones(T, dtype=bool)
            train_mask[test_start:test_end] = False
            
            # Train on training data
            Y_train = Y[train_mask]
            if len(Y_train) < 3:  # Need at least 3 points for second differences
                continue
                
            try:
                prob, X, D2 = l1_trendfilter_2d(Y_train, lam=lam, solver=solver, verbose=False)
                if prob.status == 'optimal':
                    # Predict on test data
                    X_pred = X.value
                    # Interpolate predictions to test time points
                    train_indices = np.where(train_mask)[0]
                    test_indices = np.arange(test_start, test_end)
                    
                    # Simple linear interpolation for missing points
                    Y_test_pred = np.zeros((len(test_indices), 2))
                    for j in range(2):  # For each component
                        Y_test_pred[:, j] = np.interp(test_indices, train_indices, X_pred[:, j])
                    
                    # Calculate prediction error
                    Y_test_true = Y[test_start:test_end]
                    mse = np.mean((Y_test_true - Y_test_pred) ** 2)
                    fold_scores.append(mse)
            except:
                continue
        
        if fold_scores:
            cv_scores.append(np.mean(fold_scores))
        else:
            cv_scores.append(np.inf)
    
    # Find optimal lambda
    optimal_idx = np.argmin(cv_scores)
    optimal_lambda = lambda_range[optimal_idx]
    
    results = {
        'lambda_range': lambda_range,
        'cv_scores': cv_scores,
        'optimal_lambda': optimal_lambda,
        'method': 'cross_validation'
    }
    
    if verbose:
        print(f"Optimal lambda (CV): {optimal_lambda:.4f}")
    
    return optimal_lambda, results


def _aic_lambda(Y, lambda_range, solver, verbose):
    """AIC (Akaike Information Criterion) approach for lambda selection."""
    T = Y.shape[0]
    aic_scores = []
    
    for lam in lambda_range:
        if verbose:
            print(f"Testing lambda = {lam:.4f}")
        
        try:
            prob, X, D2 = l1_trendfilter_2d(Y, lam=lam, solver=solver, verbose=False)
            if prob.status == 'optimal':
                # Calculate residuals
                residuals = Y - X.value
                rss = np.sum(residuals ** 2)
                
                # Count effective degrees of freedom (approximate)
                # For L1 trend filtering, this is roughly the number of breakpoints + 2
                D2X = D2 @ X.value
                n_breakpoints = np.sum(np.abs(D2X) > 1e-6)
                df = n_breakpoints + 2  # +2 for intercept and slope
                
                # AIC = 2k - 2ln(L), where k is number of parameters, L is likelihood
                # For Gaussian errors: AIC = 2k + n*ln(RSS/n)
                aic = 2 * df + T * np.log(rss / T)
                aic_scores.append(aic)
            else:
                aic_scores.append(np.inf)
        except:
            aic_scores.append(np.inf)
    
    # Find optimal lambda
    optimal_idx = np.argmin(aic_scores)
    optimal_lambda = lambda_range[optimal_idx]
    
    results = {
        'lambda_range': lambda_range,
        'aic_scores': aic_scores,
        'optimal_lambda': optimal_lambda,
        'method': 'aic'
    }
    
    if verbose:
        print(f"Optimal lambda (AIC): {optimal_lambda:.4f}")
    
    return optimal_lambda, results


def _bic_lambda(Y, lambda_range, solver, verbose):
    """BIC (Bayesian Information Criterion) approach for lambda selection."""
    T = Y.shape[0]
    bic_scores = []
    
    for lam in lambda_range:
        if verbose:
            print(f"Testing lambda = {lam:.4f}")
        
        try:
            prob, X, D2 = l1_trendfilter_2d(Y, lam=lam, solver=solver, verbose=False)
            if prob.status == 'optimal':
                # Calculate residuals
                residuals = Y - X.value
                rss = np.sum(residuals ** 2)
                
                # Count effective degrees of freedom
                D2X = D2 @ X.value
                n_breakpoints = np.sum(np.abs(D2X) > 1e-6)
                df = n_breakpoints + 2
                
                # BIC = k*ln(n) - 2ln(L)
                # For Gaussian errors: BIC = k*ln(n) + n*ln(RSS/n)
                bic = df * np.log(T) + T * np.log(rss / T)
                bic_scores.append(bic)
            else:
                bic_scores.append(np.inf)
        except:
            bic_scores.append(np.inf)
    
    # Find optimal lambda
    optimal_idx = np.argmin(bic_scores)
    optimal_lambda = lambda_range[optimal_idx]
    
    results = {
        'lambda_range': lambda_range,
        'bic_scores': bic_scores,
        'optimal_lambda': optimal_lambda,
        'method': 'bic'
    }
    
    if verbose:
        print(f"Optimal lambda (BIC): {optimal_lambda:.4f}")
    
    return optimal_lambda, results


def _gcv_lambda(Y, lambda_range, solver, verbose):
    """Generalized Cross-Validation approach for lambda selection."""
    T = Y.shape[0]
    gcv_scores = []
    
    for lam in lambda_range:
        if verbose:
            print(f"Testing lambda = {lam:.4f}")
        
        try:
            prob, X, D2 = l1_trendfilter_2d(Y, lam=lam, solver=solver, verbose=False)
            if prob.status == 'optimal':
                # Calculate residuals
                residuals = Y - X.value
                rss = np.sum(residuals ** 2)
                
                # Estimate effective degrees of freedom
                # This is a simplified approximation
                D2X = D2 @ X.value
                n_breakpoints = np.sum(np.abs(D2X) > 1e-6)
                df = n_breakpoints + 2
                
                # GCV = RSS / (n - df)^2
                gcv = rss / (T - df) ** 2
                gcv_scores.append(gcv)
            else:
                gcv_scores.append(np.inf)
        except:
            gcv_scores.append(np.inf)
    
    # Find optimal lambda
    optimal_idx = np.argmin(gcv_scores)
    optimal_lambda = lambda_range[optimal_idx]
    
    results = {
        'lambda_range': lambda_range,
        'gcv_scores': gcv_scores,
        'optimal_lambda': optimal_lambda,
        'method': 'gcv'
    }
    
    if verbose:
        print(f"Optimal lambda (GCV): {optimal_lambda:.4f}")
    
    return optimal_lambda, results


def plot_lambda_selection(results, title=None):
    """
    Plot lambda selection results.
    
    Parameters
    ----------
    results : dict
        Results from find_optimal_lambda
    title : str, optional
        Title for the plot
    """
    import matplotlib.pyplot as plt
    
    lambda_range = results['lambda_range']
    method = results['method']
    
    if method == 'cross_validation':
        scores = results['cv_scores']
        ylabel = 'Cross-Validation MSE'
    elif method == 'aic':
        scores = results['aic_scores']
        ylabel = 'AIC'
    elif method == 'bic':
        scores = results['bic_scores']
        ylabel = 'BIC'
    elif method == 'gcv':
        scores = results['gcv_scores']
        ylabel = 'GCV'
    
    plt.figure(figsize=(10, 6))
    plt.semilogx(lambda_range, scores, 'b-', linewidth=2)
    plt.axvline(results['optimal_lambda'], color='r', linestyle='--', 
                label=f'Optimal λ = {results["optimal_lambda"]:.4f}')
    plt.xlabel('Lambda (λ)')
    plt.ylabel(ylabel)
    plt.title(title or f'Lambda Selection using {method.upper()}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
