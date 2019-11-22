# SBCK (Statistical Bias Correction Kit)

## Features
- python3 and R version
- c++ independent files for Sparse Histogram
- Implement classic methods of bias correction (see [1,6] for the definition of bias correction)
- Quantile Mapping [2,3,4] and CDFt methods [5] 
- OTC and dOTC methods [6]

## Python instruction

Requires:
- python3
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- numpy
- scipy
- pybind11

For python, just use the command:
```
python3 setup.py install --user
```

If the Eigen library is not found, use:
```
python3 setup.py install --user eigen="path-to-eigen"
```

## R instruction

Requires:
- R
- roxygen2
- devtools
- Rcpp
- RcppEigen
- methods
- R6

Just run:
```
Rscript build.R -c -v -i
```


## Examples

For bias correction example, X0 and X1 are respectively the random variable to correct in calibration and
projection period. Y0 is the reference during calibration period. Z0 and Z1 are the corrections during calibration
and projection period.


### Univariate bias correction with QM and CDFt
![Alt](/figures/qm_cdft.png)

### Bivariate bias correction with OTC and dOTC
![Alt](/figures/otc_dotc.png)

### Sparse histogram of a 2d-normal law with 10^6 samples
![Alt](/figures/SparseHist.png)


## Licence

Copyright Yoann Robin, 2019

This software is a computer program that is part of the SBCK (Statistical
Bias Correction Kit). This library makes it possible to perform bias
correction with non parametric methods, and give some metrics between Sparse
Histogram is high dimensions.

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.


## References
- [1] Piani, C., Weedon, G., Best, M., Gomes, S., Viterbo, P., Hagemann, S., and Haerter, J.: Statistical bias correction of global simulated daily precipitation and temperature for the application of hydrological models, J. Hydrol., 395, 199–215, https://doi.org/10.1016/j.jhydrol.2010.10.024, 2010.
- [2] Panofsky, H. A. and Brier, G. W.: Some applications of statistics to meteorology, Mineral Industries Extension Services, College of Mineral Industries, Pennsylvania State University, 103 pp., 1958.
- [3] Wood, A. W., Leung, L. R., Sridhar, V., and Lettenmaier, D. P.: Hydrologic Implications of Dynamical and Statistical Approaches to Downscaling Climate Model Outputs, Clim. Change, 62, 189–216, https://doi.org/10.1023/B:CLIM.0000013685.99609.9e, 2004.
- [4] Déqué, M.: Frequency of precipitation and temperature extremes over France in an anthropogenic scenario: Model results and statistical correction according to observed values, Global Planet. Change, 57, 16–26, https://doi.org/10.1016/j.gloplacha.2006.11.030, 2007.
- [5] Michelangeli, P.-A., Vrac, M., and Loukos, H.: Probabilistic downscaling approaches: Application to wind cumulative distribution functions, Geophys. Res. Lett., 36, L11708, https://doi.org/10.1029/2009GL038401, 2009.
- [6] Robin, Y., Vrac, M., Naveau, P., Yiou, P.: Multivariate stochastic bias corrections with optimal transport, Hydrol. Earth Syst. Sci., 23, 773–786, 2019, https://doi.org/10.5194/hess-23-773-2019
- [8] Bazaraa, M. S., Jarvis, J. J., and Sherali, H. D.: Linear Programming and Network Flows, 4th edn., John Wiley & Sons, 2009.
- [9] Sinkhorn Distances: Lightspeed Computation of Optimal Transportation Distances. arXiv, https://arxiv.org/abs/1306.0895
- [10] Wasserstein, L. N. (1969). Markov processes over denumerable products of spaces describing large systems of automata. Problems of Information Transmission, 5(3), 47-52.

