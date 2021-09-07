# SBCK (Statistical Bias Correction Kit)

## Features
- python3 and R version
- c++ independent files for Sparse Histogram
- Implement classic methods of bias correction (see [8,9] for the definition of bias correction)
- Quantile Mapping [5,7,14], parametric and non parametric version
- CDFt methods [6] 
- OTC and dOTC methods [9]
- R2D2 method [11]
- MBCn method [4]
- QDM method [3]
- MRec method [1]
- ECBC method [12]
- TSMBC method [15], for autocorrelations.

## How to select a bias correction method ?

This summary of ability of each method to perform a bias correction is proposed by François, (2020). Please refer to
this article for further interpretation.

| Characteristics                             | CDF-t              | R2D2               | dOTC               | MBCn               | MRec               |
|---------------------------------------------| :----------------: | :----------------: | :----------------: | :----------------: | :----------------: |
| Correction of univariate dist. prop.        | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| Modification of correlations of the model   | :x:                | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| Capacity to correct inter-var. prop.        | :x:                | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| Capacity to correct spatial prop.           | :x:                | :heavy_check_mark: | :heavy_check_mark: | :warning:          | :warning:          |
| Capacity to correct temporal prop.          | :x:                | :x:                | :x:                | :x:                | :x:                |
| Preserve the rank structure of the model    | :heavy_check_mark: | :warning:          | :warning:          | :warning:          | :warning:          |
| Capacity to correct small geographical area | n.a.               | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| Capacity to correct large geographical area | n.a.               | :warning:          | :warning:          | :warning:          | :x:                |
| Allow for change of the multi-dim. prop.    | :heavy_check_mark: | :x:                | :heavy_check_mark: | :warning:          | :heavy_check_mark: |


## R instruction

SBCK is available in [CRAN](https://cran.r-project.org/web/packages/SBCK/index.html).

Requires:
- R
- roxygen2 (>= 7.0.0)
- devtools
- Rcpp
- RcppEigen
- methods
- R6
- [ROOPSD](https://github.com/yrobink/ROOPSD)

Just run:
```
Rscript build.R -c -v -i
```

## License

Copyright(c) 2021 Yoann Robin

This file is part of SBCK.

SBCK is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SBCK is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SBCK.  If not, see <https://www.gnu.org/licenses/>.


## References
- [[1]](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011WR011524) Bárdossy, A. and Pegram, G.: Multiscale spatial recorrelation of RCM precipitation to produce unbiased climate change scenarios over large areas and small, Water Resources Research, 48, 9502–, https://doi.org/10.1029/2011WR011524, 2012.
- [2] Bazaraa, M. S., Jarvis, J. J., and Sherali, H. D.: Linear Programming and Network Flows, 4th edn., John Wiley & Sons, 2009.
- [[3]](https://doi.org/10.1175/JCLI-D-14-00754.1) Cannon, A. J., Sobie, S. R., and Murdock, T. Q.: Bias correction of simulated precipitation by quantile mapping: how well do methods preserve relative changes in quantiles and extremes?, J. Climate, 28, 6938–6959, https://doi.org/10.1175/JCLI-D-14-00754.1, 2015.
- [[4]](https://link.springer.com/article/10.1007/s00382-017-3580-6) Cannon, Alex J.: Multivariate quantile mapping bias correction: an N-dimensional probability density function transform for climate model simulations of multiple variables, Climate Dynamics, nb. 1, vol. 50, p. 31-49, 10.1007/s00382-017-3580-6, 2018.
- [[5]](https://doi.org/10.1016/j.gloplacha.2006.11.030) Déqué, M.: Frequency of precipitation and temperature extremes over France in an anthropogenic scenario: Model results and statistical correction according to observed values, Global Planet. Change, 57, 16–26, https://doi.org/10.1016/j.gloplacha.2006.11.030, 2007.
- [[6]](https://doi.org/10.1029/2009GL038401) Michelangeli, P.-A., Vrac, M., and Loukos, H.: Probabilistic downscaling approaches: Application to wind cumulative distribution functions, Geophys. Res. Lett., 36, L11708, https://doi.org/10.1029/2009GL038401, 2009.
- [7] Panofsky, H. A. and Brier, G. W.: Some applications of statistics to meteorology, Mineral Industries Extension Services, College of Mineral Industries, Pennsylvania State University, 103 pp., 1958.
- [[8]](https://doi.org/10.1016/j.jhydrol.2010.10.024) Piani, C., Weedon, G., Best, M., Gomes, S., Viterbo, P., Hagemann, S., and Haerter, J.: Statistical bias correction of global simulated daily precipitation and temperature for the application of hydrological models, J. Hydrol., 395, 199–215, https://doi.org/10.1016/j.jhydrol.2010.10.024, 2010.
- [[9]](https://doi.org/10.5194/hess-23-773-2019) Robin, Y., Vrac, M., Naveau, P., Yiou, P.: Multivariate stochastic bias corrections with optimal transport, Hydrol. Earth Syst. Sci., 23, 773–786, 2019, https://doi.org/10.5194/hess-23-773-2019
- [[10]](https://arxiv.org/abs/1306.0895) Sinkhorn Distances: Lightspeed Computation of Optimal Transportation Distances. arXiv, https://arxiv.org/abs/1306.0895
- [[11]](https://doi.org/10.5194/hess-22-3175-2018) Vrac, M.: Multivariate bias adjustment of high-dimensional climate simulations: the Rank Resampling for Distributions and Dependences (R2 D2 ) bias correction, Hydrol. Earth Syst. Sci., 22, 3175–3196, https://doi.org/10.5194/hess-22-3175-2018, 2018.
- [[12]](https://doi.org/10.1175/JCLI-D-14-00059.1) Vrac, M. and P. Friederichs, 2015: Multivariate—Intervariable, Spatial, and Temporal—Bias Correction. J. Climate, 28, 218–237, https://doi.org/10.1175/JCLI-D-14-00059.1
- [13] Wasserstein, L. N. (1969). Markov processes over denumerable products of spaces describing large systems of automata. Problems of Information Transmission, 5(3), 47-52.
- [[14]](https://doi.org/10.1023/B:CLIM.0000013685.99609.9e) Wood, A. W., Leung, L. R., Sridhar, V., and Lettenmaier, D. P.: Hydrologic Implications of Dynamical and Statistical Approaches to Downscaling Climate Model Outputs, Clim. Change, 62, 189–216, https://doi.org/10.1023/B:CLIM.0000013685.99609.9e, 2004.
- [[15]](https://doi.org/10.5194/esd-2021-12) Robin, Y. and Vrac, M.: Is time a variable like the others in multivariate statistical downscaling and bias correction?, Earth Syst. Dynam. Discuss. [preprint], https://doi.org/10.5194/esd-2021-12, in review, 2021.
- François, B., Vrac, M., Cannon, A., Robin, Y., and Allard, D.: Multivariate bias corrections of climate simulations: Which benefits for which losses?, Earth Syst. Dyn., 11, 537–562, https://doi.org/10.5194/esd-11-537-2020, https://esd.copernicus.org/articles/11/537/2020/, 2020.

