# Rough Heston
This project contains a Matlab implementation of the rough Heston model of (). Letting S(t) denote the time t price of the underlying asset the model states the following under the risk-neutral measure 

![$dS(t) = S(t)(r-q)dt + S(t)\sqrt{V(t)}dW_1(t)$](https://render.githubusercontent.com/render/math?math=%24dS(t)%20%3D%20S(t)(r-q)dt%20%2B%20S(t)%5Csqrt%7BV(t)%7DdW_1(t)%24)

![$V(t) = V(0) + \int_0^t K(t-s) \lambda(\bar{v} - V(s))ds + \int_0^t K(t-s) \xi \sqrt{V(s)}dW_2(s)$](https://render.githubusercontent.com/render/math?math=%24V(t)%20%3D%20V(0)%20%2B%20%5Cint_0%5Et%20K(t-s)%20%5Clambda(%5Cbar%7Bv%7D%20-%20V(s))ds%20%2B%20%5Cint_0%5Et%20K(t-s)%20%5Cxi%20%5Csqrt%7BV(s)%7DdW_2(s)%24)

where r is the risk-free interest rate, q the dividend yield and where

![K(t) = \frac{1}{\Gamma(\alpha)} t^{1-\alpha}](https://render.githubusercontent.com/render/math?math=K(t)%20%3D%20%5Cfrac%7B1%7D%7B%5CGamma(%5Calpha)%7D%20t%5E%7B1-%5Calpha%7D), ![dW_1dW_2=\rho dt](https://render.githubusercontent.com/render/math?math=dW_1dW_2%3D%5Crho%20dt), ![$\alpha \in \[1/2,1\], \rho < 0, (\lambda,\bar{v},v_0)$](https://render.githubusercontent.com/render/math?math=%24%5Calpha%20%5Cin%20%5B1%2F2%2C1%5D%2C%20%5Crho%20%3C%200%2C%20(%5Clambda%2C%5Cbar%7Bv%7D%2Cv_0)%24)>0.

The implementation is based on Fourier pricing methods as suggested in e.g. (Gerhold et al., 2019) with an optimal integration contour
as proposed in (Lord & Kahl, 2006). To compute the characteristic function we solve the Volterra integral equation which appears in e.g. () and that using the scheme from (Diethelm, 2004).



Below we illustrate a few smiles under the model:

![Rough Heston Smiles](https://github.com/sigurdroemer/rough_heston/blob/master/smiles2.jpg)

The parameters are 

![v_0=\bar{v}=0.15^2,\alpha=0.6,\lambda=2,\xi=0.4,\rho=-0.6](https://render.githubusercontent.com/render/math?math=v_0%3D%5Cbar%7Bv%7D%3D0.15%5E2%2C%5Calpha%3D0.6%2C%5Clambda%3D2%2C%5Cxi%3D0.4%2C%5Crho%3D-0.6)

and we have defined log-moneyness := log(strike/forward).



See the file 'get_started.m' to get started using the Matlab code yourself. In the folder '' you will find further scripts validiting the code and experimenting with various settings.



References:
  - Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment explosions in the rough Heston model, Decisions in Economics and Finance (2019) 42:575-608.
  - Roger lord and Christian Kahl, Optimal Fourier inversion in semi-analytical option pricing, available at SRRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=921336, 2006.
  - Kai Diethelm, Neville J. Ford and Alan D. Freed, Detailed error analysis for a fractional ADAMs method, Numerical Algorithms 36:31-52, 2004.
