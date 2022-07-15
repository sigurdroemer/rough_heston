# Rough Heston
This project implements the pricing of European calls and puts under the rough Heston model of (El Euch & Rosenbaum, 2018) and (El Euch & Rosenbaum, 2019). Let S(t) denote the time t price of the underlying asset. The model then assumes the following under the risk-neutral measure,

![$dS(t) = S(t)(r-q)dt + S(t)\sqrt{V(t)}dW_1(t)$](https://render.githubusercontent.com/render/math?math=%24dS(t)%20%3D%20S(t)(r-q)dt%20%2B%20S(t)%5Csqrt%7BV(t)%7DdW_1(t)%24)

![$V(t) = V(0) + \int_0^t K(t-s) \lambda(\bar{v} - V(s))ds + \int_0^t K(t-s) \xi \sqrt{V(s)}dW_2(s)$](https://render.githubusercontent.com/render/math?math=%24V(t)%20%3D%20V(0)%20%2B%20%5Cint_0%5Et%20K(t-s)%20%5Clambda(%5Cbar%7Bv%7D%20-%20V(s))ds%20%2B%20%5Cint_0%5Et%20K(t-s)%20%5Cxi%20%5Csqrt%7BV(s)%7DdW_2(s)%24)

where r is the risk-free interest rate, q the dividend yield, and

![K(t)=\frac{1}{\Gamma(\alpha)}t^{\alpha-1}, dW_1dW_2=\rho dt,  \frac{1}{2} \textless \alpha \leq 1,  -1 \leq \rho \leq 1 , (\lambda,\xi,\bar{v},V(0)) > 0](https://render.githubusercontent.com/render/math?math=K(t)%3D%5Cfrac%7B1%7D%7B%5CGamma(%5Calpha)%7Dt%5E%7B%5Calpha-1%7D%2C%20dW_1dW_2%3D%5Crho%20dt%2C%20%20%5Cfrac%7B1%7D%7B2%7D%20%5Ctextless%20%5Calpha%20%5Cleq%201%2C%20%20-1%20%5Cleq%20%5Crho%20%5Cleq%201%20%2C%20(%5Clambda%2C%5Cxi%2C%5Cbar%7Bv%7D%2CV(0))%20%3E%200)



## Implementation
The Matlab implementation is based on numerical integration with the Fourier transform as suggested in e.g. (Gerhold et al., 2019) and with an optimal integration contour as proposed in (Lord & Kahl, 2006). To compute the characteristic function, we solve the Volterra integral equation that appears in e.g. (El Euch & Rosenbaum, 2019) with the scheme of (Diethelm, 2004). The implementation assumes ![\rho < 0](https://render.githubusercontent.com/render/math?math=%5Crho%20%3C%200).

See the file 'get_started.m' for a few examples on how to use the code. In the folder 'validation_and_test_scripts' you will find further scripts that validate the code and experiment with different settings.

Remark: The code was developed in Matlab 2019a and may not work in older versions.

## Illustration

Below we illustrate a few smiles under the model:

![Rough Heston Smiles](https://github.com/sigurdroemer/rough_heston/blob/master/smile2.jpg)

The parameters are 

![v_0=\bar{v}=0.15^2,\alpha=0.6,\lambda=2,\xi=0.4,\rho=-0.6](https://render.githubusercontent.com/render/math?math=v_0%3D%5Cbar%7Bv%7D%3D0.15%5E2%2C%5Calpha%3D0.6%2C%5Clambda%3D2%2C%5Cxi%3D0.4%2C%5Crho%3D-0.6)

and we make the definition log-moneyness := log(strike/forward).

## References
- Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment explosions in the rough Heston model, Decisions in Economics and Finance (2019) 42:575-608.
- Omar El Euch and Mathieu Rosenbaum, Perfect hedging in rough Heston models, Mathematical Finance 28(6), 3813-3856, 2018.
- Omar El Euch and Mathieu Rosenbaum, The characteristic function of rough Heston models, Mathematical Finance 29(1), 3-38, 2019.
- Roger lord and Christian Kahl, Optimal Fourier inversion in semi-analytical option pricing, Tinbergen Institute Discussion Paper No. 2006-066/2. Available at SRRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=921336.
- Kai Diethelm, Neville J. Ford and Alan D. Freed, Detailed error analysis for a fractional ADAMs method, Numerical Algorithms 36:31-52, 2004.
