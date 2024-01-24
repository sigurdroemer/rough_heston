# Rough Heston
This project implements the pricing of European calls and puts under the rough Heston model of (El Euch & Rosenbaum, 2018) and (El Euch & Rosenbaum, 2019).

## Implementation
The Matlab implementation is based on numerical integration with the Fourier transform as suggested in e.g. (Gerhold et al., 2019) and with an optimal integration contour as proposed in (Lord & Kahl, 2006). To compute the characteristic function, we solve the Volterra integral equation that appears in e.g. (El Euch & Rosenbaum, 2019) with the scheme of (Diethelm, 2004). The implementation assumes negative spot-vol correlation.

See the file 'get_started.m' for a few examples on how to use the code. In the folder 'validation_and_test_scripts' you will find further scripts that validate the code and experiment with different settings.

Remark: The code was developed in Matlab 2019a and may not work in older versions.

## Illustration

Below we illustrate a few smiles under the model:

![Rough Heston Smiles](https://github.com/sigurdroemer/rough_heston/blob/master/smile2.jpg)

## References
- Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment explosions in the rough Heston model, Decisions in Economics and Finance (2019) 42:575-608.
- Omar El Euch and Mathieu Rosenbaum, Perfect hedging in rough Heston models, Mathematical Finance 28(6), 3813-3856, 2018.
- Omar El Euch and Mathieu Rosenbaum, The characteristic function of rough Heston models, Mathematical Finance 29(1), 3-38, 2019.
- Roger lord and Christian Kahl, Optimal Fourier inversion in semi-analytical option pricing, Tinbergen Institute Discussion Paper No. 2006-066/2. Available at SRRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=921336.
- Kai Diethelm, Neville J. Ford and Alan D. Freed, Detailed error analysis for a fractional ADAMs method, Numerical Algorithms 36:31-52, 2004.
