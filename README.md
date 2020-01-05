# Rough Heston
Implements the rough Heston model of (). Letting x

![$dS(t) = S(t)(r-q)dt + S(t)\sqrt{V(t)}dW_1(t)$](https://render.githubusercontent.com/render/math?math=%24dS(t)%20%3D%20S(t)(r-q)dt%20%2B%20S(t)%5Csqrt%7BV(t)%7DdW_1(t)%24)

![$V(t) = V(0) + \int_0^t K(t-s) \lambda(\bar{v} - V(s))ds + \int_0^t K(t-s) \xi \sqrt{V(s)}dW_2(s)$](https://render.githubusercontent.com/render/math?math=%24V(t)%20%3D%20V(0)%20%2B%20%5Cint_0%5Et%20K(t-s)%20%5Clambda(%5Cbar%7Bv%7D%20-%20V(s))ds%20%2B%20%5Cint_0%5Et%20K(t-s)%20%5Cxi%20%5Csqrt%7BV(s)%7DdW_2(s)%24)

where 

![K(t) = \frac{1}{\Gamma(\alpha)} t^{1-\alpha}](https://render.githubusercontent.com/render/math?math=K(t)%20%3D%20%5Cfrac%7B1%7D%7B%5CGamma(%5Calpha)%7D%20t%5E%7B1-%5Calpha%7D)

and where we restrict 

![$\alpha \in \[1/2,1\], \rho < 0, (\lambda,\bar{v},v_0)$](https://render.githubusercontent.com/render/math?math=%24%5Calpha%20%5Cin%20%5B1%2F2%2C1%5D%2C%20%5Crho%20%3C%200%2C%20(%5Clambda%2C%5Cbar%7Bv%7D%2Cv_0)%24)>0.


