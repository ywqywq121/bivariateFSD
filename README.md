# A consistent semi-parametric test for bivariate stochastic dominance

This repository contains R programs for the article “A consistent semi-parametric test for bivariate stochastic dominance.”
<!-- This article has been submitted for publication. -->
## Part 1. Reproducing simulation results in Section 5.1 of the manuscript
This part compares three methods for estimating the bivariate distribution function: the empirical estimator, kernel density estimator, and our semiparametric estimator introduced in Section 2.2 of the manuscript. The ```R``` codes to reproduce the results are attached: [Estimation_coparision.R](https://github.com/ywqywq121/bivariateFSD/blob/master/Estimation_coparision.R).
### 1.1. Common copulas
We first generated data using common copulas to evaluate the performance of these methods across different dependency structures. We considered four different types of copulas: the Frank copula, Clayton copula, Gumbel copula, and Joe copula. The results are shown in Tables 1-4.

![image](https://github.com/user-attachments/assets/19f591f3-9c13-4e78-9ea2-31bc5e3ca2e5)
![image](https://github.com/user-attachments/assets/87357f63-75bd-4c36-9b78-6980589f25d0)
![image](https://github.com/user-attachments/assets/ec64a0ef-c09f-4365-afc5-6df69103a51a)
![image](https://github.com/user-attachments/assets/53b6de87-e825-4c24-96ec-609a55f2c1be)
### 1.2. Mixed copula
Secondly, we accessed the performance of our approach when working with data that doesn't originate from a candidate bivariate copula included in our model fitting, which means the wrong copulas are assumed. The results are shown in Table 5.
![image](https://github.com/user-attachments/assets/3b4d963a-d22d-4781-adfa-2904c61732bb)

### 1.3. Bivariate normal
Lastly, we consider the underlying distribution to be a bivariate normal distribution. The results are shown in Table 6.
![image](https://github.com/user-attachments/assets/26b846de-b76a-471e-a8af-5375364e217f)

## Part 2. Reproducing simulation results in Section 5.2 of the manuscript
We now illustrate the power of our proposed test.
### 2.1 Choose the order
Our sensitivity analysis aimed to investigate how the test is influenced by different orders of the tunning parameter.
![image](https://github.com/user-attachments/assets/d1c44d44-7585-4b5e-bd52-7c007ffa3ba1)
![image](https://github.com/user-attachments/assets/5d6930a5-b140-4a5b-987a-ea508d8265f7)

### 2.2 Power comparison
With the order chosen, we proceeded to compare the performance of our method with that of the method proposed in [Mccaig and Yatchew (2007)](../master/README.md#reference).

![image](https://github.com/user-attachments/assets/bfcefa90-9b28-40f0-88e4-4c613dc5ec66)

## Reference: 
1. Mccaig, B. and A. Yatchew (2007). [](https://onlinelibrary.wiley.com/doi/10.1002/jae.964)
   *Journal of Applied Econometrics*.  22(5), 951–969.
