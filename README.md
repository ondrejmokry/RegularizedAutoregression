# RegularizedAutoregression

This is an accompanying repository for the article *Regularized autoregressive modeling and its application to audio signal declipping*, which is to be submitted to IEEE/ACM Transactions on Audio, Speech, and Language Processing.

> Autoregressive (AR) modeling is invaluable in signal processing, in particular in speech and audio fields. Attempts in the literature can be found that regularize or constrain either the time-domain signal values or the AR coefficients, which is done for various reasons, including the incorporation of prior information or numerical stabilization. Although these attempts are appealing, an encompassing and generic modeling framework is still missing. We propose such a framework and the related optimization problem and algorithm. We discuss the computational demands of the algorithm and explore the effects of various improvements on its convergence speed. In the experimental part, we demonstrate the usefulness of our approach on the audio declipping problem. We compare its performance against the state-of-the-art methods and demonstrate the competitiveness of the proposed method, especially for mildly clipped signals. The evaluation is extended by considering a heuristic algorithm of generalized linear prediction (GLP), a strong competitor which has only been presented as a patent and is new in the scientific community.

## Contents of the repository

The repository contains MATLAB implementation of all the methods and experiments described in the article.

It is organized as follows:

### Subfolders
- **results** – numerical results of the experiments presented in the paper, as well as the scripts used to plot the results
- **signals** – audio signals used in the experiments which do not use the full set from **survey toolbox**
- **survey toolbox** – clone of selected parts of the repository [declipping2020_codes](https://github.com/rajmic/declipping2020_codes), which is used for comparison of the proposed method with the state-of-the-art optimization-based audio declipping methods
- **utils** – all the functions implementing the proposed framework and functions used by the plotting scripts in **results**

### Scripts
- `acceleration_test.m` tests different acceleration options for the DRA and ACS
- `demo.m` is a demonstrative script which runs the declipping experiment with reconstruction using GLP and the proposed ACS approach
- `iteration_tradeoff.m` tests the proposed method for different combinations of the ACS (outer) and DRA (inner) iterations
- `oracle_test.m` tests the inpainting / declipping using Janssen algorithm or GLP and compares the progression of AR coefficients to the coefficients of the ground truth signal
- `survey_test.m` performs the declipping experiment from the article [A Survey and an Extensive Evaluation of Popular Audio Declipping Methods](https://ieeexplore.ieee.org/document/9281027) using the proposed method and GLP
- `survey_test_add_CR.m` performs the post-processing of the results from `survey_test.m` as described in the article [Audio Declipping Performance Enhancement via Crossfading](https://www.sciencedirect.com/science/article/pii/S0165168421004023); note that this code *cannot* be run without first running `survey_test.m` and generating all the declipped waveforms

## Dependencies

The codes were tested in MATLAB R2024a. They depend on the following toolboxes:
- Parallel Computing Toolbox,
- Signal Processing Toolbox,
- Statistics and Machine Learning Toolbox.

## Acknowledgment

Special thank you goes to:

1. the authors of the declipping survey [1] and the follow-up article [2] for sharing publically the code and numerical results, which allowed to build on their work,
2. the authors of [3] and the [Audio Inpainting Toolbox](http://small.inria.fr/keyresults/audio-inpainting/) for sharing the implementation of the basic Janssen algorithm,
3. İlker Bayram for sharing [the codes](https://web.itu.edu.tr/ibayram/Structured/) for [4].

---

[1] P. Záviška, P. Rajmic, A. Ozerov and L. Rencker, “A Survey and an Extensive Evaluation of Popular Audio Declipping Methods,” *IEEE Journal of Selected Topics in Signal Processing*, vol. 15, no. 1, pp. 5–24, 2021, doi: [10.1109/JSTSP.2020.3042071](https://doi.org/10.1109/JSTSP.2020.3042071).
[2] P. Záviška, P. Rajmic and O. Mokrý, “Audio declipping performance enhancement via crossfading,” *Signal Processing*, vol. 192, 2022, doi: [10.1016/j.sigpro.2021.108365](https://doi.org/10.1016/j.sigpro.2021.108365).
[3] A. Adler, V. Emiya, M. G. Jafari, M. Elad, R. Gribonval and M. D. Plumbley, “Audio Inpainting,” *IEEE Transactions on Audio, Speech, and Language Processing*, vol. 20, no. 3, pp. 922-932, 2012, doi: [10.1109/TASL.2011.2168211](https://doi.org/10.1109/TASL.2011.2168211).
[4] İ. Bayram, “Proximal Mappings Involving Almost Structured Matrices,” *IEEE Signal Processing Letters*, vol. 22, no. 12, pp. 2264-2268, 2015, doi: [10.1109/LSP.2015.2476381](https://doi.org/10.1109/LSP.2015.2476381).
