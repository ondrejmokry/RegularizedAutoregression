This folder contains `.mat` files with numerical results of the experiments presented in the paper, as well as the `.m` scripts used to plot the results.

The following list describes the individual plotting scripts:
- `acceleration_test_plot_01.m` reproduces Fig. 9 in the article (scatter plots),
- `acceleration_test_plot_02.m` reproduces Fig. 10 and 11 in the article (bar plots),
- `demo_times_plot.m` reproduces Fig. 13 in the article,
- `iteration_tradeoff_plot.m` reproduces Fig. 2 in the article,
- `oracle_test_plot.m` reproduces Fig. 3 in the article,
- `survey_test_plot.m` reproduces Fig. 4, 5, 6 and 8 in the article, depending on the choice of metrics and data,
- `survey_test_plot_CR.m` reproduces the supplementary figures on the [accompanying webpage](https://ondrejmokry.github.io/RegularizedAutoregression/), featuring also the replace-reliable post-processing.

For the axes annotations/legend entries to be displayed correctly, set the interpreter to `'latex'`, e.g. using the following function:

```matlab
function setInterpreter(interpreter)
    set(groot, 'defaultAxesTickLabelInterpreter', interpreter);
    set(groot, 'defaultColorbarTickLabelInterpreter', interpreter);
    set(groot, 'defaultLegendInterpreter', interpreter);
    set(groot, 'defaultTextInterpreter', interpreter);
end
```