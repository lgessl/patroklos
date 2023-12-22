# lymphomaSurvivalPipeline

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/lgessl/lymphomaSurvivalPipeline/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lgessl/lymphomaSurvivalPipeline/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

## Pipelining supervised learning for gene-expression-based lymphoma survival analysis made easy

The lymphomaSurvivalPipeline provides the infrastructure to learn multiple supervised models predicting progression-free survival for cancer patients and finally pick a best models. It helps you

- preprocess your data making sure it has the right format for this pipeline,
- split it into a training and test cohort,
- prepare the data for model training and predicting,
- fit a variety of models (at the moment typically Cox proportional-hazards and binomial regression models at the core), and
- assess an arbitrary subset of your models.

Fitting a model includes storing it so you can later use it in a late-integration paradigm. This implies you can add new models at any time and again run (new) assessments.

The last step, assessment, aims for finding a model that best characterizes a high-risk and low-risk group of patients in terms of progression-free survival. You can also add a benchmark like the International Prognostic Index (IPI) for non-Hogdekin lymphoma.

Adding new data and doing exactly the same as for previous data sets, i.e., fitting the same models and doing the same assessments, can be done with minimal effort.

## Minimal coding effort

Using the lymphomaSurvivalPipeline means minimal coding effort for you since it takes care of all the routine steps of the machine-learning project. You just need to call a small number of functions–and tell it what to do in a quite abstract way, namely using S3 classes specifying

- a data set (like the names of certain columns/features): `DataSpec`,
- a model (like the fitter function, which features to add to the predictor variables from the pheno data): `ModelSpec`, and
- a plot assessing a model's performance (e.g. that you want to plot the rate of positive predictions against the precision): `PerfPlotSpec`.

## Read the docs

We will add a quick start and examples on how to use the lymphomaSurvivalPipeline (as vignettes and hopefully on a static webpage) in the future once there is time. For now, browsing the exported functions below [`R`](R) and reading their documentation is the best way to start. 