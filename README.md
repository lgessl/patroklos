# patroklos <img src="man/figures/logo.png" align="right" height="139" alt="" />

  <!-- badges: start -->

  [![R-CMD-check](https://github.com/lgessl/patroklos/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lgessl/patroklos/actions/workflows/R-CMD-check.yaml)

  <!-- badges: end -->

## An R package pipelining binary classification of survival

patroklos is pipeline that allows you to define, train, validate, test and analyze models whose 
final goal it is to predict if an event occurs before or after a pivotal time. It does so for 
the famous train-validate-test split of your data, where you train and test your models on a 
training cohort and pick the best validated model and take it over to the test cohort to assess 
its generalization error on new data.

### Installation

Make sure you have the `devtools` package installed, start an R session and then type:
```R
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("lgessl/patroklos")
```

### What's the scenario this R package is tailored for?

patroklos is a good choice for the just mentioned scenario and it is very good choice for the 
scenario we tailored it for.

Imagine you are given a data set holding information on cancer patients. This data includes a many, 
many features from a high-throughput measurement like gene-expression levels for hundreds or even 
thousands of genes ("expression data") and some other features ("pheno data"). Among the pheno data
is survival information in the form of time to event and right censoring. You want to build a 
binary classifier that takes a subset of your features (excluding the two survival features) as 
input to predict for a given sample if it faces the event (like progression or death) before or 
after a certain time. 

This tasks involves

- preprocessing the data,
- deciding on an error measure for validation and testing,
- training a whole bunch of models and validating them, including tuning their hyperparameters,
- assessing the best validated model on the test cohort and 
- unlocking the test set for all trained models to do some analysis, e.g., how reliably 
- cross-validated errors estimate the test error.

At any point in this sequence, you or your boss may come up with new ideas: training a new model, 
training the same models on a different data set, testing the models on a different test set or 
changing the error measure for validation. But you are the boss of this machine learning project
and, as such, you need to be flexible and able to integrate new ideas quickly. You may hire a 
secretary. Or you resort to patroklos.

patroklos will do a big part of the repetitive, administrative work and lets you focus on the 
statistical and biological part of your work. If there is still bureaucratic work left for you, 
patroklos at least tells you how to do it. 

All in all, patroklos is not the star of the shows. But it plays a key role in it — just like 
Πάτροκλος in the Iliad: being Achilles' best friend, his death brings Achilles back to the 
battlefield helping the Greeks score decisive victories.

### How does this R package make your life easier?

patroklos has three major goals: making integrating 

- new data,
- new models with new hyperparameters and 
- your custom validation and testing

into your workflow as effortless as possible and at the same leaving you enough freedom 
to adapt the pipeline to the needs of your project.

#### Data

The R6 class `Data` tells you how to preprocess your data. patroklos also gives you some tools 
you often need for preprocessing.

#### Models

patroklos abstracts the data from the models, which you specify in another R6 class, `Model`. 
A model typically specifies more than just one combination of hyperparameters that need to be tuned.
Concerning model classes and fitting functions, you have full flexibility — you can provide any 
fitting function as an attribute to a `Model` as long as its function interface fulfills patroklos's 
expectations. If the fitting function generates hyperparameters automatically or tunes 
hyperparameters internally, that's totally fine as long it returns a single fit object with 
validated predictions and a predict method. Again, it is left 100% to the model how it calculates 
the validated predictions, e.g. by cross validation or out-of-bag predictions. patroklos also 
provides patroklos-compliant fitting functions out of the box, namely for Gauss, logistic and Cox 
models by means of the `zeroSum` package and random forests by means of the `ranger` package. 
For every `Model`, patroklos can tune some model-agnostic hyperparameters itself.

#### Validation and testing

Validate your models on their training data and test them on some test data for a variety of 
metrics including (scalar) error measures — precision (under the constraint of a minimum 
rate of positive predictions a.k.a prevalence), hazard ratio, ROC-AUC, lower and upper bound of the 
$\gamma$-confidence interval of the precision or hazard ratio — with the R6 classes `AssScalar`.
For many models you first need to threshold their output before you have a binary classifier at 
hand. The R6 class `Ass2d`, which, for every possible threshold, plots one metric of binary 
classifier against another — e.g. prevalence against precision — guides in thresholding. In testing, 
you often want to compare your picked model to a benchmark model — no problem, just include its 
output into your data and specify a `Model` for it with the help of the `projection_on_feature()` 
fitter.

 #### Meta analysis

After testing the picked model on a test cohort, unlock it for more models, plot validation versus 
test error, group models by their hyperparameters and find out which hyperparameters come with 
systematic flaws in validation and notoriously high test errors. Exclude problematic models from 
the set of `Model`s in the future. `val_vs_test()` is your friend here.

### More tweaks by patroklos

#### Storing and reading in models

Training and assessing so many models takes a lot of time. So as your project evolves over time, 
you don't want to fit your models again and again. To this end, patroklos heavily resorts to 
storing the models and reading them in on demand. Typically you would also "store" the `Model`, 
`AssScalar` and `Ass2d` objects, either as an .rds file or, more typically, their initialization in 
an R script. 

#### Integrating models into another model

Some of the features of a `Model` may be the output of another model, e.g., from a gene-expression
signature. We call the latter model, whose output becomes the new feature, an *early* and the 
former model, which experiences the output of the early model as a feature, the *late* model. 
If the early model was trained on another data set, this is just an ordinary feature 
(which might have issues with batch effects) — easy. If, however, you want to train the early 
model on the same data set as the late model, things become tricky and you might want to have a 
look at `greedy_nestor()` and `long_nestor()`.


### Usage

See patroklos in action in the 
[repository of my master thesis](https://github.com/lgessl/master-thesis).

### Using your own models

To fit into the pipeline, models and their related functions must meet certain 
requirements. Typically in R, training and assessing a model involves 

- a function that fit_obj the model to training data, optionally tunes hyperparameters
  in a cross validation, e.g., and returns 
- an S3 object representing the fit. Often this object holds some kind of validated 
  predictions like those from a cross validation or out-of-bag predictions in the 
  case of a random forest. This S3 object has a 
- `predict()` method that predicts on new data.

We will refer to them as *fitter*, *fit object* and *predict method* in the 
following. 

For many packages, these three components already fulfill the expectations 
patroklos has or you can at least use part of patroklos's functionality with them.
To be able to access patroklos's full power for a model, you may need to wrap 
the fitter and the predict method and to modify the fit object to make them 
patroklos-compatible. There are multiple definitions of patroklos-compatibility 
we will elaborate on below.

#### patroklos-compliant predict method

A *patroklos-compliant predict method* declares three parameters,

- `object`: the fit object,
- `newx`: prediction data as a named numeric matrix, samples correspond to named
  rows, features to named columns, and
- `...`: additional parameters (often passed on to the wrapped predict method).

It returns a named numeric vector of predictions from `newx`.

#### patroklos-compliant fitter

A *patroklos-compliant fitter* declares three parameters,

- `x`: training data as a named numeric matrix, samples correspond to named rows,
  features to named columns. It has an attribute `li_var_suffix`, i.e. columns in 
  `x` with the suffix `li_var_suffix` are designated for late integration and the 
  fitter may use this information.
- `y`: binary response as a named numeric vector. Its names must match 
  the row names of `x`.
- `...`: additional parameters (often passed on to the wrapped fitter).

It returns an S3 object with a patroklos-compliant predict method.

#### patroklos-compliant fitter with validated predictions

A *patroklos-compliant fitter with validated predictions* is a patroklos-compliant
fitter whose S3 return value has an attribute `val_predict`, a 
numeric vector holding some form of validated predictions. Validated predictions 
are predictions made on independent data like cross-validated (CV) predictions 
or out-of-bag (OOB) predictions in the case of random forests.

#### patroklos-compliant fitter with CV tuning 

A *patroklos-compliant fitter with cross-validation tuning* is a
patroklos-compliant fitter whose return value has the attributes

- `val_predict_list`: list of numeric vectors holding cross-validated 
predictions for every hyperparameter $\lambda$,
- `lambda`: numeric or character vector of the hyperparameters $\lambda$.

While a *patroklos-compliant fitter with validated predictions* only performs 
a cross validation for one combination of hyperparameters, a *patroklos-compliant
fitter with CV tuning* performs a cross validation for multiple values of a 
hyperparameter $\lambda$ (which may be a tuple, not just a scalar). The embraced 
word "scalar" is a clear constraint here and it may be worth relaxing the definition 
from a search line to a search grid for hyperparameters in the future (which of 
course involves modifying those functions requiring a *patroklos-compliant 
fitter with CV tuning*).