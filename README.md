# patroklos <img src="man/figures/logo.png" align="right" height="139" alt="" />

  <!-- badges: start -->

  [![R-CMD-check](https://github.com/lgessl/patroklos/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lgessl/patroklos/actions/workflows/R-CMD-check.yaml)

  <!-- badges: end -->

## A pipeline to learn omics-based survival models for cancer patients

### What's the scenario this R package is made for?

Imagine you are given a data set holding information on cancer patients. This data includes a high-dimensional 
measurement like bulk RNA-seq data ("expression data") and some other features ("pheno data"). Among the pheno data
is survival information in the form of time to event and right censoring. You want to build a binary classifier 
that takes a subset of your features (excluding the two survival features) as input to predict for a given sample 
if it faces the event (like progression or death) before or after a certain time.

Finding such a classifier involves a lot of trial and error: you will try out a whole bunch of models. Every model 
family again involves some hyperparameters you need to tune. The high-dimensional part of your data requires you 
to apply regularization to your models, in other words: tuning a regularization parameter. After you've trained 
the models, you need to assess them and finally come up with the "best" one, where you first need to decide what 
you mean by "best".

Even after you've picked a winner, this is just the winner on a specific data set. Somebody (maybe even yourself) 
may come up with another data set and you need to do an analogous analysis on the new data set.

patroklos will do all the repetitive, administrative work and lets you focus on the mathematical 
and biological part of your work.

All in all, patroklos is not the star of the show, but it plays a key role in it—just like Πάτροκλος in the
Iliad: being Achilles' best friend, his death brings Achilles back to the battlefield helping the Greeks 
score decisive victories.

### How does this R package make your life easier?

patroklos has three major goals: making integrating 

- new data,
- new models and 
- new kinds of assessments 

into your workflow as easy and effortless as possible without forcing you too much to use a certain paradigm.

#### Data

We give you the tools to preprocess new data and bring it into the format patroklos works with. 
You specify and store everything in the R6 class `Data`.

#### Models

You specify models in another R6 class, `Model`; "specifying" means fixing all hyperparameters except for the 
regularization parameter $\lambda$. Why? To tune $\lambda$ we split the data into a train and test cohort 
and determine its optimal choice in a cross validation on the train cohort. The function optimizing $\lambda$ and 
all the other model parameters is an attribute of the `Model` object and you can easily provide your own ones.

#### Assessment

On the test cohort (which we haven't touched so far) we can now assess the trained models. Keep in mind that most 
models used as binary classifiers don't output the final classification, but only continuous score that 
we need to threshold (this threshold is yet another hyperparameter of the model). To find the model best 
fulfilling its job and threshold it, we proceed as follows:

1. In a pre-selection step, we map every trained model to a single real number: a score indicating its goodness 
   on the test cohort (like the ROC-AUC). Of course, you can define your metric here. This allows us to reduce 
   the models of interest to an arbitrary small number if not one. The R6 class `AssScalar` does this.
2. For the remaining models, we plot one metric against another metric (e.g. rate of positive predictions 
   versus precision). Keep in mind: we essentially get as many points in the two-dimensional space as the 
   model has different output values here as we can use every output value as a threshold. This plot or 
   these plots will guide our decision for a threshold and therefore for a final binary classifier. The 
   R6 class `Ass2d` does this.

#### In a nut shell

<img src="man/figures/README-/test-train.jpeg" alt="drawing" width="400"/>

### More tweaks powered by patroklos

#### Repeated splits into training and test cohort

Splitting your data into a training and test cohort is very close to the real-world scenario, in which the entity 
generating the data withholds the test cohort in the first place and you are only given the training cohort. 
Sometimes this split may be lucky for you, sometimes it may be really tough, but this randomness is not in your 
heads. What is in your hands though is to submit the model with the best expected (or average) performance in 
such a scenario. patroklos lets you easily repeat the training part for 
multiple splits into training and test data and averages in a reasonable manner in the assessment step.

#### Storing and reading in models

Training and assessing so many models takes a lot of time. So as your project evolves over time, you don't 
want to fit your models again and again. To this end, patroklos heavily resorts 
to storing the models and reading them in on demand. Typically you would also "store" the `Model`, `AssScalar` and 
`Ass2d` objects, either as an .rds file or, more typically, their initialization in an R script. Adding a new model 
or assessment is then equivalent to adding some lines of code to these R scripts. Typically, you would store your 
models below `models/<data-set name>`, patroklos will then store the assessments in the analogous 
file tree below `results/<data-set name>`. 

#### Early integration

For a set of diverse features (keep in mind: high-dimensional data and pheno data), it is tempting and simple 
to do early integration, i.e. providing them to a well known model with a well implemented fitting algorithm. 
While patroklos assumes you always give the model the high-dimensional omics data, 
you have the freedom to add features from pheno data, continuous and categorical ones alike. When I started 
designing this package, I had the [`zeroSum`](https://github.com/rehbergT/zeroSum) package delivering the key 
models that deal with the high-dimensional data in my mind. The `zeroSum` package provides all the functionality 
of the [`glmnet`](https://github.com/cran/glmnet) with each (generalized linear) model optionally endowed with 
the zero-sum constraint that causes the model to become scale-invariant. patroklos provides 
a wrapper around `zeroSum::zeroSum` named `zeroSumLI` to make `zeroSum` ready for early integration.

Schematically this renders as follows:

<img src="man/figures/README-/early-int.jpeg" alt="drawing" width="400"/>

#### Late integration

Early integration seems easy at first glance, but sometimes the fitting algorithm might not be able to handle 
the data you provide: features on vastly different scales. In such a case, you may want to prefer late 
integration, i.e. first training an "early" model on the high-dimensional data and then using its (scalar) output 
together with some more pheno features as the input of another, "late" model. 

I'm currently working to include late integration in patroklos and provide fitting 
algorithms to do this. 

Nesting a regularized model whose regularization is tuned in a cross validation into a another model means 
training the resulting model involves a cross validation. We should continue the cross validation of the 
early model into the late model. I made `zeroSum` report enough details on the (many, many) models fit in 
a cross validation in a fork called [`zeroSumLI`](https://github.com/lgessl/zeroSumLI). 

Nesting has no limits as the below picture shows.

<img src="man/figures/README-/late-int.jpeg" alt="drawing" width="400"/>

### Installation

Make sure you have the `devtools` package installed, start an R session and then type:
```R
devtools::install_github("lgessl/patroklos")
```

### Quick start

There are four R6 classes a user of patroklos needs to know about: `Data`, `Model`, 
`AssScalar` and `Ass2d`. To understand how they work together over the course of the pipeline and 
how to use patroklos for your project, you should read their documentation in the same 
order.

As for preprocessing, you might be interested in 

- [`discretize_tbl_columns()`] to turn continuous features in your pheno data into 
   binary ones by thresholding them,
- [`ensure_patients_math()`] to ensure expression and pheno data contain exactly 
  the same samples,
- [`qc_preprocess()`] to make sure the readily preprocessed data is in the format
   patroklos expects, and
- [`write_data_info()`] to store meta information on the data set in .json file.

`prepend_to_directory()` helps you use `Model`s you only initialized once to you 
for multiple data sets. 

### Using your own models

To fit into the pipeline, models and their related functions must meet certain 
requirements. Typically in R, training and assessing a model involves 

- a function that fits the model to training data, optionally tunes hyperparameters
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
fitter whose S3 return value has an attribute `cv_predict` or `oob_predict`, a 
numeric vector holding cross-validated or out-of-bag (OOB) predictions, 
respectively. 

#### patroklos-compliant fitter with CV tuning   

A *patroklos-compliant fitter with cross-validation tuning* is a
patroklos-compliant fitter whose return value has the attributes

- `cv_predict_list`: list of numeric vectors holding cross-validated 
predictions for every hyperparameter $\lambda$,
- `lambda`: numeric vector of the hyperparameters $\lambda$.

While a *patroklos-compliant fitter with validated predictions* only performs 
a cross validation for one combination of hyperparameters, a *patroklos-compliant
fitter with CV tuning* performs a cross validation for multiple values of a 
(scalar) hyperparameter $\lambda$. The embraced word "scalar" is a clear 
constraint here and it may be worth relaxing the definition from a search line 
to a search grid for hyperparameters in the future (which of course involves 
modifying those functions requiring a *patroklos-compliant fitter with CV tuning*).