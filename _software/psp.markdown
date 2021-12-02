---
layout: post
title: psp [aut, cre]
description:  Implements a parameter space partitioning algorithm for evaluating the global behaviour of computational models as described by Pitt, Kim, Navarro and Myung (2006)
img: /assets/img/software_markus-spiske-cvBBO4PzWPg-unsplash.jpg
importance: 1
category: cognitive-science
---
<a href="https://cran.r-project.org/package=psp"><img src="https://cranlogs.r-pkg.org/badges/grand-total/psp" height = 20px alt="" /></a>
<a href="https://cran.r-project.org/package=psp"><img src="https://img.shields.io/cran/v/psp" height = 20px alt="" /></a>
<a href="https://cran.r-project.org/package=psp"><img src="https://img.shields.io/cran/l/psp" height = 20px alt="" /></a>


## Implements an n-dimensional parameter space partitioning algorithm for evaluating the global behaviour of formal computational models as described by [Pitt, Kim, Navarro and Myung (2006)](https://psycnet.apa.org/doiLanding?doi=10.1037/0033-295X.113.1.57)

<br>

<blockquote class="twitter-tweet" data-dnt="true" data-theme="dark"><p lang="en" dir="ltr">My first <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a> <a href="https://twitter.com/hashtag/rpackage?src=hash&amp;ref_src=twsrc%5Etfw">#rpackage</a> is out on CRAN!<br>The package implements a global qualitative <a href="https://twitter.com/hashtag/modelevaluation?src=hash&amp;ref_src=twsrc%5Etfw">#modelevaluation</a><br>tool as described by Pitt, Kim, Navarro and Myung (2006). It also works in n dimensions. :)<a href="https://t.co/XFVyyNgUoi">https://t.co/XFVyyNgUoi</a><br><br>Give it a whirl with install.packages(&quot;psp&quot;)!</p>&mdash; Lénárd Döme (@lenarddome) <a href="https://twitter.com/lenarddome/status/1407269362667560960?ref_src=twsrc%5Etfw">June 22, 2021</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

<br>

This is a brief manual for those without time. A more extensive introduction is
in the works and will be hosted on the Project's GitHub Wiki.

## CODE DEVELOPMENT

A big influence on this implementation is an instantiation of the Open Models
Initiative, [catlearn](https://github.com/ajwills72/catlearn).

Watch the talk of [Andy Wills: “The OpenModels project”](https://youtu.be/SfqkqEYagJU)
from Open Research Working Group (ORWG) virtual meeting 08/09/20.

The project's architecture is also influenced by [DEoptim](https://github.com/ArdiaD/DEoptim).
`DEoptim` implements a Differential Evolutionary Optimization algorithm for
model-fitting.

We are completely open-source and free. Anyone can contribute. If you would
like to raise an issue or contribute code, use Github, message or email me
(@lenarddome).

**There is a great and concise [blog post](https://www.andywills.info/2021-06-23-psp/)
about it by [Andy Wills](https://www.andywills.info/)**, who helped me quite a bit
in understanding what psp is and how it works. He is also an author
on the package. His post makes some great points about the essence of parameter
space partitioning. I thought that I would complement that post with a
non-exhaustive manual.

In the remainder of this post, I will walk through the steps of using `psp`.
This walk-through will use a two parameter model. PSP will need to find 10 distinct regions.

## INSTALL

For the stable version:

```r
install.packages("psp")
```

For the developmental version:

```r
devtools::install_github("lenarddome/psp")
```


## MODEL

Then we put together a model (essentially a model of a polytope), that
calculates the euclidean distance from a selected number of points. These
points are selected at random at the beginning and kept constant throughout
the simulation. In the original paper,
[Pitt, Kim, Navarro and Myung (2006)](https://psycnet.apa.org/doi/10.1037/0033-295X.113.1.57)
used a hypercube to test the algorithm, here I will use a polytope. I choose
to use a polytope, because I want the regions to vary in size, compared to a hypercube
where the space is uniformly partitioned.

```r
#' euclidean distance
#'
#' @param a vector coordinate 1
#' @param b vector coordinate 2
#' @return euclidean distance between coordinates
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# define center points for the 10 regions in a two-dimensional space
positions <- NULL
for (i in seq_len(2)) positions <- cbind(positions, sample(500, 10))
```

If we have those two, we could put together our model of a polytope.
We need to code our model, so it take sin a vector of parameters and
outputs a character vector. It doesn't have to be a character vector.
It can also be a Boolean, or integers, but I am yet to test it with those.

```r
#' dummy polytope model to test the PSP function
#' The model takes in a set of coordinates, calculates its distance from all
#' all of available coordinates, then return closest region number.
#' This model generalizes to n-dimensions
#'
#' @param x a vector of coordinates
#' @return The number of the region as character
#' @examples
#' model(runif(5))
model <- function(par) {
    areas <- NULL
    for (i in seq_along(par)) {
        range <- c(1, 0)
        if (i %% 2 == 0) {
            range <- c(0, 1)
        }
        areas <- cbind(areas,
                       seq(range[1], range[2], length.out = 500)[positions[,i]])
    }
    dist <- apply(areas, 1, function(x) euclidean(par, x))
    return(as.character(which.min(dist)))
}
```

## SIMULATION

```r
library(psp)
```

Now we can let psp do its job.
Here we run the MCMC for 400 iterations, but the partitioning
will stop if the population of all regions reach 300.
Note that we have to load our utility function into
the clusters, because psp_global will run parallel.

```r
# run Parameter Space Partitioning with some default settings
out <- psp_global(model, psp_control(lower = rep(0, 2),
                                   upper = rep(1, 2),
                                   init = rep(0.5, 2),
                                   radius = rep(0.25, 2),
                                   pop = 300,
                                   cluster_names = c("positions",
                                                     "euclidean"),
                                   iterations = 500))
```

This process produces us the following result:

```r
$ps_patterns

  1   2   3   4   5   6   7   8   9  10
300 344 317 306 359 358 307 396 416 397
```

In this case, psp_global stopped before it reached the 500th iteration, because
all regions reached at least 300 `pop`. We can also see that some regions have a
population larger than 300. This is because even though the sampling from that
regions stopped, new points can still be classed as members of those regions.

So PSP partitioned the parameter space into distinct disjointed regions, according
to what the model predicted.

This is how it looks under the hood in real time:

<iframe title="vimeo-player" src="https://player.vimeo.com/video/548351899?h=097ecfdea2" width="640" height="360" frameborder="0" allowfullscreen></iframe>

Each colour is a separate region.

Admittedly, this is a very simple use-case. Below is a more apt illustration of
what `psp` is doing behind the scenes for a 3D polytope model:

<iframe title="vimeo-player" src="https://player.vimeo.com/video/652405415?h=ffafca4001" width="640" height="426" frameborder="0" allowfullscreen></iframe>

### AMAZING, BUT WHAT NOW?

This was a simple model bearing no psychological commitment - unless you are
one of the few believing that the brain has a geometry module or that learning
happens according to laws of geometry.

The question is what we can do with the output now? A simple next step could
be **calculating volume of the regions**. Behaviours that occupy large regions
in the parameter space might be more frequent behaviours, therefore might need
more of our attentiont. Maybe the behaviour with the largest region is something
that we haven't even observed yet.

Alternatively, you can try to discover clusters of points in the parameter space
within regions and ordinal patterns. There might be many ways a model might
operate but still outputs the same result.

You can compare how many different qualitative outputs the model produces and
how many of those have been observed in humans. You might also try to figure out
when unobserved qualitative outputs occur according to the model.

## NOTES FOR THE CURIOUS

I had some thoughts while trying to implement the algorithm. Turns out,
some of my concern was picked up by people other than me.

### VOLUMES [will not be included]

Not feasible to implement a method that generalizes to n-dimensional polyhedra
or convex polytope. There are already packages out there that can do it. I would
leave it for the user. The method of calculating the volume/area of each region
should be an explicit choice the modeller makes.

### BURN-IN [will not be implemented]

If you have a decent starting point (e.g. parameters EXIT uses to produce inverse
base-rate effect, or ALCOVE best-fitting parameters for the Type I-VI problems),
burn-in is unnecessary.

I am also not sure why burn-in is necessary for parameter space partitioning.
It seems counter-intuitive to discard areas you explored in the parameter space
if you'd like to explore said parameter space. Here we have no target to reach
other than to fill in the whole space - unlike scenarios where you want
to optimize some point estimate like the mean of a distribution

One problem we might encounter is that *regions further away from our starting
jumping distribution will be under-sampled*. This could be avoided by increasing the number
of `iterations`, so the MCMC will sample long enough to adequately populate
those regions as well. One might also choose to decrease the radius to
sample from smaller areas surrounding the jumping distributions.

Resources to look through:

*   https://stats.stackexchange.com/questions/88819/mcmc-methods-burning-samples
*   http://users.stat.umn.edu/%7Egeyer/mcmc/burn.html
*   https://www.johndcook.com/blog/2016/01/25/mcmc-burn-in/
