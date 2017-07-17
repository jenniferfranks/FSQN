\name{Feature Specific Quantile Normalization}
\alias{quantileNormalizeByFeature}
\title{Feature Specific Quantile Normalization - FSQN}
\description{Returns the normalized matrix.}
\examples{
target <- matrix(rnorm(100*150, mean = 1, sd = 1), nrow = 100, ncol = 150)
test <- matrix(rnorm(30*150, mean = 2, sd = 2), nrow = 30, ncol = 150)

normalized_test <- quantileNormalizeByFeature(test, target)

par(mfrow=c(1,2))
plot(density(target), main = "Pre-normalization", xlim = c(-6,6),
    ylim = c(0,.6), lwd=2)
lines(density(test), col = "red", lwd=2)

plot(density(target), main = "Post-normalization", xlim = c(-6,6),
    ylim = c(0,1), lwd=2)
lines(density(normalized_test), col = "red", lwd=2)
}
