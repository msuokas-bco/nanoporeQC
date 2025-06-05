# R/utils-globals.R

# Suppress R CMD check "no visible binding for global variable" notes
# related to non-standard evaluation (e.g., in ggplot2::aes)
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "accuracy"
    # Add any other variables here that cause similar warnings
    # e.g., "frequency", "value", etc. if they appear in your ggplot calls
  ))
}
