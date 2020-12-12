# Bret Larget
# October 28, 2010
# September 30, 2011


# Use source() to read this file into R to use it.
# Example:
#   # First, change the working directory in R to the location where the file gtest.R resides.
#   # Then, use source() to read in this file.
#   source("gtest.R")
#   # Create a matrix with data for this test.
#   fish = matrix(c(1,49,10,35,37,9),nrow=2,ncol=3)
#   rownames(fish) = c("Eaten","Not Eaten")
#   colnames(fish) = c("Uninfected","Lightly Infected","Highly Infected")
#   # Call the function.
#   g.test(fish)


# Function to carryout G-test for contingency table
# x is a matrix with (non-negative) integer counts

g.test = function(x,print=T) {
  # find the row and column sums of x with apply
  # apply(x,1,f) applies function f to each row (first dimension) of x
  # apply(x,2,f) applies function f to each column (second dimension) of x
  row.sum = apply(x,1,sum)
  col.sum = apply(x,2,sum)
  # n is the sum of all elements of the table
  n = sum(x)
  # x.expected[i,j] = row.sum[i] * col.sum[j] / n
  x.expected = row.sum %o% col.sum / n
  # g = 2 * sum_i sum_j ( x[i,j] * log(x[i,j]/x.expected[i,j] )
  g = 2*sum( x*log(x/x.expected) )
  # degrees of freedom is (#rows-1)*(#columns-1)
  degf = (nrow(x)-1)*(ncol(x)-1)
  # If conditions are right, p.value is from chi^2 distribution (area to right of test statistic)
  p.g = 1 - pchisq(g,degf)
  # return a list with expected counts, test statistic, and p-value
  if ( print ) {
    cat("G-Test for Contingency Tables\n\n")
    cat("Data:\n")
    print(x)
    cat("\n")
    cat("The test statistic is ",signif(g,6),".\n")
    cat("There are ",degf," degrees of freedom.\n")
    cat("The p-value is ",p.g,".\n\n\n")
  }
  return(invisible(list(expected=x.expected,test.stat=g,p.value=p.g)))
}
