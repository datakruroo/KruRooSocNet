% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ego_metrics.R
\name{analyze_categorical_similarity}
\alias{analyze_categorical_similarity}
\title{Analyze homophily (ego-alter similarity) for categorical attributes}
\usage{
analyze_categorical_similarity(
  adj_matrix,
  node_attributes,
  category_col,
  directed = TRUE,
  self_ties = FALSE
)
}
\arguments{
\item{adj_matrix}{An adjacency matrix representing network ties}

\item{node_attributes}{A data frame containing node attributes}

\item{category_col}{The name of the categorical attribute column to analyze (e.g., "department")}

\item{directed}{Logical, whether the network is directed. Default is TRUE}

\item{self_ties}{Logical, whether to include self ties. Default is FALSE}
}
\value{
A data frame with homophily measures including counts, proportions, EI index,
odds ratio, log odds ratio, and Yule's Q
}
\description{
Analyze homophily (ego-alter similarity) for categorical attributes
}
\examples{
# Create sample adjacency matrix
friendship_matrix <- matrix(
  c(0,1,1,0,1,0,0,1,1,0,0,1),
  nrow = 4, byrow = TRUE
)
rownames(friendship_matrix) <- colnames(friendship_matrix) <- paste0("P", 1:4)

# Create sample attributes
attributes <- data.frame(
  node = paste0("P", 1:4),
  department = c("HR", "IT", "HR", "Finance")
)

# Analyze homophily based on department
analyze_categorical_similarity(
  friendship_matrix,
  attributes,
  "department"
)
}
