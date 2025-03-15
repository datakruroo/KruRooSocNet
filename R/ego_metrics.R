#' Calculate degree for each node in multiple network types
#'
#' @param network_list A list of adjacency matrices, each representing a different type of network
#' @param directed Logical, whether the networks are directed. Default is TRUE
#' @param mode Character, the type of degree to calculate: "in", "out", or "all". Default is "out"
#' @return A data frame with nodes as rows and different network types as columns
#' @export
#' @importFrom igraph graph_from_adjacency_matrix degree
#' @importFrom dplyr bind_cols
#'
#' @examples
#' # Create sample adjacency matrices
#' friendship <- matrix(c(0,1,1,0,0,0,1,0,1), nrow=3)
#' advice <- matrix(c(0,1,0,0,0,1,1,0,0), nrow=3)
#'
#' # Calculate outdegree for each network
#' networks <- list(friendship = friendship, advice = advice)
#' calculate_degree(networks)
calculate_degree <- function(network_list, directed = TRUE, mode = "out") {

  # Validate inputs
  if (!is.list(network_list)) {
    stop("network_list must be a list of adjacency matrices")
  }

  if (!all(sapply(network_list, function(x) is.matrix(x)))) {
    stop("All elements in network_list must be matrices")
  }

  # Initialize an empty list to store degree data frames
  degree_list <- list()

  # Process each network
  for (net_name in names(network_list)) {
    # Create igraph object
    g <- graph_from_adjacency_matrix(network_list[[net_name]],
                                     mode = ifelse(directed, "directed", "undirected"))

    # Calculate degree
    node_degree <- igraph::degree(g, mode = mode)

    # Create a data frame with a descriptive column name
    degree_df <- data.frame(node_degree)
    colnames(degree_df) <- paste0("degree_", net_name)

    # Add to list
    degree_list[[net_name]] <- degree_df
  }

  # Combine all degree data frames
  result <- bind_cols(degree_list)

  # Add node names as a column if available
  if (!is.null(rownames(network_list[[1]]))) {
    result$node <- rownames(network_list[[1]])
    result <- result[, c(ncol(result), 1:(ncol(result)-1))]  # Reorder columns
  } else {
    result$node <- 1:nrow(result)
    result <- result[, c(ncol(result), 1:(ncol(result)-1))]  # Reorder columns
  }

  return(result)
}


#' Calculate Agresti's Index of Qualitative Variation (IQV) for network ties
#'
#' @param degree_df A data frame containing node degrees for different relation types
#' @param node_col The name of the column containing node identifiers. Default is "node"
#'
#' @return A data frame with node identifiers and corresponding IQV values
#' @export
#' @importFrom dplyr select group_by mutate summarize ungroup
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' # Using the output from calculate_degree
#' degree_data <- data.frame(
#'   node = c("A", "B", "C"),
#'   degree_friendship = c(2, 1, 3),
#'   degree_advice = c(1, 1, 0),
#'   degree_conflict = c(0, 2, 1)
#' )
#' calculate_iqv(degree_data)
calculate_iqv <- function(degree_df, node_col = "node") {

  # Validate inputs
  if (!is.data.frame(degree_df)) {
    stop("degree_df must be a data frame")
  }

  if (!(node_col %in% colnames(degree_df))) {
    stop(paste("Column", node_col, "not found in degree_df"))
  }

  # Identify degree columns (all except the node column)
  degree_cols <- setdiff(colnames(degree_df), node_col)

  if (length(degree_cols) < 2) {
    stop("At least two degree columns are required to calculate IQV")
  }

  # Calculate IQV
  result <- degree_df %>%
    # Select relevant columns
    select(all_of(c(node_col, degree_cols))) %>%
    # Reshape to long format
    pivot_longer(cols = -all_of(node_col),
                 names_to = "relation_type",
                 values_to = "degree") %>%
    # Group by node
    group_by(.data[[node_col]]) %>%
    # Calculate proportion and IQV
    mutate(total_degree = sum(degree),
           proportion = degree / total_degree,
           Q = n()) %>%
    # Handle cases with zero total degree
    group_by(.data[[node_col]]) %>%
    summarize(
      IQV = if(sum(degree) > 0) {
        (1 - sum(proportion^2, na.rm = TRUE)) / (1 - 1/Q[1])
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  return(result)
}

#' Summarize valued tie data for each node in a network
#'
#' @param graph An igraph object with weighted edges
#' @param weight_attr The name of the edge attribute containing tie weights. Default is "weight"
#' @param mode Character, the type of edges to consider: "in", "out", or "all". Default is "all"
#'
#' @return A data frame with summary statistics for each node's tie strengths
#' @export
#' @importFrom igraph V E incident strength degree edge_attr_names
#' @importFrom stats mean sd median min max
#'
#' @examples
#' # Create a weighted network
#' library(igraph)
#' g <- graph_from_edgelist(matrix(c(1,2,1,3,2,3,2,4,3,4), ncol=2, byrow=TRUE), directed=FALSE)
#' E(g)$weight <- c(3, 1, 2, 4, 5)
#' V(g)$name <- paste0("A", 1:4)
#' summarize_valued_ties(g)
summarize_valued_ties <- function(graph, weight_attr = "weight", mode = "all") {
  # ตรวจสอบว่ามี attribute weight_attr หรือไม่
  if (!(weight_attr %in% igraph::edge_attr_names(graph))) {
    warning(paste("Edge attribute", weight_attr, "not found in the graph.",
                  "Using uniform weights (1) for all edges."))
    igraph::edge_attr(graph, weight_attr) <- rep(1, igraph::ecount(graph))
  }
  # ตรวจสอบและดึงชื่อ node
  if ("name" %in% igraph::vertex_attr_names(graph)) {
    node_names <- igraph::V(graph)$name
  } else {
    node_names <- paste0("Node", 1:igraph::vcount(graph))
  }

  # ดึงค่าน้ำหนักของ edge ทั้งหมดและคำนวณ strength ของแต่ละ node
  edge_weights <- as.numeric(igraph::edge_attr(graph, weight_attr))
  node_strengths <- igraph::strength(graph, mode = mode, weights = edge_weights)

  # คำนวณ degree (จำนวน alter) สำหรับแต่ละ node
  num_alters <- igraph::degree(graph, mode = mode)

  # เตรียม vector สำหรับเก็บ summary statistics
  n_nodes <- igraph::vcount(graph)
  avg_tie_strength    <- numeric(n_nodes)
  sd_tie_strength     <- numeric(n_nodes)
  median_tie_strength <- numeric(n_nodes)
  min_tie_strength    <- numeric(n_nodes)
  max_tie_strength    <- numeric(n_nodes)
  range_tie_strength  <- numeric(n_nodes)

  # คำนวณสถิติสำหรับแต่ละ node
  for (i in seq_len(n_nodes)) {
    # ดึง edge ที่ incident กับ node i
    incident_edges <- igraph::incident(graph, igraph::V(graph)[i], mode = mode)

    if (length(incident_edges) > 0) {
      # ใช้ edge_attr() เพื่อดึงค่าน้ำหนักของ edge เหล่านั้น
      edge_weights_i <- igraph::edge_attr(graph, weight_attr, incident_edges)

      avg_tie_strength[i]    <- mean(edge_weights_i)
      sd_tie_strength[i]     <- if (length(edge_weights_i) > 1) sd(edge_weights_i) else 0
      median_tie_strength[i] <- median(edge_weights_i)
      min_tie_strength[i]    <- min(edge_weights_i)
      max_tie_strength[i]    <- max(edge_weights_i)
      range_tie_strength[i]  <- max_tie_strength[i] - min_tie_strength[i]
    } else {
      # ถ้า node ไม่มี edge ใดๆ ให้กำหนดเป็น NA
      avg_tie_strength[i]    <- NA_real_
      sd_tie_strength[i]     <- NA_real_
      median_tie_strength[i] <- NA_real_
      min_tie_strength[i]    <- NA_real_
      max_tie_strength[i]    <- NA_real_
      range_tie_strength[i]  <- NA_real_
    }
  }

  # สร้าง data frame สรุปผลสำหรับแต่ละ node
  result <- data.frame(
    node                 = node_names,
    number_of_alters     = num_alters,
    sum_of_tie_strengths = node_strengths,
    average_tie_strength = avg_tie_strength,
    sd_tie_strength      = sd_tie_strength,
    median_tie_strength  = median_tie_strength,
    min_tie_strength     = min_tie_strength,
    max_tie_strength     = max_tie_strength,
    range_tie_strength   = range_tie_strength
  )

  return(result)
}
#' Analyze composition of alters for each ego in a network
#'
#' @param adj_matrix An adjacency matrix representing network ties
#' @param node_attributes A data frame containing node attributes
#' @param attribute_name The name of the attribute column to analyze
#' @param directed Logical, whether the network is directed. Default is TRUE
#' @param self_ties Logical, whether to include self ties. Default is FALSE
#'
#' @return A data frame with alter composition statistics for each ego
#' @export
#' @importFrom dplyr left_join group_by summarize mutate rename
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @examples
#' # Create sample data
#' adj_mat <- matrix(c(0,1,1,0,0,0,1,0,1), nrow=3)
#' rownames(adj_mat) <- colnames(adj_mat) <- c("A", "B", "C")
#' attributes <- data.frame(
#'   node = c("A", "B", "C"),
#'   department = c("HR", "IT", "Finance")
#' )
#' analyze_alter_composition(adj_mat, attributes, "department")
#' Analyze composition of alters for each ego in a network
#'
#' @param adj_matrix An adjacency matrix representing network ties
#' @param node_attributes A data frame containing node attributes
#' @param attribute_name The name of the attribute column to analyze
#' @param directed Logical, whether the network is directed. Default is TRUE
#' @param self_ties Logical, whether to include self ties. Default is FALSE
#'
#' @return A data frame with alter composition statistics for each ego
#' @export
#' @importFrom dplyr left_join group_by summarize mutate rename
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @examples
#' # Create sample data
#' adj_mat <- matrix(c(0,1,1,0,0,0,1,0,1), nrow=3)
#' rownames(adj_mat) <- colnames(adj_mat) <- c("A", "B", "C")
#' attributes <- data.frame(
#'   node = c("A", "B", "C"),
#'   department = c("HR", "IT", "Finance")
#' )
#' analyze_alter_composition(adj_mat, attributes, "department")
analyze_alter_composition <- function(adj_matrix, node_attributes, attribute_name,
                                      directed = TRUE, self_ties = FALSE) {

  # Validate inputs
  if (!is.matrix(adj_matrix)) {
    stop("adj_matrix must be a matrix")
  }

  if (!is.data.frame(node_attributes)) {
    stop("node_attributes must be a data frame")
  }

  if (!(attribute_name %in% colnames(node_attributes))) {
    stop(paste("Attribute", attribute_name, "not found in node_attributes"))
  }

  # Get node names from matrix
  if (is.null(rownames(adj_matrix))) {
    rownames(adj_matrix) <- colnames(adj_matrix) <- 1:nrow(adj_matrix)
  }

  # Prepare node id column name in attributes
  node_id_col <- colnames(node_attributes)[1]

  # Check if node names match between matrix and attributes
  if (!all(rownames(adj_matrix) %in% node_attributes[[node_id_col]])) {
    stop("Node names in adjacency matrix do not match node_attributes")
  }

  # Convert matrix to edge list format
  edges <- data.frame(
    ego = rep(rownames(adj_matrix), each = ncol(adj_matrix)),
    alter = rep(colnames(adj_matrix), times = nrow(adj_matrix)),
    tie = as.vector(adj_matrix)
  )

  # Remove self-ties if specified
  if (!self_ties) {
    edges <- edges[edges$ego != edges$alter, ]
  }

  # Keep only existing ties
  edges <- edges[edges$tie > 0, ]

  # Join with attribute data
  result <- edges %>%
    # Join attributes for ego
    left_join(node_attributes, by = c("ego" = node_id_col)) %>%
    rename(ego_attribute = !!attribute_name) %>%
    # Join attributes for alter
    left_join(node_attributes, by = c("alter" = node_id_col)) %>%
    rename(alter_attribute = !!attribute_name)

  # For categorical attribute, calculate frequency distribution
  if (is.character(result$alter_attribute) || is.factor(result$alter_attribute)) {

    # Calculate frequencies and proportions
    alter_composition <- result %>%
      group_by(ego, ego_attribute, alter_attribute) %>%
      summarize(
        frequency = sum(tie, na.rm = TRUE),
        .groups = "drop_last"
      ) %>%
      mutate(proportion = frequency / sum(frequency))

    # Calculate IQV (Index of Qualitative Variation) for each ego separately
    # Using group_by + summarize with a custom function to avoid the if condition issue
    iqv_by_ego <- alter_composition %>%
      group_by(ego) %>%
      summarize(
        num_categories = n_distinct(alter_attribute),
        sum_prop_squared = sum(proportion^2),
        IQV = ifelse(num_categories > 1,
                     (1 - sum_prop_squared) / (1 - 1/num_categories),
                     0),
        .groups = "drop"
      )

    # Join the IQV back to the main data
    alter_composition <- alter_composition %>%
      left_join(select(iqv_by_ego, ego, IQV), by = "ego")

    # Create wide format for easier viewing
    alter_composition_wide <- alter_composition %>%
      pivot_wider(
        id_cols = c(ego, ego_attribute, IQV),
        names_from = alter_attribute,
        values_from = c(frequency, proportion),
        values_fill = list(frequency = 0, proportion = 0)
      )

    return(alter_composition_wide)

  } else if (is.numeric(result$alter_attribute)) {
    # For continuous attribute, calculate summary statistics

    alter_composition <- result %>%
      group_by(ego, ego_attribute) %>%
      summarize(
        num_alters = n(),
        sum_attribute = sum(alter_attribute * tie, na.rm = TRUE),
        mean_attribute = mean(alter_attribute, na.rm = TRUE),
        sd_attribute = sd(alter_attribute, na.rm = TRUE),
        min_attribute = min(alter_attribute, na.rm = TRUE),
        max_attribute = max(alter_attribute, na.rm = TRUE),
        range_attribute = max_attribute - min_attribute,
        weighted_mean = sum(alter_attribute * tie, na.rm = TRUE) / sum(tie, na.rm = TRUE),
        .groups = "drop"
      )

    return(alter_composition)
  } else {
    stop("Attribute must be either categorical (character/factor) or continuous (numeric)")
  }
}

#' Calculate E-I index for ego-alter similarity based on categorical attributes
#'
#' @param adj_matrix An adjacency matrix representing network ties
#' @param node_attributes A data frame containing node attributes
#' @param attribute_name The name of the categorical attribute column to analyze
#' @param directed Logical, whether the network is directed. Default is TRUE
#' @param self_ties Logical, whether to include self ties. Default is FALSE
#'
#' @return A data frame with E-I index values for each ego
#' @export
#' @importFrom dplyr left_join group_by summarize mutate filter
#'
#' @examples
#' # Create sample data
#' adj_mat <- matrix(c(0,1,1,0,0,0,1,0,1), nrow=3)
#' rownames(adj_mat) <- colnames(adj_mat) <- c("A", "B", "C")
#' attributes <- data.frame(
#'   node = c("A", "B", "C"),
#'   department = c("HR", "IT", "HR")
#' )
#' calculate_ei_index(adj_mat, attributes, "department")
calculate_ei_index <- function(adj_matrix, node_attributes, attribute_name,
                               directed = TRUE, self_ties = FALSE) {

  # Validate inputs
  if (!is.matrix(adj_matrix)) {
    stop("adj_matrix must be a matrix")
  }

  if (!is.data.frame(node_attributes)) {
    stop("node_attributes must be a data frame")
  }

  if (!(attribute_name %in% colnames(node_attributes))) {
    stop(paste("Attribute", attribute_name, "not found in node_attributes"))
  }

  # Get node names from matrix
  if (is.null(rownames(adj_matrix))) {
    rownames(adj_matrix) <- colnames(adj_matrix) <- 1:nrow(adj_matrix)
  }

  # Prepare node id column name in attributes
  node_id_col <- colnames(node_attributes)[1]

  # Check if node names match between matrix and attributes
  if (!all(rownames(adj_matrix) %in% node_attributes[[node_id_col]])) {
    stop("Node names in adjacency matrix do not match node_attributes")
  }

  # Convert matrix to edge list format
  edges <- data.frame(
    ego = rep(rownames(adj_matrix), each = ncol(adj_matrix)),
    alter = rep(colnames(adj_matrix), times = nrow(adj_matrix)),
    tie = as.vector(adj_matrix)
  )

  # Remove self-ties if specified
  if (!self_ties) {
    edges <- edges[edges$ego != edges$alter, ]
  }

  # Keep only existing ties
  edges <- edges[edges$tie > 0, ]

  # Join with attribute data
  result <- edges %>%
    # Join attributes for ego
    left_join(node_attributes, by = c("ego" = node_id_col)) %>%
    rename(ego_attribute = !!attribute_name) %>%
    # Join attributes for alter
    left_join(node_attributes, by = c("alter" = node_id_col)) %>%
    rename(alter_attribute = !!attribute_name)

  # Calculate E and I values for each ego
  ei_index <- result %>%
    mutate(
      same_group = ego_attribute == alter_attribute,
      # For valued ties, multiply by tie strength
      internal_ties = ifelse(same_group, tie, 0),
      external_ties = ifelse(!same_group, tie, 0)
    ) %>%
    group_by(ego, ego_attribute) %>%
    summarize(
      E = sum(external_ties, na.rm = TRUE),
      I = sum(internal_ties, na.rm = TRUE),
      total_ties = E + I,
      # Calculate E-I index
      EI_index = ifelse(total_ties > 0, (E - I) / total_ties, NA_real_),
      # Additional proportions
      prop_same = I / total_ties,
      prop_different = E / total_ties,
      .groups = "drop"
    )

  return(ei_index)
}

#' Calculate similarity between ego and alters based on continuous attributes
#'
#' @param adj_matrix An adjacency matrix representing network ties
#' @param node_attributes A data frame containing node attributes
#' @param attribute_name The name of the continuous attribute column to analyze
#' @param similarity_method Method to calculate similarity: "absolute_diff", "squared_diff",
#'                         "normalized" (default), or "correlation"
#' @param directed Logical, whether the network is directed. Default is TRUE
#' @param self_ties Logical, whether to include self ties. Default is FALSE
#'
#' @return A data frame with similarity measures between ego and alters
#' @export
#' @importFrom dplyr left_join group_by summarize mutate filter
#' @importFrom stats cor
#'
#' @examples
#' # Create sample data
#' adj_mat <- matrix(c(0,1,1,0,0,0,1,0,1), nrow=3)
#' rownames(adj_mat) <- colnames(adj_mat) <- c("A", "B", "C")
#' attributes <- data.frame(
#'   node = c("A", "B", "C"),
#'   tenure = c(5, 3, 8)
#' )
#' calculate_similarity(adj_mat, attributes, "tenure")
calculate_similarity <- function(adj_matrix, node_attributes, attribute_name,
                                 similarity_method = "normalized",
                                 directed = TRUE, self_ties = FALSE) {

  # Validate inputs
  if (!is.matrix(adj_matrix)) {
    stop("adj_matrix must be a matrix")
  }

  if (!is.data.frame(node_attributes)) {
    stop("node_attributes must be a data frame")
  }

  if (!(attribute_name %in% colnames(node_attributes))) {
    stop(paste("Attribute", attribute_name, "not found in node_attributes"))
  }

  # Check if attribute is numeric
  if (!is.numeric(node_attributes[[attribute_name]])) {
    stop(paste("Attribute", attribute_name, "must be numeric for similarity calculation"))
  }

  # Get node names from matrix
  if (is.null(rownames(adj_matrix))) {
    rownames(adj_matrix) <- colnames(adj_matrix) <- 1:nrow(adj_matrix)
  }

  # Prepare node id column name in attributes
  node_id_col <- colnames(node_attributes)[1]

  # Convert matrix to edge list format (including all possible pairs for correlation)
  all_edges <- data.frame(
    ego = rep(rownames(adj_matrix), each = ncol(adj_matrix)),
    alter = rep(colnames(adj_matrix), times = nrow(adj_matrix)),
    tie = as.vector(adj_matrix)
  )

  # Remove self-ties if specified
  if (!self_ties) {
    all_edges <- all_edges[all_edges$ego != all_edges$alter, ]
  }

  # Join with attribute data
  all_pairs <- all_edges %>%
    # Join attributes for ego
    left_join(node_attributes, by = c("ego" = node_id_col)) %>%
    rename(ego_attribute = !!attribute_name) %>%
    # Join attributes for alter
    left_join(node_attributes, by = c("alter" = node_id_col)) %>%
    rename(alter_attribute = !!attribute_name)

  # Keep only pairs with existing ties for non-correlation methods
  tie_pairs <- all_pairs[all_pairs$tie > 0, ]

  # Calculate attribute range for normalization
  attr_min <- min(node_attributes[[attribute_name]], na.rm = TRUE)
  attr_max <- max(node_attributes[[attribute_name]], na.rm = TRUE)
  attr_range <- attr_max - attr_min

  # Calculate similarity based on selected method
  if (similarity_method == "absolute_diff") {

    result <- tie_pairs %>%
      mutate(
        difference = abs(ego_attribute - alter_attribute),
        similarity = -difference  # Negative difference (higher = more similar)
      ) %>%
      group_by(ego) %>%
      summarize(
        mean_difference = mean(difference, na.rm = TRUE),
        sd_difference = sd(difference, na.rm = TRUE),
        min_difference = min(difference, na.rm = TRUE),
        max_difference = max(difference, na.rm = TRUE),
        mean_similarity = mean(similarity, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (similarity_method == "squared_diff") {

    result <- tie_pairs %>%
      mutate(
        squared_difference = (ego_attribute - alter_attribute)^2,
        similarity = -squared_difference  # Negative squared difference
      ) %>%
      group_by(ego) %>%
      summarize(
        mean_squared_difference = mean(squared_difference, na.rm = TRUE),
        sum_squared_difference = sum(squared_difference, na.rm = TRUE),
        mean_similarity = mean(similarity, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (similarity_method == "normalized") {

    result <- tie_pairs %>%
      mutate(
        # Normalized similarity index (0-1 range, 1 = identical)
        similarity = 1 - (abs(ego_attribute - alter_attribute) / attr_range)
      ) %>%
      group_by(ego) %>%
      summarize(
        mean_similarity = mean(similarity, na.rm = TRUE),
        sd_similarity = sd(similarity, na.rm = TRUE),
        min_similarity = min(similarity, na.rm = TRUE),
        max_similarity = max(similarity, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (similarity_method == "correlation") {

    # Calculate point-biserial correlation for each ego
    result <- data.frame(ego = unique(all_pairs$ego))
    correlations <- numeric(nrow(result))

    for (i in 1:nrow(result)) {
      ego_i <- result$ego[i]

      # Get pairs for this ego
      ego_pairs <- all_pairs[all_pairs$ego == ego_i, ]

      # Calculate absolute difference
      ego_pairs$difference <- abs(ego_pairs$ego_attribute - ego_pairs$alter_attribute)

      # Calculate point-biserial correlation between tie existence and difference
      if (length(unique(ego_pairs$tie)) > 1 && length(unique(ego_pairs$difference)) > 1) {
        correlations[i] <- cor(ego_pairs$tie, ego_pairs$difference,
                               method = "pearson", use = "complete.obs")
      } else {
        correlations[i] <- NA_real_
      }
    }

    result$correlation <- correlations
    result$interpretation <- ifelse(result$correlation < 0,
                                    "Tends to connect with similar alters",
                                    "Tends to connect with different alters")

  } else {
    stop("Invalid similarity_method. Choose from: 'absolute_diff', 'squared_diff', 'normalized', or 'correlation'")
  }

  # Join with ego attribute values
  result <- result %>%
    left_join(
      node_attributes %>%
        select(!!node_id_col, !!attribute_name) %>%
        rename(ego_attribute = !!attribute_name),
      by = c("ego" = node_id_col)
    )

  return(result)
}


#' Analyze categorical composition of alters for each ego in a network
#'
#' @param adj_matrix An adjacency matrix representing network ties
#' @param node_attributes A data frame containing node attributes
#' @param category_col The name of the categorical attribute column to analyze (e.g., "Department")
#' @param directed Logical, whether the network is directed. Default is TRUE
#' @param self_ties Logical, whether to include self ties. Default is FALSE
#' @param values_fill Value to fill for missing combinations in the wide format. Default is 0
#'
#' @return A data frame with categorical composition of alters for each ego in wide format
#' @export
#' @importFrom dplyr left_join group_by summarize mutate rename filter select
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @examples
#' # Create sample adjacency matrix
#' friendship_matrix <- matrix(
#'   c(0,1,1,0,1,0,0,1,1,0,0,1),
#'   nrow = 4, byrow = TRUE
#' )
#' rownames(friendship_matrix) <- colnames(friendship_matrix) <- paste0("P", 1:4)
#'
#' # Create sample attributes
#' attributes <- data.frame(
#'   node = paste0("P", 1:4),
#'   department = c("HR", "IT", "HR", "Finance")
#' )
#'
#' # Analyze department composition
#' analyze_categorical_composition(
#'   friendship_matrix,
#'   attributes,
#'   "department"
#' )
analyze_categorical_composition <- function(adj_matrix, node_attributes, category_col,
                                            directed = TRUE, self_ties = FALSE,
                                            values_fill = 0) {

  # Check if adj_matrix is a matrix
  if (!is.matrix(adj_matrix)) {
    stop("adj_matrix must be a matrix")
  }

  # Check if node_attributes is a data frame
  if (!is.data.frame(node_attributes)) {
    stop("node_attributes must be a data frame")
  }

  # Check if category_col exists in node_attributes
  if (!(category_col %in% colnames(node_attributes))) {
    stop(paste("Column", category_col, "not found in node_attributes"))
  }

  # Determine node ID column (first column)
  node_id_col <- colnames(node_attributes)[1]

  # Ensure rownames and colnames in adj_matrix
  if (is.null(rownames(adj_matrix))) {
    if (nrow(adj_matrix) == nrow(node_attributes)) {
      rownames(adj_matrix) <- node_attributes[[node_id_col]]
      colnames(adj_matrix) <- node_attributes[[node_id_col]]
    } else {
      stop("Cannot determine node names for the adjacency matrix")
    }
  }

  # Check if all nodes in the adjacency matrix exist in node_attributes
  if (!all(rownames(adj_matrix) %in% node_attributes[[node_id_col]])) {
    stop("Some nodes in the adjacency matrix are not found in node_attributes")
  }

  # Convert adjacency matrix to edge list (long format)
  edge_list <- data.frame(
    ego = rep(rownames(adj_matrix), each = ncol(adj_matrix)),
    alter = rep(colnames(adj_matrix), times = nrow(adj_matrix)),
    tie = as.vector(adj_matrix)
  )

  # Remove self-ties if specified
  if (!self_ties) {
    edge_list <- edge_list[edge_list$ego != edge_list$alter, ]
  }

  # Keep only existing ties
  edge_list <- edge_list[edge_list$tie > 0, ]

  # Prepare category data for both ego and alter
  result <- edge_list %>%
    # Join category data for ego
    left_join(node_attributes, by = c("ego" = node_id_col)) %>%
    rename(category_ego = !!category_col) %>%
    # Join category data for alter
    left_join(node_attributes, by = c("alter" = node_id_col)) %>%
    rename(category_alter = !!category_col)

  # Calculate frequencies by ego, ego_category, and alter_category
  composition <- result %>%
    group_by(ego, category_ego, category_alter) %>%
    summarize(
      frequency = sum(tie),
      .groups = "drop_last"
    ) %>%
    # Calculate proportions within each ego
    group_by(ego, category_ego) %>%
    mutate(proportion = frequency / sum(frequency))

  # Calculate IQV (Index of Qualitative Variation) for each ego
  iqv_data <- composition %>%
    group_by(ego) %>%
    summarize(
      num_categories = n_distinct(category_alter),
      sum_prop_squared = sum(proportion^2),
      IQV = ifelse(num_categories > 1,
                   (1 - sum_prop_squared) / (1 - 1/num_categories),
                   0),
      .groups = "drop"
    )

  # Join IQV back to composition data
  composition <- composition %>%
    left_join(iqv_data[, c("ego", "IQV")], by = "ego")

  # Create wide format
  composition_wide <- composition %>%
    pivot_wider(
      id_cols = c(ego, category_ego, IQV),
      names_from = category_alter,
      values_from = c(frequency, proportion),
      values_fill = list(frequency = values_fill, proportion = values_fill)
    )

  return(composition_wide)
}

#' Analyze continuous attribute composition of alters for each ego
#'
#' @param adj_matrix An adjacency matrix representing network ties
#' @param node_attributes A data frame containing node attributes
#' @param attribute_col The name of the continuous attribute column to analyze (e.g., "Tenure")
#' @param category_col Optional. If provided, results will be grouped by this categorical variable
#' @param directed Logical, whether the network is directed. Default is TRUE
#' @param self_ties Logical, whether to include self ties. Default is FALSE
#'
#' @return A data frame with summary statistics of the continuous attribute for each ego's alters
#' @export
#' @importFrom dplyr left_join group_by summarize mutate filter select n
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats mean sd median min max var
#'
#' @examples
#' # Create sample adjacency matrix
#' friendship_matrix <- matrix(
#'   c(0,1,1,0,1,0,0,1,1,0,0,1),
#'   nrow = 4, byrow = TRUE
#' )
#' rownames(friendship_matrix) <- colnames(friendship_matrix) <- paste0("P", 1:4)
#'
#' # Create sample attributes
#' attributes <- data.frame(
#'   node = paste0("P", 1:4),
#'   department = c("HR", "IT", "HR", "Finance"),
#'   tenure = c(5, 3, 8, 2)
#' )
#'
#' # Analyze tenure composition
#' analyze_continuous_composition(
#'   friendship_matrix,
#'   attributes,
#'   "tenure"
#' )
#'
#' # Analyze tenure composition grouped by department
#' analyze_continuous_composition(
#'   friendship_matrix,
#'   attributes,
#'   "tenure",
#'   category_col = "department"
#' )
analyze_continuous_composition <- function(adj_matrix, node_attributes, attribute_col,
                                           category_col = NULL, directed = TRUE,
                                           self_ties = FALSE) {

  # Check if adj_matrix is a matrix
  if (!is.matrix(adj_matrix)) {
    stop("adj_matrix must be a matrix")
  }

  # Check if node_attributes is a data frame
  if (!is.data.frame(node_attributes)) {
    stop("node_attributes must be a data frame")
  }

  # Check if attribute_col exists in node_attributes
  if (!(attribute_col %in% colnames(node_attributes))) {
    stop(paste("Column", attribute_col, "not found in node_attributes"))
  }

  # Check if the attribute is numeric
  if (!is.numeric(node_attributes[[attribute_col]])) {
    stop(paste("Column", attribute_col, "must be numeric"))
  }

  # Check if optional category_col exists (if provided)
  if (!is.null(category_col) && !(category_col %in% colnames(node_attributes))) {
    stop(paste("Column", category_col, "not found in node_attributes"))
  }

  # Determine node ID column (first column)
  node_id_col <- colnames(node_attributes)[1]

  # Ensure rownames and colnames in adj_matrix
  if (is.null(rownames(adj_matrix))) {
    if (nrow(adj_matrix) == nrow(node_attributes)) {
      rownames(adj_matrix) <- node_attributes[[node_id_col]]
      colnames(adj_matrix) <- node_attributes[[node_id_col]]
    } else {
      stop("Cannot determine node names for the adjacency matrix")
    }
  }

  # Check if all nodes in the adjacency matrix exist in node_attributes
  if (!all(rownames(adj_matrix) %in% node_attributes[[node_id_col]])) {
    stop("Some nodes in the adjacency matrix are not found in node_attributes")
  }

  # Convert adjacency matrix to edge list (long format)
  edge_list <- data.frame(
    ego = rep(rownames(adj_matrix), each = ncol(adj_matrix)),
    alter = rep(colnames(adj_matrix), times = nrow(adj_matrix)),
    tie = as.vector(adj_matrix)
  )

  # Remove self-ties if specified
  if (!self_ties) {
    edge_list <- edge_list[edge_list$ego != edge_list$alter, ]
  }

  # Keep only existing ties
  edge_list <- edge_list[edge_list$tie > 0, ]

  # Join attribute data
  result <- edge_list

  # Add ego's attributes
  selected_ego_cols <- c(node_id_col, attribute_col)
  if (!is.null(category_col)) {
    selected_ego_cols <- c(selected_ego_cols, category_col)
  }

  result <- result %>%
    # Join ego attributes
    left_join(node_attributes[, selected_ego_cols], by = c("ego" = node_id_col))

  # Rename ego columns
  if (!is.null(category_col)) {
    result <- result %>%
      rename(
        ego_attribute = !!attribute_col,
        ego_category = !!category_col
      )
  } else {
    result <- result %>%
      rename(ego_attribute = !!attribute_col)
  }

  # Join alter attributes
  result <- result %>%
    left_join(node_attributes[, c(node_id_col, attribute_col)], by = c("alter" = node_id_col)) %>%
    rename(alter_attribute = !!attribute_col)

  # Prepare grouping variables for summarization
  if (!is.null(category_col)) {
    grouping_vars <- c("ego", "ego_category")
  } else {
    grouping_vars <- "ego"
  }

  # Calculate summary statistics
  summary_stats <- result %>%
    group_by_at(grouping_vars) %>%
    summarize(
      ego_attribute_value = first(ego_attribute),
      num_alters = n(),
      sum_attribute = sum(alter_attribute * tie, na.rm = TRUE),
      mean_attribute = mean(alter_attribute, na.rm = TRUE),
      median_attribute = median(alter_attribute, na.rm = TRUE),
      sd_attribute = sd(alter_attribute, na.rm = TRUE),
      var_attribute = var(alter_attribute, na.rm = TRUE),
      min_attribute = min(alter_attribute, na.rm = TRUE),
      max_attribute = max(alter_attribute, na.rm = TRUE),
      range_attribute = max_attribute - min_attribute,
      weighted_mean = sum(alter_attribute * tie, na.rm = TRUE) / sum(tie, na.rm = TRUE),
      abs_diff_from_ego = mean(abs(alter_attribute - ego_attribute_value), na.rm = TRUE),
      .groups = "drop"
    )

  return(summary_stats)
}
