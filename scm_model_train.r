library(bnlearn)
library(Rgraphviz)
library(mygene)
library(dplyr)
library(purrr)
load("feature_gene.rda")
load("sig_gene.rda")
load("cancer_exp_mt.rda")
load("cluster_assign_tumor.rda")
load("cancer_sample_df.rda")
gene_ls<-as.vector(unlist(gene.df$SYMBOL))
Protein_IDs <- queryMany(gene_ls, scopes = "symbol", 
                              fields = c("name", "uniprot",  "ensemblgene"), 
                              species = "human", as_dataframe = "True")

ppi_inter<-read.table("string_interactions_short1129.tsv",header=F,sep="\t")            
gene_list<-sig_g
ppi_list<-ppi_inter[,c(1,2)]
colnames(ppi_list)<-c("from","to")
ppi_list%>%dplyr::filter(from%in%gene_list)%>%dplyr::filter(to%in%gene_list)->ppi_list_filter
TF_target_data<-read.table("trrust_rawdata.human.tsv",header=T,sep="\t")
#gene_list<-gene.df$SYMBOL
TF_target_data%>%dplyr::filter(TF%in%gene_list)%>%dplyr::filter(target%in%gene_list)->TF_target_relation
TF_target_relation[,c(1,2)]->tf_target_list
colnames(tf_target_list)<-c("from","to")
white_list<-rbind(ppi_list_filter,tf_target_list)
sample_cluster<-cluster_assign



if (!all(colnames(white_list) == c("from", "to"))) {
  colnames(white_list) <- c("from", "to")
}
set.seed(123)
for(sub in unique(sample_cluster)){
    train_data<-cbind(t(cancer_exp_mt[gene_list,names(which(sample_cluster==sub))]),sample_cluster)
    colnames(train_data)[ncol(train_data)] <- "Label"
    learned_dag <- hc(train_data, whitelist = white_list)
    fitted_bn <- bn.fit(learned_dag, train_data)
    graphviz.plot(learned_dag, layout = "fdp")
    save(fitted_bn,file=(paste(sub,"_bnn_res.rda",sep="")))
}


# Load necessary libraries


# Set random seed for reproducibility
set.seed(123)

# Prior perturbation analysis function
prior_perturbation_analysis <- function(train_data, initial_white_list, n_perturbations = 500, perturbation_rate = 0.2) {
  
  # List to store all results
  results <- list()
  
  # Get information about initial white list
  initial_edges <- nrow(initial_white_list)
  nodes <- unique(c(initial_white_list$from, initial_white_list$to, colnames(train_data)))
  nodes <- nodes[nodes %in% colnames(train_data)] # Ensure nodes exist in data
  
  cat("Starting prior perturbation analysis...\n")
  cat("Initial white list edges:", initial_edges, "\n")
  cat("Perturbation rate:", perturbation_rate, "\n")
  cat("Number of perturbations:", n_perturbations, "\n")
  
  # Generate all possible edges (excluding self-loops)
  all_possible_edges <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE) %>%
    filter(from != to) %>%
    distinct()
  
  # Convert edges to string sets for easier comparison
  initial_edges_set <- paste(initial_white_list$from, initial_white_list$to, sep = "|")
  all_possible_edges_set <- paste(all_possible_edges$from, all_possible_edges$to, sep = "|")
  
  # Find edges not in initial white list
  non_white_edges <- all_possible_edges[!all_possible_edges_set %in% initial_edges_set, ]
  
  # Perform perturbation analysis
  for(i in 1:n_perturbations) {
    cat("Processing perturbation", i, "/", n_perturbations, "\n")
    
    # Create perturbed white list by copying initial list
    perturbed_white_list <- initial_white_list
    
    # Calculate number of edges to add/remove
    n_perturb <- max(1, round(initial_edges * perturbation_rate))
    
    # Randomly remove edges (if there are edges to remove)
    if(nrow(perturbed_white_list) > n_perturb) {
      edges_to_remove <- sample(1:nrow(perturbed_white_list), n_perturb)
      perturbed_white_list <- perturbed_white_list[-edges_to_remove, ]
    }
    
    # Randomly add edges from non-white list edges (if available)
    if(nrow(non_white_edges) > 0 && n_perturb > 0) {
      edges_to_add <- sample(1:nrow(non_white_edges), 
                            min(n_perturb, nrow(non_white_edges)))
      
      new_edges <- non_white_edges[edges_to_add, ]
      colnames(new_edges) <- c("from", "to")
      perturbed_white_list <- rbind(perturbed_white_list, new_edges)
    }
    
    # Remove any duplicate edges that might have been created
    perturbed_white_list <- perturbed_white_list %>%
      distinct(from, to)
    
    # Learn SCM structure with perturbed white list
    tryCatch({
      # Determine score function based on data type
      data_types <- sapply(train_data, class)
      score_func <- ifelse("factor" %in% data_types, "bic-cg", "bic-g")
      
      # Learn structure with perturbed constraints
      learned_dag <- hc(train_data, 
                       whitelist = perturbed_white_list,
                       score = score_func,
                       restart = 5,    # Increase restarts for stability
                       perturb = 3)    # Increase perturbations for better exploration
      
      # Learn parameters
      fitted_bn <- bn.fit(learned_dag, train_data)
      
      # Extract edge information from learned structure
      learned_arcs <- arcs(learned_dag)
      
      # Store results
      results[[i]] <- list(
        perturbation_id = i,
        perturbed_white_list = perturbed_white_list,
        learned_dag = learned_dag,
        fitted_bn = fitted_bn,
        learned_arcs = learned_arcs,
        n_edges_learned = nrow(learned_arcs),
        success = TRUE
      )
      
    }, error = function(e) {
      # Store error information if learning fails
      cat("Error in perturbation", i, ":", e$message, "\n")
      results[[i]] <- list(
        perturbation_id = i,
        perturbed_white_list = perturbed_white_list,
        error = e$message,
        success = FALSE
      )
    })
  }
  
  return(results)
}

# Function to analyze robustness of learned edges
analyze_edge_robustness <- function(perturbation_results) {
  
  # Extract successful results only
  successful_results <- perturbation_results[sapply(perturbation_results, function(x) x$success)]
  
  if(length(successful_results) == 0) {
    stop("No successful perturbations found")
  }
  
  # Collect all learned edges across perturbations
  all_edges <- lapply(successful_results, function(result) {
    if(!is.null(result$learned_arcs)) {
      edges_df <- as.data.frame(result$learned_arcs)
      edges_df$edge_id <- paste(edges_df$from, edges_df$to, sep = "|")
      return(edges_df)
    }
    return(NULL)
  })
  
  all_edges_df <- bind_rows(all_edges)
  
  # Calculate edge frequency across perturbations
  edge_frequency <- all_edges_df %>%
    group_by(from, to) %>%
    summarise(
      frequency = n() / length(successful_results),
      .groups = 'drop'
    ) %>%
    arrange(desc(frequency))
  
  return(edge_frequency)
}

# Function to visualize robust edges
visualize_robust_network <- function(edge_frequency, frequency_threshold = 0.7) {
  
  # Filter edges above frequency threshold
  robust_edges <- edge_frequency %>%
    filter(frequency >= frequency_threshold)
  
  if(nrow(robust_edges) == 0) {
    cat("No edges found above frequency threshold:", frequency_threshold, "\n")
    return(NULL)
  }
  
  # Create a new BN object with robust edges
  robust_arcs <- as.matrix(robust_edges[, c("from", "to")])
  robust_dag <- empty.graph(nodes = unique(c(robust_arcs)))
  arcs(robust_dag) <- robust_arcs
  
  # Plot the robust network
  graphviz.plot(robust_dag, 
               layout = "fdp",
               main = paste("Robust Causal Network (Frequency >=", frequency_threshold, ")"))
  
  return(robust_dag)
}

# Main execution function
run_robust_scm_analysis <- function(train_data, initial_white_list, 
                                   n_perturbations = 500, perturbation_rate = 0.2,
                                   frequency_threshold = 0.7) {
  
  # Data preparation
  cat("=== Data Preparation ===\n")
  
  # Ensure the last column is named "Label"
  colnames(train_data)[ncol(train_data)] <- "Label"
  
  # Ensure white list has correct column names
  if (!all(colnames(initial_white_list) == c("from", "to"))) {
    colnames(initial_white_list) <- c("from", "to")
  }
  
  # Convert factors if needed
  train_data[] <- lapply(train_data, function(x) {
    if(is.character(x)) as.factor(x) else x
  })
  
  # Run perturbation analysis
  cat("\n=== Running Perturbation Analysis ===\n")
  perturbation_results <- prior_perturbation_analysis(
    train_data = train_data,
    initial_white_list = initial_white_list,
    n_perturbations = n_perturbations,
    perturbation_rate = perturbation_rate
  )
  
  # Analyze edge robustness
  cat("\n=== Analyzing Edge Robustness ===\n")
  edge_frequency <- analyze_edge_robustness(perturbation_results)
  
  # Print summary statistics
  cat("\n=== Robustness Analysis Summary ===\n")
  cat("Total perturbations:", length(perturbation_results), "\n")
  cat("Successful perturbations:", sum(sapply(perturbation_results, function(x) x$success)), "\n")
  cat("Unique edges learned:", nrow(edge_frequency), "\n")
  cat("Edges with frequency >= 0.9:", sum(edge_frequency$frequency >= 0.9), "\n")
  cat("Edges with frequency >= 0.7:", sum(edge_frequency$frequency >= 0.7), "\n")
  cat("Edges with frequency >= 0.5:", sum(edge_frequency$frequency >= 0.5), "\n")
  
  # Display most robust edges
  cat("\nTop 10 most robust edges:\n")
  print(head(edge_frequency, 10))
  
  # Visualize robust network
  cat("\n=== Visualizing Robust Network ===\n")
  robust_dag <- visualize_robust_network(edge_frequency, frequency_threshold)
  
  # Return comprehensive results
  return(list(
    perturbation_results = perturbation_results,
    edge_frequency = edge_frequency,
    robust_network = robust_dag,
    summary = list(
      total_perturbations = length(perturbation_results),
      successful_perturbations = sum(sapply(perturbation_results, function(x) x$success)),
      unique_edges_learned = nrow(edge_frequency)
    )
  ))
}


results <- run_robust_scm_analysis(train_data, white_list, 
                                   n_perturbations = 500, 
                                   perturbation_rate = 0.2)

