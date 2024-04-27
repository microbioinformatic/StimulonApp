#Read again export mathes GO#####

matches_df <- read_csv("matches_df.csv")
subset_GO_df_full_BP <- subset(matches_df, matches_df$`GO.ASPECT` == "C")
annotations = subset_GO_df_full_BP[,c(8,3,4)]
colnames(annotations) <- c("SYMBOL", "GO_TERM", "GO_NAME")

set1 <- colnames(Salmonella_expr_As)

all_genes <- colnames(Salmonella_expr)

genes_of_interest <- set1
# # Create a named vector with all genes and mark genes of interest
allGenes <- factor(ifelse(all_genes %in% genes_of_interest, "genes_of_interest", "background_genes"))

# # Assign names to the vector
names(allGenes) <- all_genes

# Remove duplicate rows based on SYMBOL and GO_TERM columns
unique_annotations <- annotations[!duplicated(annotations[c("SYMBOL", "GO_TERM")], fromLast = TRUE), ]

colnames(unique_annotations) <- c("SYMBOL", "GO_TERM", "GO_NAME")

# Example data (replace this with your actual data)
#set1 <- rownames(mynoiseq2.deg1)
genes_of_interest <- set1

# Create a named vector with all genes and mark genes of interest
allGenes <- factor(ifelse(all_genes %in% genes_of_interest, "genes_of_interest", "background_genes"))
names(allGenes) <- all_genes

# Remove duplicate rows based on SYMBOL and GO_TERM columns
unique_annotations <- annotations[!duplicated(annotations[c("SYMBOL", "GO_TERM")], fromLast = TRUE), ]
colnames(unique_annotations) <- c("SYMBOL", "GO_TERM", "GO_NAME")

# Initialize a list to store enrichment results
enrichment_results <- list()
p_values <- c()  # Store p-values for multiple testing correction

# Iterate over each unique GO term
unique_go_terms <- unique(unique_annotations$GO_TERM)
for (go_term in unique_go_terms) {
  # Subset annotations for the current GO term
  go_annotations <- unique_annotations[unique_annotations$GO_TERM == go_term, ]
  
  go_count_genes_of_interest <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL] == "genes_of_interest")
  go_count_background <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL]  == "background_genes")
  
  total_genes_of_interest <- sum(allGenes[names(allGenes)] == "genes_of_interest")
  total_background_genes <- sum(allGenes[names(allGenes)]  == "background_genes")
  
  # Perform Fisher's exact test if there are sufficient counts for comparison
  if (go_count_genes_of_interest >= 3 && go_count_background >= 0) {
    fisher_result <- fisher.test(matrix(c(go_count_genes_of_interest, go_count_background, 
                                          total_genes_of_interest - go_count_genes_of_interest,
                                          total_background_genes - go_count_background), 
                                        nrow = 2), alternative = "greater")
    
    # Extract p-value
    p_value <- fisher_result$p.value
    p_values <- c(p_values, p_value)  # Store p-value
    
    # Store enrichment result in the list
    enrichment_results[[length(enrichment_results) + 1]] <- list(GO_TERM = go_term, GO_NAME = unique(go_annotations$GO_NAME))
    
    # Check if enrichment is significant (e.g., p-value < 0.05)
    if (p_value < 0.05) {
      # Print out the GO term and GO name
      cat("Significant enrichment for GO term:", go_term, "\n")
      cat("GO name:", unique(go_annotations$GO_NAME), "\n\n")
    }
  }
}

# Adjust p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Filter significant enrichments based on adjusted p-values
significant_enrichments <- enrichment_results[adjusted_p_values < 0.05]

# If you want to have a combined list of significant enrichments
enrichment_results_list <- do.call(rbind, significant_enrichments)

enrichment_df <- as.data.frame(enrichment_results_list)

# Create a function to calculate proportions for each GO term
calculate_proportions <- function(go_term) {
  # Subset annotations for the current GO term
  go_annotations <- unique_annotations[unique_annotations$GO_TERM == go_term, ]
  
  # Count occurrences of the current GO term in genes_of_interest and background_genes
  go_count_genes_of_interest <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL] == "genes_of_interest")
  go_count_background <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL]  == "background_genes")
  
  # Calculate proportions
  proportions <- data.frame(GO_TERM = go_term,
                            GO_NAME = unique(go_annotations$GO_NAME),
                            Group = c("Genes_of_Interest", "Background_Genes"),
                            Proportion = c(go_count_genes_of_interest / total_genes_of_interest,
                                           go_count_background / total_background_genes))
  
  return(proportions)
}

# Create a dataframe to store combined plot data
combined_plot_data <- data.frame(GO_TERM = character(),
                                 GO_NAME = character(),
                                 Group = character(),
                                 Proportion = numeric())

# Create a bar plot for each significant GO term
for (go_term in unique(enrichment_df$GO_TERM)) {
  plot_data <- calculate_proportions(go_term)
  combined_plot_data <- bind_rows(combined_plot_data, plot_data)
}

# Plot stacked bars for all significant GO terms
p <- ggplot(combined_plot_data, aes(x = Proportion, y = GO_NAME, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Enrichment for Significant GO Terms",
       x = "Proportion", y = "GO Name") +
  theme_minimal()

# Print the plot
print(p)

#####More ways#####

set1= colnames(Salmonella_expr_As)

all_genes <- colnames(Salmonella_expr)
genes_of_interest <- set1
# # Create a named vector with all genes and mark genes of interest
allGenes <- factor(ifelse(all_genes %in% genes_of_interest, "genes_of_interest", "background_genes"))
# # Assign names to the vector
names(allGenes) <- all_genes

# Remove duplicate rows based on SYMBOL and GO_TERM columns
unique_annotations <- annotations[!duplicated(annotations[c("SYMBOL", "GO_TERM")], fromLast = TRUE), ]

colnames(unique_annotations) <- c("SYMBOL", "GO_TERM", "GO_NAME")
# Initialize a list to store enrichment results
enrichment_results <- list()
p_values <- c()  # Store p-values for multiple testing correction

# Iterate over each unique GO term
unique_go_terms <- unique(unique_annotations$GO_TERM)
for (go_term in unique_go_terms) {
  # Subset annotations for the current GO term
  go_annotations <- unique_annotations[unique_annotations$GO_TERM == go_term, ]
  
  go_count_genes_of_interest <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL] == "genes_of_interest")
  go_count_background <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL]  == "background_genes")
  
  total_genes_of_interest <- sum(allGenes[names(allGenes)] == "genes_of_interest")
  total_background_genes <- sum(allGenes[names(allGenes)]  == "background_genes")
  
  # Perform Fisher's exact test if there are sufficient counts for comparison
  if (go_count_genes_of_interest >= 3 && go_count_background >= 0) {
    fisher_result <- fisher.test(matrix(c(go_count_genes_of_interest, go_count_background, 
                                          total_genes_of_interest - go_count_genes_of_interest,
                                          total_background_genes - go_count_background), 
                                        nrow = 2), alternative = "greater")
    
    # Extract p-value
    p_value <- fisher_result$p.value
    p_values <- c(p_values, p_value)  # Store p-value
    
    # Store enrichment result in the list
    enrichment_results[[length(enrichment_results) + 1]] <- list(GO_TERM = go_term, GO_NAME = unique(go_annotations$GO_NAME))
    
    # Check if enrichment is significant (e.g., p-value < 0.05)
    if (p_value < 0.05) {
      # Print out the GO term and GO name
      cat("Significant enrichment for GO term:", go_term, "\n")
      cat("GO name:", unique(go_annotations$GO_NAME), "\n\n")
    }
  }
}

# Adjust p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Filter significant enrichments based on adjusted p-values
significant_enrichments <- enrichment_results[adjusted_p_values < 0.05]

# If you want to have a combined list of significant enrichments
enrichment_results_list <- do.call(rbind, significant_enrichments)

enrichment_df <- as.data.frame(enrichment_results_list)

# Create a function to calculate proportions for each GO term
calculate_proportions <- function(go_term) {
  # Subset annotations for the current GO term
  go_annotations <- unique_annotations[unique_annotations$GO_TERM == go_term, ]
  
  # Count occurrences of the current GO term in genes_of_interest and background_genes
  go_count_genes_of_interest <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL] == "genes_of_interest")
  go_count_background <- sum(allGenes[names(allGenes) %in% go_annotations$SYMBOL]  == "background_genes")
  
  # Calculate proportions
  proportions <- data.frame(GO_TERM = go_term,
                            GO_NAME = unique(go_annotations$GO_NAME),
                            Group = c("Genes_of_Interest", "Background_Genes"),
                            Proportion = c(go_count_genes_of_interest / total_genes_of_interest,
                                           go_count_background / total_background_genes))
  
  return(proportions)
}

# Create a dataframe to store combined plot data
combined_plot_data <- data.frame(GO_TERM = character(),
                                 GO_NAME = character(),
                                 Group = character(),
                                 Proportion = numeric())

# Create a bar plot for each significant GO term
for (go_term in unique(enrichment_df$GO_TERM)) {
  plot_data <- calculate_proportions(go_term)
  combined_plot_data <- bind_rows(combined_plot_data, plot_data)
}


# Calculate fold change within each GO term
combined_plot_data <- combined_plot_data %>%
  group_by(GO_NAME) %>%
  mutate(
    FoldChange = ifelse(
      Proportion[Group == "Background_Genes"] == 0,
      log2(Proportion[Group == "Genes_of_Interest"]),
      log2(Proportion[Group == "Genes_of_Interest"] / (Proportion[Group == "Background_Genes"] + 1e-10))
    )
  )

# Plot grouped bar plot with fold change
p <- ggplot(combined_plot_data, aes(x = reorder(GO_NAME, Proportion), y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_text(data = filter(combined_plot_data, Group == "Genes_of_Interest"), aes(label = round(FoldChange, 2)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3, color = "black") +
  labs(title = "Enrichment for Significant GO Terms",
       x = "GO Name", y = "Proportion", fill = "Group") +
  scale_fill_manual(values = c("Genes_of_Interest" = "blue", "Background_Genes" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  coord_flip()  # Flip the coordinates to make bars horizontal

# Print the plot
print(p)

# Plot scatter plot with fold change as x-axis and GO terms as y-axis
p <- ggplot(combined_plot_data, aes(x = FoldChange, y = GO_NAME)) +
  geom_point(aes(size = Proportion, fill = Group), shape = 21) +
  scale_size_area(max_size = 20) +  # Adjust maximum circle size
  scale_fill_manual(values = c("Genes_of_Interest" = "blue", "Background_Genes" = "red")) +
  labs(title = "Enrichment for Significant GO Terms",
       x = "Fold Change", y = "GO Name", size = "Proportion", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "right", axis.text.y = element_text(size = 8))  # Adjust font size for GO terms

# Print the plot
print(p)

#######other#####

# Calculate fold change within each GO term
combined_plot_data <- combined_plot_data %>%
  group_by(GO_NAME) %>%
  mutate(
    FoldChange = ifelse(
      Proportion[Group == "Background_Genes"] == 0,
      ifelse(Proportion[Group == "Genes_of_Interest"] == 0, 0, Inf),  # Set to 0 if both are 0, else Inf
      log2(Proportion[Group == "Genes_of_Interest"] / Proportion[Group == "Background_Genes"])
    )
  )

# Replace Inf with a large value
combined_plot_data$FoldChange[combined_plot_data$FoldChange == Inf] <- 10

# Plot scatter plot with fold change as x-axis and GO terms as y-axis
p <- ggplot(combined_plot_data, aes(x = FoldChange, y = GO_NAME)) +
  geom_point(aes(size = Proportion, fill = Group), shape = 21) +
  scale_size_area(max_size = 20) +  # Adjust maximum circle size
  scale_fill_manual(values = c("Genes_of_Interest" = "blue", "Background_Genes" = "red")) +
  labs(title = "Enrichment for Significant GO Terms",
       x = "Fold Change", y = "GO Name", size = "Proportion", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "right", axis.text.y = element_text(size = 8))  # Adjust font size for GO terms

# Print the plot
print(p)

#####other2####

# Reorder GO_NAME by FoldChange
combined_plot_data$GO_NAME <- factor(combined_plot_data$GO_NAME, levels = unique(combined_plot_data$GO_NAME[order(combined_plot_data$FoldChange)]))

# Plot scatter plot with fold change as x-axis and GO terms as y-axis, sorted by fold change
p <- ggplot(combined_plot_data, aes(x = FoldChange, y = GO_NAME)) +
  geom_point(aes(size = Proportion, fill = Group), shape = 21) +
  scale_size_area(max_size = 20) +  # Adjust maximum circle size
  scale_fill_manual(values = c("Genes_of_Interest" = "blue", "Background_Genes" = "red")) +
  labs(title = "Enrichment for Significant GO Terms",
       x = "Fold Change", y = "GO Name", size = "Proportion", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "right", axis.text.y = element_text(size = 8))  # Adjust font size for GO terms

# Print the plot
print(p)


