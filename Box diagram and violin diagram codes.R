library(pcutils)
library(ggstatsplot)
library(PMCMRplus)
library(ggplot2)
library(cowplot)

set.seed(123)

# Contains two columns for the gene name and group
data <- read.table("expression_matrix.txt", sep = "\t", check.names = F, stringsAsFactors = F, header = T)

{
  # Adjust according to your actual data structure
  plist = list()
  for (i in 1:9) {
    plist[[i]] = group_box(tab = data[c("ACSF2")], group = "group", metadata = data,
                           mode = i, p_value1 = "t.test") + # Compare using methods such as t.test, wilcox.test, anova, or kruskal.test
      ggtitle(paste0("mode", i)) +
      theme_classic() +
      theme(legend.position = "none")
  }
  
  # Draw box plots
  plot_grid(plotlist = plist, ncol = 3)
}


# Display a specific mode separately
plot(plist[[4]])