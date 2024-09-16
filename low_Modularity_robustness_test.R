
#LIBRARIES
library(igraph)
library(netUtils)
library(networkD3)
library(ggplot2)


#FUNCTIONS
## TO GENERATE NULL MODEL FOR GIVEN GRAPH
random_self <- function(graph, z=NULL,verbose=FALSE)
{
  print("Randomizing the graph edges.\n")
  if (is.null(z)) {
    
    z <- igraph::gsize(graph) ## number of edges
    graphRandom <- igraph::rewire(graph, 
                                  with=igraph::keeping_degseq(loops=FALSE, niter=z))
    if (!is.null(E(graph)$weight)){
      #rewiring for z all the edges
      print("works")
      original_weights <- E(graph)$weight
      mean_weight <- mean(original_weights)
      sd_weight <- sd(original_weights)
      
      # Generate new random weights with the same statistical properties
      new_weights <- rnorm(length(original_weights), mean = mean_weight, sd = sd_weight)
      
      # Ensure non-negative weights (if necessary)
      new_weights <- round(pmax(new_weights, 0))
      
      E(graphRandom)$weight <- new_weights
    }
    return(graphRandom)
    
  }else{
    graphRandom <- igraph::rewire(graph, 
                                  with=igraph::keeping_degseq(loops=FALSE, niter=z))
    #rewiring for z all the edges
    original_weights <- E(graph)$weight
    mean_weight <- mean(original_weights)
    sd_weight <- sd(original_weights)
    
    # Generate new random weights with the same statistical properties
    new_weights <- rnorm(length(original_weights), mean = mean_weight, sd = sd_weight)
    
    # Ensure non-negative weights (if necessary)
    new_weights <- round(pmax(new_weights, 0))
    
    E(graphRandom)$weight <- new_weights
    return(graphRandom)
  }
}

## SIMULATED RANDOM NETWORKS OF VARIABLE PARAMETERS
LFR<- function(n, tau1, tau2, mu,min_community,max_community,average_degree, max_degree)
{
  lfr_graph <- sample_lfr(n, tau1, tau2, mu, 
                          min_community = min_community, 
                          max_community = max_community,
                          average_degree = average_degree, 
                          max_degree = max_degree)
  mean_weight <- 10
  sd_weight <- 2
  E(lfr_graph)$weight <- abs(rnorm(ecount(lfr_graph), mean = mean_weight, sd = sd_weight))
  
  
  communities1 <- cluster_louvain(lfr_graph)
  membership_vec1 <- membership(communities1)
  
  
  modularity_value <- modularity(communities1)
  print(modularity_value)
  list(graph = lfr_graph, modularity = modularity_value)
}

## CREATING LOW MODULARITY
Graph_modular_lOW<-LFR(n=500, tau1=6, tau2=3, mu=0.6,min_community=50,max_community=100,average_degree=55, max_degree=80)
graph<-Graph_modular_lOW$graph


graph<-random_self(graph)
rand<-random_self(graph)

par(mfrow = c(1,1))
plot(graph,
     layout = layout_with_graphopt(graph),
     vertex.size = 0,             # Size of the vertices
     vertex.label.cex =0.5,      # Size of the vertex labels
     vertex.label.color ="red",
     edge.width = 1,  # Width of the edges based on weight
     edge.label = round(E(graph)$weight, 1),  # Edge labels showing weights
     edge.label.width=0.2,
     edge.label.cex = 0.7,        # Size of the edge labels
     edge.label.color = "black",   # Color of the edge labels
     edge.label.dist = 0.5,       # Distance of the labels from the edges
     main = "Graph (Modularity=0.27)")
plot(rand,
     layout = layout_with_graphopt(rand),
     vertex.size = 0,             # Size of the vertices
     vertex.label.cex =0.5,      # Size of the vertex labels
     vertex.label.color ="red",
     edge.width = 1,  # Width of the edges based on weight
     edge.label = round(E(rand)$weight, 1),  # Edge labels showing weights
     edge.label.width=0.2,
     edge.label.cex = 0.7,        # Size of the edge labels
     edge.label.color = "black",   # Color of the edge labels
     edge.label.dist = 0.5,       # Distance of the labels from the edges
     main = "Null Model ")

## ORIGINAL COMMUNITIES
og_com<-igraph::cluster_louvain(graph,E(graph)$weight)
rand_com<-igraph::cluster_louvain(rand,E(rand)$weight)

members_og <- membership(og_com)
members_rand <- membership(rand_com)

par(mfrow = c(1, 1)) # Set up the plotting area for two plots

# Plot for original graph
plot(graph, vertex.color=members_og, main="Graph Community (Modularity:0.27)", 
     vertex.label=NA, vertex.size=5, edge.arrow.size=0.5, 
     layout=layout_with_fr(graph))

# Plot for random graph
plot(rand, vertex.color=members_rand, main="Random Graph Community", 
     vertex.label=NA, vertex.size=5, edge.arrow.size=0.5, 
     layout=layout_with_fr(rand))

m1 <- list()
m2 <- list()

## INTRODUCING PERURBATION AND CALCULATION STABILITY MEASURES
for (i in seq(0, 80, by = 5)) {
  temp_m1 <- numeric(10)  
  temp_m2 <- numeric(10)
  
  counter=0
  while (counter < 10) {
    if (i==0){
      measure1 <- 0
      measure2 <- 0
    }else{
      
      
      og_perturb<-random_self(graph,round(i*igraph::gsize(graph)/100, 0))
      
      rand_perturb<-random_self(rand,round(i*igraph::gsize(rand)/100, 0))
      
      
      
      og_perturb_com<-igraph::cluster_louvain(og_perturb,E(og_perturb)$weight)
      rand_perturb_com<-igraph::cluster_louvain(rand_perturb,E(rand_perturb)$weight)
      
      members_perturb_og <- membership(og_perturb_com)
      members_perturb_rand <- membership(rand_perturb_com)
      
      
      measure1 <-  igraph::compare(members_og, members_perturb_og,method="vi")/log2(igraph::vcount(graph))
      measure2 <-  igraph::compare(members_rand, members_perturb_rand,method="vi")/log2(igraph::vcount(graph))
    }
    # Store measures in temporary vectors
    temp_m1[counter + 1] <- measure1
    temp_m2[counter + 1] <- measure2
    
    counter=counter+1
  }
  # Store temporary vectors in the list of matrices
  m1[[as.character(i)]] <- temp_m1
  m2[[as.character(i)]] <- temp_m2
  
}
# Convert lists to data frames for easier viewing
m1_df <- do.call(cbind, m1)
m2_df <- do.call(cbind, m2)

print(m1_df)
print(m2_df)

saveRDS(m1_df, file = "low_m1_df.rds")
saveRDS(m2_df, file = "low_m2_df.rds")

# Calculate the means of each column
m1_means <- colMeans(m1_df)
m2_means <- colMeans(m2_df)


# Convert the means to a data frame for easier plotting
means_df <- data.frame(
  i_values = as.numeric(colnames(m1_df)),
  m1_means = m1_means,
  m2_means = m2_means
  
)
means_df

# Plot the means
ggplot(means_df, aes(x = i_values)) +
  geom_line(aes(y = m1_means, color = "m1")) +
  geom_line(aes(y = m2_means, color = "m2")) +
  
  labs(
    title = "Stability Measure for Graph(Modularity=0.27)",
    x = "Perturbations(%)",
    y = "Mean Variational Information Score "
  ) +
  scale_color_manual(
    name = "Measures",
    values = c("m1" = "blue", "m2" = "green"),
    labels = c("Original Graph", "Null Model")
  ) +
  theme_minimal()

#FDA TEST
model1<-m1_df
model2<-m2_df
graph<-graph
muParam=0
orderParam=4
nKnots=7
BParam=10000
isPaired=TRUE
verbose=TRUE
modeled1 <- as.matrix(model1)
modeled2 <- as.matrix(model2)   
## Two populations Interval Testing Procedure with B-spline basis
ITPresult <- fdatest::ITP2bspline(modeled1, modeled2, mu=muParam, 
                                  order=orderParam, nknots=nKnots, 
                                  B=BParam, paired=isPaired)

if(verbose) cat("Computing Interval testing procedure.\n")

graph <- graph

legend <- c("Low Modularity Graph","null model")


object <- ITPresult
str(object)
#Functional Data plot
J <- dim(object$data.eval)[2]
xmin <- 0
xmax <- 0.8
Abscissa <- seq(xmin,xmax,len=J)

model1 <- cbind(as.numeric(as.vector(t(object$data.eval[1:10,]))))
model2 <- cbind(as.numeric(as.vector(t(object$data.eval[11:20,]))))

measures <- rbind(model1, model2)
model <- c(rep(legend[1],each=10000),rep(legend[2],each=10000))
percPert <- as.numeric(rep(Abscissa, times = 10))
s <- c(rep(1:10,each=1000),rep(11:20,each=1000))
dataFrame <- data.frame(measures,model,s,percPert)
plot1 <- ggplot2::ggplot(dataFrame, ggplot2::aes(x=as.numeric(percPert),
                                                 y=as.numeric(measures), color= model, group=s)) +
  ggplot2::geom_line() +
  ggplot2::xlab("Percentage of perturbation") +
  ggplot2::ylab("Measure")+
  ggplot2::ggtitle("Functional Data Analysis")+
  ggplot2::scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7,0.8))


plot1

#P value plot
p <- length(object$pval)
xmin <- 0
xmax <- 0.8
abscissa.pval <- rep(seq(xmin,xmax,len=p),time=2)
pvalue <- c(object$pval,object$corrected.pval)
type <- c(rep("pvalue",p),rep("pvalue.adj",p))
PdataFrame <- data.frame(cbind(abscissa.pval,pvalue,type))

plot2 <- ggplot2::ggplot(PdataFrame, ggplot2::aes(x = as.numeric(abscissa.pval),
                                                  y = as.numeric(pvalue), 
                                                  color = type)) +
  ggplot2::geom_point(size = 3, shape = 16, alpha = 0.7) +  # Increase size, set shape, and adjust transparency
  ggplot2::xlab("Percentage of perturbation") +
  ggplot2::ylab("p_value") +
  ggplot2::ggtitle("P-values") +
  ggplot2::scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) +
  ggplot2::geom_hline(yintercept = 0.05, color = "red") +
  ggplot2::lims(y = c(0, 1))

# Display the plot
print(plot2)

plot <- gridExtra::grid.arrange(plot1,plot2, ncol=2)
print(plot)

adj.pvalue <- object$corrected.pval
pvalue <- object$pval
output <- list(adj.pvalue=adj.pvalue,
               pvalues=pvalue)
output
