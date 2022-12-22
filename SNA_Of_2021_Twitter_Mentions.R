######Package load########

# Loads packages using the groundhog package. Groundhog enables reproducible     
# analysis by recording the date the packages were used. It then downloads those 
# exact package version upon later analysis, ensuring scripts run as intended.

library("groundhog")
pkgs <- c("tidyverse","igraph", "splitstackshape", "ggraph", "kableExtra", "gplots", "skimr", "DescTools", "FSA", "qgraph")
groundhog.library(pkgs, "2022-07-19")

# Disabling scientific notation to improve readability of certain variables.
# options(scipen=0, digits=7) # To return to default setting
options(scipen = 999)

# Setting seed for results replication
set.seed(12345)

######Data Load##########

# Load the three datasets of removed tweets acquired from Twitter and one baseline dataset. Feb 2021 removal
# of Russian IRA accounts, Feb 2021 removal of Iranian accounts, and Dec 2021
# removal of PRC Xinjiang-focused accounts.

russia_clean <- read_csv("russia_2021_cleaned.csv")
iran_clean <- read_csv("iran_2021_cleaned.csv")
china_clean <- read_csv("china_2021_cleaned.csv")

# Datasets comes from Arunava Kumar Chakraborty on Kaggle: https://www.kaggle.com/datasets/arunavakrchakraborty/covid19-twitter-dataset
# These datasets contains tweets about COVID-19 from August to September 2020 and April to June 2020, respectively
covid_2020_1 <- read_csv("Covid-19 Twitter Dataset (Apr-Jun 2020).csv")
covid_2020_2 <- read_csv("Covid-19 Twitter Dataset (Aug-Sep 2020).csv")

######Pre-processing data#######
#Transforming user_mentions into split, long data to ensure 1 to 1 matching of accounts
russia_clean <- cSplit(russia_clean, "user_mentions", sep=",", direction = "long")

iran_clean <- cSplit(iran_clean, "user_mentions", sep=",", direction = "long")

china_clean <- cSplit(china_clean, "user_mentions", sep=",", direction = "long")

covid_2020_1 <- cSplit(covid_2020_1, "user_mentions", sep=",", direction = "long")

covid_2020_2 <- cSplit(covid_2020_2, "user_mentions", sep=",", direction = "long")

######Graph Creation######

# Creating a directed and undirected network graph of users mentioned by Russian info op

russia_net_d <- russia_clean %>%
  select(c(user_screen_name, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  filter(user_mentions != 0) %>% # removes zero values present
  graph_from_data_frame # creates network graph

russia_net_u <- russia_clean %>%
  select(c(user_screen_name, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  filter(user_mentions != 0) %>% # removes zero values present
  graph_from_data_frame(directed = FALSE) # creates undirected network graph

# Writing graph to a .gml file for visualization in Gephi
write_graph(simplify(russia_net_d),  "russia_net_m.gml", format = "gml")

### Steps to recreate Gephi visualization
# 1. Use ForceAtlas 2 layout with defaults in the layout menu in the lower left of the screen.
# 2. Add community detection using Louvain algorithm, with default settings, by selecting modularity in the statistics tab on the right. 
# 3. Select modularity class from the drop-down in the partition tab in the upper left of the screen, adding color to communities above 1%.
# 4. Select apply to add color based on community to the graph.
# 5. To detect and highlight influential accounts, adjust node appearance based on degree.
# 6. Set minimum size to 7 and maximum size to 70 to ensure nodes are visible.
# 7. Omitted labels as key nodes were account IDs hashed by Twitter.


# Creating a directed and undirected network graph of users mentioned by Iranian info op
iran_net_d <- iran_clean %>%
  select(c(user_screen_name, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  filter(user_mentions != 0) %>% # removes zero values present
  graph_from_data_frame # creates network graph

iran_net_u <- iran_clean %>%
  select(c(user_screen_name, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  filter(user_mentions != 0) %>% # removes zero values present
  graph_from_data_frame(directed = FALSE) # creates undirected network graph

# Writing graph to a .gml file for visualization in Gephi
write_graph(simplify(iran_net_d),  "iran_net_m.gml", format = "gml")

### Steps to recreate Gephi visualization
# 1. Use ForceAtlas 2 layout with defaults in the layout menu in the lower left of the screen.
# 2. Add community detection using Louvain algorithm, with default settings, by selecting modularity in the statistics tab on the right. 
# 3. Select modularity class from the drop-down in the partition tab in the upper left of the screen, adding color to communities above 1%.
# 4. Select apply to add color based on community to the graph.
# 5. To detect and highlight influential accounts, adjust node appearance based on degree.
# 6. Set minimum size to 16 and maximum size to 160 to ensure nodes are visible in such a complex graph.
# 7. Omitted labels as key nodes were account IDs hashed by Twitter.


# Creating a directed and undirected network graph of users mentioned by Chinese info op
china_net_d <- china_clean %>%
  select(c(user_screen_name, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  filter(user_mentions != 0) %>% # removes zero values present
  graph_from_data_frame # creates network graph

china_net_u <- china_clean %>%
  select(c(user_screen_name, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  filter(user_mentions != 0) %>% # removes zero values present
  graph_from_data_frame(directed = FALSE) # creates undirected network graph

# Writing graph to a .gml file for visualization in Gephi
write_graph(simplify(china_net_d),  "china_net_m.gml", format = "gml")

### Steps to recreate Gephi visualization
# 1. Use ForceAtlas 2 layout with gravity of 40.0 (to ensure both sub-networks were in the same view) in the layout menu in the lower left of the screen.
# 2. Add community detection using Louvain algorithm, with default settings, by selecting modularity in the statistics tab on the right. 
# 3. Select modularity class from the drop-down in the partition tab in the upper left of the screen, adding color to communities above 1%..
# 4. Select apply to add color based on community to the graph.
# 5. To detect and highlight influential accounts, adjust node appearance based on degree.
# 6. Set minimum size to 2 and maximum size to 30 to ensure nodes are visible.
# 7. Select Labels and Configuration and choose the Name option.
# 8. Filter names by degree size (min = 0.01, max = 5.5)
# 9. Converted Twitter IDs for some of the most prominently mentioned accounts to their @ name, if available.

# Creating a directed and undirected network graph of user interactions in baseline dataset

covid_net_d1 <- covid_2020_1 %>%
  select(c(original_author, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  graph_from_data_frame # creates network graph

covid_net_u1 <- covid_2020_1 %>%
  select(c(original_author, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  graph_from_data_frame(directed = FALSE) # creates undirected network graph

covid_net_d2 <- covid_2020_2 %>%
  select(c(original_author, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  graph_from_data_frame # creates network graph

covid_net_u2 <- covid_2020_2 %>%
  select(c(original_author, user_mentions)) %>% # selects userid and those mentioned
  filter(!is.na(user_mentions)) %>% # removes numerous NA values that would negatively impact graph
  graph_from_data_frame(directed = FALSE)

# Writing graph to a .gml file for visualization in Gephi
write_graph(simplify(covid_net_d1),  "covid_net_m.gml", format = "gml")

write_graph(simplify(covid_net_d2),  "covid_net2_m.gml", format = "gml")

### Steps to recreate Gephi visualizations
# 1. Use ForceAtlas 2 layout with gravity of 50.0 (to try to fit data onto the screen) in the layout menu in the lower left of the screen.
# 2. Add community detection using Louvain algorithm, with default settings, by selecting modularity in the statistics tab on the right. 
# 3. Select modularity class from the drop-down in the partition tab in the upper left of the screen, adding color to communities above 1%..
# 4. Select apply to add color based on community to the graph.
# 5. To detect and highlight influential accounts, adjust node appearance based on degree.
# 6. Set minimum size to 10 and maximum size to 100 to ensure nodes are visible.
# 7. Select Labels and Configuration and choose the Name option.
# 8. Filter names by degree size (min = 0.01, max = 5.5)
# 9. Converted Twitter IDs for some of the most prominently mentioned accounts to their @ name, if available.

######Analysis of networks##########

#### Russia Analysis ####

# Verifying successful creation of the network graph
summary(russia_net_d)
summary(russia_net_u)

# Creating random graphs to check if network structure is unique
g1 <- vector('list', 100) # Generate 100 random graphs

for(i in 1:100){
  g1[[i]] <- erdos.renyi.game(
    n = gorder(russia_net_d),
    p.or.m = edge_density(russia_net_d),
    type = "gnp"
  )
}

# Calculate average path length for the random graphs
g1.apls <- unlist(
  lapply(g1, mean_distance, directed = TRUE)
)

g1.apls <- unlist(
  lapply(g1, transitivity)
)

# Compare random graphs to the original network
hist(g1.apls, xlab = "Mean Distance", main = paste("Comparison of Mean Distance for the Russian Network and Random Networks"), xlim = range(c(1, 6)))
abline(
  v = mean_distance(
    russia_net_d, directed = TRUE
  ),
  col = "red",
  lty = 3,
  lwd = 2
) 

hist(g1.apls, xlim = range(c(0, 0.001)))
abline(
  v = transitivity(
    russia_net_d
  ),
  col = "red",
  lty = 3,
  lwd = 2
) 

mean(g1.apls) 
r3 <- signif(mean_distance(russia_net_d), 4)
mean_distance(russia_net_u)
r4 <- signif(transitivity(russia_net_d), 4) 
transitivity(russia_net_u) # the mean distance/transitivity of the random graphs are much higher than the observed graph, unlikely network distance came about by luck


# Remove any isolated nodes before analysis
isolated_ru <- which(degree(russia_net_d) == 0)
net_clean_ru <- delete.vertices(russia_net_d, isolated_ru) # Zero values were dropped

# Create data frame with six measures of centrality, betweeness, eigenvector centrality, Pagerank, degree in, degree out, and strength.
ru_cent<- data.frame(bet = betweenness(net_clean_ru), eig = centr_eigen(net_clean_ru)$vector,
                     p_rank = (page_rank(net_clean_ru)$vector), degr_in = degree(net_clean_ru, mode = "in"),
                     degr_out = degree(net_clean_ru, mode = "out"), stg = strength(net_clean_ru))
ru_cent <- cbind(account = rownames(ru_cent), ru_cent)
dim(ru_cent)

# Centrality results
top_n(ru_cent, 10) %>% # Betweeness
  arrange(desc(bet)) %>% 
  select(bet) 

top_n(ru_cent, 10, eig) %>% # Eigenvector
  arrange(desc(eig)) %>% 
  select(eig) 

top_n(ru_cent, 10, p_rank) %>%  # Pagerank
  arrange(desc(p_rank)) %>% 
  select(p_rank)

top_n(ru_cent, 10, degr_in) %>% # Degree in
  arrange(desc(degr_in)) %>% 
  select(degr_in) 

top_n(ru_cent, 10, degr_out) %>% # Degree out
  arrange(desc(degr_out)) %>% 
  select(degr_out) 

top_n(ru_cent, 10, stg) %>% # Strength
  arrange(desc(stg)) %>% 
  select(stg)

summary(ru_cent)
summary(ru_cent$bet)

# Modularity analysis

leiden_ru <- cluster_leiden( #clustering algorithm
  russia_net_u,
  objective_function = "modularity",
  weights = E(russia_net_u)$weight
)

V(russia_net_u)$cluster <- membership(leiden_ru) #assign vertices to clusters

r5 <- signif(modularity( #modularity of the network
  russia_net_u,
  membership = V(russia_net_u)$cluster,
  weights = E(russia_net_u)$weight
), 4)

r6 <- signif(max(V(russia_net_u)$cluster), 4) #number of communities

# Russia Undirected Cluster visualization
(ru1 <- ggraph(russia_net_u, layout = "fr") +
    geom_edge_link(color = "grey", alpha = 0.7) + 
    geom_node_point(aes(color = factor(cluster)), size = 1) +
    labs(color = "Leiden community") +
    theme_void())

ru1

# Identifying most important accounts based on degree centrality
(highest_degree_ru <- tibble(community = 1:max(V(russia_net_u)$cluster)) %>% 
    rowwise() %>% 
    mutate(
      highest_degree = induced_subgraph(
        russia_net_u,
        vids = V(russia_net_u)[V(russia_net_u)$cluster == community]) %>% 
        degree() %>% 
        which.max() %>% 
        names(),
      community_size = V(russia_net_u)[V(russia_net_u)$cluster == community] %>% 
        length()
    ))

# Descriptive analysis
r1 <- gorder(net_clean_ru)

r2 <- gsize(net_clean_ru)

r7 <- farthest_vertices(net_clean_ru)
r7 <- r7$distance

get_diameter(net_clean_ru, directed = TRUE)

r8 <- signif(edge_density(net_clean_ru), 4) # higher means greater interconnection

mean_distance(net_clean_ru, directed = TRUE) # Shorter means greater interconnection

transitivity(net_clean_ru) # higher means higher probability adjacent vertices of a given vertex are connected

r9 <- reciprocity(net_clean_ru)

r10 <- signif(assortativity.degree(net_clean_ru, directed = TRUE), 4)

rc = edge.betweenness.community(net_clean_ru)
sizes(rc)

# Weighted betweenness
dist_weight_ru <- 1 / E(net_clean_ru)$weight
ebr <- edge_betweenness(net_clean_ru, weights = dist_weight_ru)
ebr
summary(ebr)
dim(ebr)

###### Iran Analysis ####

# Verifying successful creation of the network graph
summary(iran_net_d)
summary(iran_net_u)

# Creating random graphs to check if network structure is unique
g2 <- vector('list', 100) # Generate 100 random graphs

for(i in 1:100){
  g2[[i]] <- erdos.renyi.game(
    n = gorder(iran_net_d),
    p.or.m = edge_density(iran_net_d),
    type = "gnp"
  )
}

# Calculate average path length for the random graphs
g2.apls <- unlist(
  lapply(g2, mean_distance, directed = TRUE)
)

#g2.apls <- unlist(
#  lapply(g2, transitivity)
#)

# Compare random graphs to the original network
hist(g2.apls, xlim = range(c(0, 7)))
abline(
  v = mean_distance(
    iran_net_d, directed = TRUE
  ),
  col = "red",
  lty = 3,
  lwd = 2
) 

hist(g2.apls, xlim = range(c(0, 0.0002)))
abline(
  v = transitivity(
    iran_net_d
  ),
  col = "red",
  lty = 3,
  lwd = 2
) 

mean(g2.apls) 
i3 <- signif(mean_distance(iran_net_d), 4)
mean_distance(iran_net_u)
i4 <- signif(transitivity(iran_net_d), 4) 
transitivity(iran_net_u) # the mean distance/transitivity of the random graphs are lower than the observed graph, unlikely network distance came about by luck

# Remove any isolated nodes before analysis
isolated_ir <- which(degree(iran_net_d) == 0)
net_clean_ir <- delete.vertices(iran_net_d, isolated_ir) # Zero values were dropped

# Create data frame with six measures of centrality, betweeness, eigenvector centrality, Pagerank, degree in, degree out, and strength.
ir_cent<- data.frame(bet = betweenness(net_clean_ir),eig = centr_eigen(net_clean_ir)$vector,
                     p_rank = (page_rank(net_clean_ir)$vector), degr_in = degree(net_clean_ir, mode = "in"),
                     degr_out = degree(net_clean_ir, mode = "out"), stg = strength(net_clean_ir))
ir_cent <- cbind(account = rownames(ir_cent), ir_cent)

# Centrality results
top_n(ir_cent, 10) %>% # Betweeness
  arrange(desc(bet)) %>% 
  select(bet) 

top_n(ir_cent, 10, eig) %>% # Eigenvector
  arrange(desc(eig)) %>% 
  select(eig) 

top_n(ir_cent, 10, p_rank) %>%  # Pagerank
  arrange(desc(p_rank)) %>% 
  select(p_rank)

top_n(ir_cent, 10, degr_in) %>% # Degree in
  arrange(desc(degr_in)) %>% 
  select(degr_in) 

top_n(ir_cent, 10, degr_out) %>% # Degree out
  arrange(desc(degr_out)) %>% 
  select(degr_out) 

top_n(ir_cent, 10, stg) %>% # Strength
  arrange(desc(stg)) %>% 
  select(stg)

summary(ir_cent)

# Modularity analysis

leiden_ir <- cluster_leiden( #clustering algorithm
  iran_net_u,
  objective_function = "modularity",
  weights = E(iran_net_u)$weight
)

V(iran_net_u)$cluster <- membership(leiden_ir) #assign vertices to clusters

i5 <- signif(modularity( #modularity of the network
  iran_net_u,
  membership = V(iran_net_u)$cluster,
  weights = E(iran_net_u)$weight
), 4)

i6 <- signif(max(V(iran_net_u)$cluster), 4) #number of communities

# Iran Undirected Cluster visualization
(ir1 <- ggraph(iran_net_u, layout = "fr") +
    geom_edge_link(color = "grey", alpha = 0.7) + 
    geom_node_point(aes(color = factor(cluster)), size = 1) +
    labs(color = "Leiden community") +
    theme_void())

ir1

# Identifying most important accounts based on degree centrality
(highest_degree_ir <- tibble(community = 1:max(V(iran_net_u)$cluster)) %>% 
    rowwise() %>% 
    mutate(
      highest_degree = induced_subgraph(
        iran_net_u,
        vids = V(iran_net_u)[V(iran_net_u)$cluster == community]) %>% 
        degree() %>% 
        which.max() %>% 
        names(),
      community_size = V(iran_net_u)[V(iran_net_u)$cluster == community] %>% 
        length()
    ))

# Descriptive analysis
i1 <- gorder(net_clean_ir)

i2 <- gsize(net_clean_ir)

i7 <- farthest_vertices(net_clean_ir)
i7 <- i7$distance

get_diameter(net_clean_ir, directed = TRUE)

i8 <- signif(edge_density(net_clean_ir), 4) # higher means greater interconnection

mean_distance(net_clean_ir, directed = TRUE) # Shorter means greater interconnection

transitivity(net_clean_ir) # higher means higher probability adjacent vertices of a given vertex are connected

i9 <- signif(reciprocity(net_clean_ir), 4)

i10 <- signif(assortativity.degree(net_clean_ir, directed = TRUE), 4)

ic = edge.betweenness.community(net_clean_ir)
sizes(ic)

# Weighted betweenness
dist_weight_ir <- 1 / E(net_clean_ir)$weight
ebi <- edge_betweenness(net_clean_ir, weights = dist_weight_ir)
ebi
summary(ebi)


###### China Analysis ####

# Verifying successful creation of the network graph
summary(china_net_d) #Note: despite Chinese IO having the most accounts, it seems to have the fewest mentions of other accounts, by a significant amount.
summary(china_net_u)

# Creating random graphs to check if network structure is unique
g3 <- vector('list', 100) # Generate 100 random graphs

for(i in 1:100){
  g3[[i]] <- erdos.renyi.game(
    n = gorder(china_net_d),
    p.or.m = edge_density(china_net_d),
    type = "gnp"
  )
}

# Calculate average path length for the random graphs
g3.apls <- unlist(
  lapply(g3, mean_distance, directed = TRUE)
)

g3.apls <- unlist(
  lapply(g3, transitivity)
)

# Compare random graphs to the original network
hist(g3.apls, xlim = range(c(1, 6)))
abline(
  v = mean_distance(
    china_net_d, directed = TRUE
  ),
  col = "red",
  lty = 3,
  lwd = 2
) 

hist(g3.apls, xlim = range(c(0, 0.005)))
abline(
  v = transitivity(
    china_net_d
  ),
  col = "red",
  lty = 3,
  lwd = 2
) 

mean(g3.apls) 
c3 <- signif(mean_distance(china_net_d), 4)
mean_distance(china_net_u)
c4 <- signif(transitivity(china_net_d), 4) 
transitivity(china_net_u) # the mean distance of the random graphs are much higher than the observed graph, unlikely network distance came about by luck.
# Transitivity of the observed graph is closer to the random graphs but has a different mean than the mean of the random graphs.

# Remove any isolated nodes before analysis
isolated_ch <- which(degree(china_net_d) == 0)
net_clean_ch <- delete.vertices(china_net_d, isolated_ch) # Zero values were dropped

# Create data frame with six measures of centrality, betweeness, eigenvector centrality, Pagerank, degree in, degree out, and strength.
ch_cent<- data.frame(bet = betweenness(net_clean_ch),eig = centr_eigen(net_clean_ch)$vector,
                     p_rank = (page_rank(net_clean_ch)$vector), degr_in = degree(net_clean_ch, mode = "in"),
                     degr_out = degree(net_clean_ch, mode = "out"), stg = strength(net_clean_ch))
ch_cent <- cbind(account = rownames(ch_cent), ch_cent)

# Centrality results
top_n(ch_cent, 10) %>% # Betweeness
  arrange(desc(bet)) %>% 
  select(bet) 

top_n(ch_cent, 10, eig) %>% # Eigenvector
  arrange(desc(eig)) %>% 
  select(eig) 

top_n(ch_cent, 10, p_rank) %>%  # Pagerank
  arrange(desc(p_rank)) %>% 
  select(p_rank)

top_n(ch_cent, 10, degr_in) %>% # Degree in
  arrange(desc(degr_in)) %>% 
  select(degr_in) 

top_n(ch_cent, 10, degr_out) %>% # Degree out
  arrange(desc(degr_out)) %>% 
  select(degr_out) 

top_n(ch_cent, 10, stg) %>% # Strength
  arrange(desc(stg)) %>% 
  select(stg)

summary(ch_cent)

# Modularity analysis

leiden_ch <- cluster_leiden( #clustering algorithm
  china_net_u,
  objective_function = "modularity",
  weights = E(china_net_u)$weight
)

V(china_net_u)$cluster <- membership(leiden_ch) #assign vertices to clusters

c5 <- signif(modularity( #modularity of the network
  china_net_u,
  membership = V(china_net_u)$cluster,
  weights = E(china_net_u)$weight
), 4)

c6 <- signif(max(V(china_net_u)$cluster), 4) #number of communities

# China Undirected Cluster visualization
(ch1 <- ggraph(china_net_u, layout = "fr") +
    geom_edge_link(color = "grey", alpha = 0.7) + 
    geom_node_point(aes(color = factor(cluster)), size = 1) +
    labs(color = "Leiden community") +
    theme_void())

ch1

# Identifying most important accounts based on degree centrality
(highest_degree_ch <- tibble(community = 1:max(V(china_net_u)$cluster)) %>% 
    rowwise() %>% 
    mutate(
      highest_degree = induced_subgraph(
        china_net_u,
        vids = V(china_net_u)[V(china_net_u)$cluster == community]) %>% 
        degree() %>% 
        which.max() %>% 
        names(),
      community_size = V(china_net_u)[V(china_net_u)$cluster == community] %>% 
        length()
    ))

# Descriptive analysis
c1 <- gorder(net_clean_ch)

c2 <- gsize(net_clean_ch)

c7 <- farthest_vertices(net_clean_ch)
c7 <- c7$distance

get_diameter(net_clean_ch, directed = TRUE)

c8 <- signif(edge_density(net_clean_ch), 4) # higher means greater interconnection

mean_distance(net_clean_ch, directed = TRUE) # Shorter means greater interconnection

transitivity(net_clean_ch) # higher means higher probability adjacent vertices of a given vertex are connected

c9 <- reciprocity(net_clean_ch)

c10 <- signif(assortativity.degree(net_clean_ch, directed = TRUE), 4)

cc = edge.betweenness.community(net_clean_ch)
sizes(cc)

# Weighted betweenness
dist_weight_ch <- 1 / E(net_clean_ch)$weight
ebc <- edge_betweenness(net_clean_ch, weights = dist_weight_ch)
ebc
summary(ebc)

##### Baseline network metrics ####

summary(covid_net_d1)
summary(covid_net_u1)
summary(covid_net_d2)
summary(covid_net_u2)

covid1_3 <- signif(mean_distance(covid_net_d1), 4)
covid2_3 <- signif(mean_distance(covid_net_d2), 4)
mean_distance(covid_net_u1)
mean_distance(covid_net_u2)
covid1_4 <- signif(transitivity(covid_net_d1), 4)
covid2_4 <- signif(transitivity(covid_net_d2), 4)
transitivity(covid_net_u1)
transitivity(covid_net_u2) # the mean distance/transitivity of the random graphs are much higher than the observed graph, unlikely network distance came about by luck

# Remove any isolated nodes before analysis
isolated_cv1 <- which(degree(covid_net_d1) == 0)
net_clean_cv1 <- delete.vertices(covid_net_d1, isolated_cv1) # Zero values were dropped

isolated_cv2 <- which(degree(covid_net_d2) == 0)
net_clean_cv2 <- delete.vertices(covid_net_d2, isolated_cv2) # Zero values were dropped

# Create data frame with six measures of centrality, betweeness, eigenvector centrality, Pagerank, degree in, degree out, and strength.
cv1_cent<- data.frame(bet = betweenness(net_clean_cv1),eig = centr_eigen(net_clean_cv1)$vector,
                     p_rank = (page_rank(net_clean_cv1)$vector), degr_in = degree(net_clean_cv1, mode = "in"),
                     degr_out = degree(net_clean_cv1, mode = "out"), stg = strength(net_clean_cv1))
cv1_cent <- cbind(account = rownames(cv1_cent), cv1_cent)

cv2_cent<- data.frame(bet = betweenness(net_clean_cv2),eig = centr_eigen(net_clean_cv2)$vector,
                     p_rank = (page_rank(net_clean_cv2)$vector), degr_in = degree(net_clean_cv2, mode = "in"),
                     degr_out = degree(net_clean_cv2, mode = "out"), stg = strength(net_clean_cv2))
cv2_cent <- cbind(account = rownames(cv2_cent), cv2_cent)

# Centrality results
top_n(cv1_cent, 10) %>% # Betweeness
  arrange(desc(bet)) %>% 
  select(bet) 

top_n(cv1_cent, 10, eig) %>% # Eigenvector
  arrange(desc(eig)) %>% 
  select(eig) 

top_n(cv1_cent, 10, p_rank) %>%  # Pagerank
  arrange(desc(p_rank)) %>% 
  select(p_rank)

top_n(cv1_cent, 10, degr_in) %>% # Degree in
  arrange(desc(degr_in)) %>% 
  select(degr_in) 

top_n(cv1_cent, 10, degr_out) %>% # Degree out
  arrange(desc(degr_out)) %>% 
  select(degr_out) 

top_n(cv1_cent, 10, stg) %>% # Strength
  arrange(desc(stg)) %>% 
  select(stg)

summary(cv1_cent)

top_n(cv2_cent, 10) %>% # Betweeness
  arrange(desc(bet)) %>% 
  select(bet) 

top_n(cv2_cent, 10, eig) %>% # Eigenvector
  arrange(desc(eig)) %>% 
  select(eig) 

top_n(cv2_cent, 10, p_rank) %>%  # Pagerank
  arrange(desc(p_rank)) %>% 
  select(p_rank)

top_n(cv2_cent, 10, degr_in) %>% # Degree in
  arrange(desc(degr_in)) %>% 
  select(degr_in) 

top_n(cv2_cent, 10, degr_out) %>% # Degree out
  arrange(desc(degr_out)) %>% 
  select(degr_out) 

top_n(cv2_cent, 10, stg) %>% # Strength
  arrange(desc(stg)) %>% 
  select(stg)

summary(cv2_cent)

# Modularity analysis

leiden_cv1 <- cluster_leiden( #clustering algorithm
  covid_net_u1,
  objective_function = "modularity",
  weights = E(covid_net_u1)$weight
)

V(covid_net_u1)$cluster <- membership(leiden_cv1) #assign vertices to clusters

covid1_5 <- signif(modularity( #modularity of the network
  covid_net_u1,
  membership = V(covid_net_u1)$cluster,
  weights = E(covid_net_u1)$weight
), 4)

covid1_6 <- signif(max(V(covid_net_u1)$cluster), 4) #number of communities

leiden_cv2 <- cluster_leiden( #clustering algorithm
  covid_net_u2,
  objective_function = "modularity",
  weights = E(covid_net_u2)$weight
)

V(covid_net_u2)$cluster <- membership(leiden_cv2) #assign vertices to clusters

covid2_5 <- signif(modularity( #modularity of the network
  covid_net_u2,
  membership = V(covid_net_u2)$cluster,
  weights = E(covid_net_u2)$weight
), 4)

covid2_6 <- signif(max(V(covid_net_u2)$cluster), 4) #number of communities

# Covid Undirected Cluster visualization
(covid1 <- ggraph(covid_net_u1, layout = "fr") +
    geom_edge_link(color = "grey", alpha = 0.7) + 
    geom_node_point(aes(color = factor(cluster)), size = 1) +
    labs(color = "Leiden community") +
    theme_void())

covid1

(covid2 <- ggraph(covid_net_u2, layout = "fr") +
    geom_edge_link(color = "grey", alpha = 0.7) + 
    geom_node_point(aes(color = factor(cluster)), size = 1) +
    labs(color = "Leiden community") +
    theme_void())

covid2

# Identifying most important accounts based on degree centrality
(highest_degree_cv1 <- tibble(community = 1:max(V(covid_net_u1)$cluster)) %>% 
    rowwise() %>% 
    mutate(
      highest_degree = induced_subgraph(
        covid_net_u1,
        vids = V(covid_net_u1)[V(covid_net_u1)$cluster == community]) %>% 
        degree() %>% 
        which.max() %>% 
        names(),
      community_size = V(covid_net_u1)[V(covid_net_u1)$cluster == community] %>% 
        length()
    ))

(highest_degree_cv2 <- tibble(community = 1:max(V(covid_net_u2)$cluster)) %>% 
    rowwise() %>% 
    mutate(
      highest_degree = induced_subgraph(
        covid_net_u2,
        vids = V(covid_net_u2)[V(covid_net_u2)$cluster == community]) %>% 
        degree() %>% 
        which.max() %>% 
        names(),
      community_size = V(covid_net_u2)[V(covid_net_u2)$cluster == community] %>% 
        length()
    ))

# Descriptive analysis
covid1_1 <- gorder(net_clean_cv1)

covid2_1 <- gorder(net_clean_cv2)

covid1_2 <- gsize(net_clean_cv1)

covid2_2 <- gsize(net_clean_cv2)

covid1_7 <- farthest_vertices(net_clean_cv1)
covid1_7 <- covid1_7$distance

covid2_7 <- farthest_vertices(net_clean_cv2)
covid2_7 <- covid2_7$distance

get_diameter(net_clean_cv1, directed = TRUE)
get_diameter(net_clean_cv2, directed = TRUE)

covid1_8 <- signif(edge_density(net_clean_cv1), 4) # higher means greater interconnection

covid2_8 <- signif(edge_density(net_clean_cv2), 4) # higher means greater interconnection

mean_distance(net_clean_cv1, directed = TRUE) # Shorter means greater interconnection

mean_distance(net_clean_cv2, directed = TRUE) # Shorter means greater interconnection

transitivity(net_clean_cv1) # higher means higher probability adjacent vertices of a given vertex are connected

transitivity(net_clean_cv2) # higher means higher probability adjacent vertices of a given vertex are connected

covid1_9 <- signif(reciprocity(net_clean_cv1), 4)

covid2_9 <- signif(reciprocity(net_clean_cv2), 4)

covid1_10 <- signif(assortativity.degree(net_clean_cv1, directed = TRUE), 4)

covid2_10 <- signif(assortativity.degree(net_clean_cv2, directed = TRUE), 4)

cvc1 = edge.betweenness.community(net_clean_cv1)
sizes(cvc1)

cvc2 = edge.betweenness.community(net_clean_cv1)
sizes(cvc2)

# Weighted betweenness
dist_weight_cv1 <- 1 / E(net_clean_cv1)$weight
ebcv1 <- edge_betweenness(net_clean_cv1, weights = dist_weight_cv1)
ebcv1
summary(ebcv1)

dist_weight_cv2 <- 1 / E(net_clean_cv2)$weight
ebcv2 <- edge_betweenness(net_clean_cv2, weights = dist_weight_cv2)
ebcv2
summary(ebcv2)

##### Comparisons of similarity across networks #######

# function for jaccard similarity of edge sets
jaccard_edgeset_similarity <- function(G1, G2) {
  inter <- length(E(G1 %s% G2))
  un <- length(E(G1 %u% G2))
  
  if (un == 0) {
    0
  } else {
    inter/un
  }
}

# Test China/Russia
jaccard_edgeset_similarity(china_net_d, russia_net_d)

# Test China/Iran
jaccard_edgeset_similarity(china_net_d, iran_net_d)

# Test Russia/Iran
jaccard_edgeset_similarity(russia_net_d, iran_net_d)

# All three tests show zero similar edgesets between the datasets, suggesting no cross-interaction within the networks.

##### Creating tables of results for each network #####

# Create results dataframe for network dimensions
network_structures <- data.frame(Descriptors = c('Vertices', 'Edges', 'Mean Distance', 
                                       'Transitivity', 'Modularity', 'Communities', 
                                       'Farthest Vertices', 'Edge Density', 'Reciprocity', 
                                       'Assortativity'),
                           Russian_Network = c(as.character(r1), r2, r3, r4, r5, r6, r7, r8, r9, r10),
                           Iranian_Network = c(as.character(i1), i2, i3, i4, i5, i6, i7, i8, i9, i10),
                           Chinese_Network = c(as.character(c1), c2, c3, c4, c5, c6, c7, c8, c9, c10),
                           Control_1 = c(as.character(covid1_1), covid1_2, covid1_3, covid1_4,
                                                     covid1_5, covid1_6, covid1_7, covid1_8, covid1_9,
                                                     covid1_10),
                           Control_2 = c(as.character(covid2_1), covid2_2, covid2_3, covid2_4,
                                                     covid2_5, covid2_6, covid2_7, covid2_8, covid2_9,
                                                     covid2_10)) 
                           # coerced lists to characters for table formatting purposes

write_csv(network_structures, "mentions_network_results.csv") #to produce table for report

kbl(network_structures, format = "html", align = "c") %>%
  kable_classic(full_width = F, html_font = "Garamond") %>%
  row_spec(0, bold = T)

network_results <- read_csv("mentions_network_results2.csv", show_col_types = FALSE) #transposed values from columns to rows for easier analysis

grouped_results <- read_csv("mentions_network_results3.csv", show_col_types = FALSE) # transposed values with a group column added for "IO" and "control"

### Testing significance of key network statistics

# Comparing results visually to see if there appear to be significant differences between the Twitter networks and the controls.
# In each case there appears to be at least some noticeable difference between values
plotmeans(vertices_ ~ group, data = network_results)
plotmeans(edges_ ~ group, data = network_results)
plotmeans(mean_distance_ ~ group, data = network_results)
plotmeans(transitivity_ ~ group, data = network_results)
plotmeans(modularity_ ~ group, data = network_results)
plotmeans(community_size_ ~ group, data = network_results)
plotmeans(farthest_vertex_ ~ group, data = network_results)
plotmeans(edge_density_ ~ group, data = network_results)
plotmeans(reciprocity_ ~ group, data = network_results)
plotmeans(assortativity_ ~ group, data = network_results)

# But...those values also need to be tested for statistical significance.
# Many of the values were quite close, so I suspected at least some would not pass statistical significance testing.
# However, because each igraph metric above only outputs a single value, trying to test for significance was a challenge, so multiple metrics were used to double check.

# Comparing mean values of each key network statistics from the info op networks against the mean values from the control networks using Welch Two Sample t-test
# Pass criteria p-value <= 0.05
t.test(network_results$mean_distance_[1:3], network_results$mean_distance_[4:5], alternative = "two.sided") # fail 
t.test(network_results$transitivity_[1:3], network_results$transitivity_[4:5], alternative = "two.sided") # fail 
t.test(network_results$modularity_[1:3], network_results$modularity_[4:5], alternative = "two.sided") # fail 
t.test(network_results$community_size_[1:3], network_results$community_size_[4:5], alternative = "two.sided") # fail 
t.test(network_results$farthest_vertex_[1:3], network_results$farthest_vertex_[4:5], alternative = "two.sided") # pass 
t.test(network_results$edge_density_[1:3], network_results$edge_density_[4:5], alternative = "two.sided") # fail 
t.test(network_results$reciprocity_[1:3], network_results$reciprocity_[4:5], alternative = "two.sided") # fail 
t.test(network_results$assortativity_[1:3], network_results$assortativity_[4:5], alternative = "two.sided") # fail 

# Creating one-way ANOVA models to conduct grouped tests for the "IO" group versus the "control" group to validate farthest vertex pass amid all others failing. Pass criteria <= 0.05
model1 <- aov(mean_distance_ ~ group, data = grouped_results)
summary(model1) # fail

model2 <- aov(transitivity_ ~ group, data = grouped_results)
summary(model2) #fail

model3 <- aov(modularity_ ~ group, data = grouped_results)
summary(model3) # fail

model4 <- aov(community_size_ ~ group, data = grouped_results)
summary(model4) # pass

model5 <- aov(farthest_vertex_ ~ group, data = grouped_results)
summary(model5) # pass

model6 <- aov(edge_density_ ~ group, data = grouped_results)
summary(model6) # fail

model7 <- aov(reciprocity_ ~ group, data = grouped_results)
summary(model7) # fail

model8 <- aov(assortativity_ ~ group, data = grouped_results)
summary(model8) # fail

# Using Conover's test of multiple comparisons to try to attempt to make pair-wise comparisons
# Pass criteria p-value <= 0.05
ConoverTest(network_results$mean_distance_, as.factor(network_results$group)) # NA- p-value unable to be established
ConoverTest(network_results$transitivity_, as.factor(network_results$group)) # NA- p-value unable to be established
ConoverTest(network_results$modularity_, as.factor(network_results$group)) # NA- p-value unable to be established
ConoverTest(network_results$community_size_, as.factor(network_results$group)) # NA- p-value unable to be established 
ConoverTest(network_results$farthest_vertex_, as.factor(network_results$group)) # NA- p-value unable to be established
ConoverTest(network_results$edge_density_, as.factor(network_results$group)) # NA- p-value unable to be established 
ConoverTest(network_results$reciprocity_, as.factor(network_results$group)) # NA- p-value unable to be established
ConoverTest(network_results$assortativity_, as.factor(network_results$group)) # NA- p-value unable to be established 
# Suspect small sample size (n=5) may be at play. 

### Testing significance of centrality measurements

# Adding group assignments to each centrality dataframe
ir_cent <- cbind(group = rep("IO"), ir_cent)
ru_cent <- cbind(group = rep("IO"), ru_cent)
ch_cent <- cbind(group = rep("IO"), ch_cent)
cv1_cent <- cbind(group = rep("control"), cv1_cent)
cv2_cent <- cbind(group = rep("control"), cv2_cent)

# Adding network labels so the values can be differentiated in the large dataframe
ir_cent <- cbind(network = rep("Iran"), ir_cent)
ru_cent <- cbind(network = rep("Russia"), ru_cent)
ch_cent <- cbind(network = rep("China"), ch_cent)
cv1_cent <- cbind(network = rep("Control_1"), cv1_cent)
cv2_cent <- cbind(network = rep("Control_2"), cv2_cent)

# Combining the smaller centrality dataframes into one super-sized dataframe
df_list <- list(ir_cent, ru_cent, ch_cent, cv1_cent, cv2_cent)
super_cent <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)

skim_without_charts(super_cent)
summary(super_cent)

#Creating ANOVA models between IO networks and Controls
# Betweenness test  
model_bet <- super_cent %>%
  aov(bet ~ group, data = .)
summary(model_bet) #pass
# Eigenvector test
model_eig <- super_cent %>%
  aov(eig ~ group, data = .)
summary(model_eig) #pass
# Pagerank test
model_p_rank <- super_cent %>%
  aov(p_rank ~ group, data = .)
summary(model_p_rank) #pass
# Degree in test
model_deg_in <- super_cent %>%
  aov(degr_in ~ group, data = .)
summary(model_deg_in) #pass
# Degree out test
model_deg_out <- super_cent %>%
  aov(degr_out ~ group, data = .)
summary(model_deg_out) #pass
# Strength test
model_stg <- super_cent %>%
  aov(stg ~ group, data = .)
summary(model_stg) #pass

# Creating network-specific models
# Betweenness test  
network_bet <- super_cent %>%
  aov(bet ~ network, data = .)
summary(network_bet) #pass
# Eigenvector test
network_eig <- super_cent %>%
  aov(eig ~ network, data = .)
summary(network_eig) #pass
# Pagerank test
network_p_rank <- super_cent %>%
  aov(p_rank ~ network, data = .)
summary(network_p_rank) #pass
# Degree in test
network_deg_in <- super_cent %>%
  aov(degr_in ~ network, data = .)
summary(network_deg_in) #pass
# Degree out test
network_deg_out <- super_cent %>%
  aov(degr_out ~ network, data = .)
summary(network_deg_out) #pass
# Strength test
network_stg <- super_cent %>%
  aov(stg ~ network, data = .)
summary(network_stg) #pass

# Running a Tukey multiple comparisons of means test to see which networks are different
# Pairing pass criteria adjusted p-value <= 0.05
TukeyHSD(network_bet) # Iran-Control_1, Iran-Control_2, and Russia-Iran Pass, all else fail
TukeyHSD(network_eig) # All pass except Control_2-Control_1
TukeyHSD(network_p_rank) # All pass except Control_2-Control_1
TukeyHSD(network_deg_in) # Russia-Iran, Russia-China, Iran-China, and Control_2-China fail, all else pass
TukeyHSD(network_deg_out) # Russia-Control_2, Iran-Control_2, Russia-Control_1, Iran-Control_1 pass, all else fail
TukeyHSD(network_stg) # Control_2-China, Iran-China, Russia-China, Russia-Iran fail, all else pass

# Plotting the 95% confidence interval of the family-wise confidence for each pairing
with(par(mai=c(1,2.5,1,1)),{plot(TukeyHSD(network_bet, conf.level=.95), las=1,cex.axis=0.8)})

with(par(mai=c(1,2.5,1,1)),{plot(TukeyHSD(network_eig, conf.level=.95), las=1,cex.axis=0.8)})

with(par(mai=c(1,2.5,1,1)),{plot(TukeyHSD(network_p_rank, conf.level=.95), las=1,cex.axis=0.8)})

with(par(mai=c(1,2.5,1,1)),{plot(TukeyHSD(network_deg_in, conf.level=.95), las=1,cex.axis=0.8)})

with(par(mai=c(1,2.5,1,1)),{plot(TukeyHSD(network_deg_out, conf.level=.95), las=1,cex.axis=0.8)})

with(par(mai=c(1,2.5,1,1)),{plot(TukeyHSD(network_stg, conf.level=.95), las=1,cex.axis=0.8)})

