## Load required packages
library(FKmL)
library(dplyr); library(ggplot2); library(ggpubr); library(reshape2); library(flexclust)

## -------------------------------------------------------------------------- ##
## Generate trajectories
## -------------------------------------------------------------------------- ##
## Randomly choose number of trajectories (between 20 and 30) in cluster
## Randomly select time points (between 5 and 20) per trajectory within a time range of 1 to 20

## cluster 1
set.seed(100)
cl1_traj <- NULL
cnt <- 0

for (i in 1:sample(c(20:30), 1)) {
  cnt <- cnt + 1
  set.seed(100 + cnt)
  time <- sample(c(5:20), 1)
  
  gen_traj <- data.frame(
    id = rep(i, time),
    time = sample(c(1:20), size = time) %>% sort(),
    group = rep(1, time),
    
    var1 = seq(2, 18, length = time) + rnorm(time, runif(1, -1, 1), 2),
    var2 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(1, 5, 3) + rnorm(time, 10, 2),
    
    var3 = (1 - pnorm(1:time, runif(1, 3, 7), 1)) * rnorm(1, 5, 2),
    var4 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(1, 15, 2),
    
    var5 = (1 - pnorm(1:time, runif(1, 3, 7), 1)) * rnorm(1, 10, 6),
    var6 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(1, 40, 10),
    
    var7 = dnorm(1:time, runif(1, 3, 5), 1) * rnorm(1, 35, 2),
    
    var8 = rnorm(time, 5, 3),
    var9 = rnorm(time, 5, 3)
  )
  cl1_traj <- rbind(cl1_traj, gen_traj)
}

cl1_traj %>% head(20)
cl1_traj %>% group_by(id) %>%  summarise(n = n())

## cluster 2
n_cl1 <- unique(cl1_traj$id) %>% length
set.seed(200)
cl2_traj <- NULL
cnt <- 0

for (i in 1:sample(c(20:30), 1)) {
  cnt <- cnt + 1
  set.seed(200 + cnt)
  time <- sample(c(5:20), 1)
  
  gen_traj <- data.frame(
    id = rep(n_cl1 + i, time),
    time = sample(c(1:20), size = time) %>% sort(),
    group = rep(2, time),
    
    var1 = 10 + rnorm(time, runif(1, -1, 1), 2),
    var2 = (dnorm(1:time, runif(1, 2, 8), 1) + dnorm(1:time, runif(1, 7, 13), 1) + dnorm(1:time, runif(1, 12, 18), 1)) * rnorm(time, 30, 3),
    
    var3 = (1 - pnorm(1:time, runif(1, 3, 7), 1)) * rnorm(1, 10, 2),
    var4 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(1, 30, 2),
    
    var5 = (1 - pnorm(1:time, runif(1, 3, 7), 1)) * rnorm(1, 20, 6),
    var6 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(1, 80, 10),
    
    var7 = dnorm(1:time, runif(1, 6, 8), 1) * rnorm(1, 45, 2),
    
    var8 = rnorm(time, 5, 3),
    var9 = rnorm(time, 5, 3)
  )
  cl2_traj <- rbind(cl2_traj, gen_traj)
}

cl2_traj %>% group_by(id) %>%  summarise(n = n())

## cluster 3
n_cl2 <- unique(cl2_traj$id) %>% length
set.seed(300)
cl3_traj <- NULL
cnt <- 0

for (i in 1:sample(c(20:30), 1)) {
  cnt <- cnt + 1
  set.seed(300 + cnt)
  time <- sample(c(5:20), 1)
  
  gen_traj <- data.frame(
    id = rep(n_cl1 + n_cl2 + i, time) %>% factor,
    time = sample(c(1:20), size = time) %>% sort(),
    group = rep(3, time),
    
    var1 = seq(15, 5, length= time) + rnorm(time, runif(1, -1, 1), 2),
    var2 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(time, 45, 3) + dnorm(1:time, runif(1, 12, 18), 1) * rnorm(time, 35, 3),
    
    var3 = (1 - pnorm(1:time, runif(1, 3, 7), 1)) * rnorm(1, 15, 2),
    var4 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(1, 45, 2),
    
    var5 = (1 - pnorm(1:time, runif(1, 3, 7), 1)) * rnorm(1, 30, 6),
    var6 = dnorm(1:time, runif(1, 2, 8), 1) * rnorm(1, 120, 10),
    
    var7 = (0.4 * dnorm(1:time, runif(1, 3, 5), 1) + 0.6 * dnorm(1:time, runif(1, 6, 8), 1)) * rnorm(1, 40, 2),
    
    var8 = rnorm(time, 5, 3),
    var9 = rnorm(time, 5, 3)
  )
  cl3_traj <- rbind(cl3_traj, gen_traj)
}

cl3_traj %>% group_by(id) %>%  summarise(n = n())

## Combine all trajectories
df_traj <- rbind(cl1_traj, cl2_traj, cl3_traj)
df_traj$id %>% unique %>% length  ## 83 trajectories
df_traj %>% head

n_traj <- df_traj %>% group_by(id) %>% summarise(n = n())
n_traj$n %>% table

## Base trajectories for each cluster
generate_base_traj <- function(group, time_size) {
  data.frame(
    time = seq(1, 20, length = time_size),
    group = group,
    
    var1 = if (group == 1) seq(2, 18, length = time_size) else if (group == 2) rep(10, time_size) else seq(15, 5, length = time_size),
    
    var2 = if (group == 1) dnorm(1:time_size, 5, 1) * 5 + 10 else if (group == 2) (dnorm(1:time_size, 5, 1) + dnorm(1:time_size, 15, 1)) * 30 else (dnorm(1:time_size, 5, 1) * 45 + dnorm(1:time_size, 15, 1) * 35),
    
    var3 = (1 - pnorm(1:time_size, 5, 1)) * c(5, 10, 15)[group],
    var4 = dnorm(1:time_size, 5, 1) * c(15, 30, 45)[group],
    
    var5 = (1 - pnorm(1:time_size, 5, 1)) * c(10, 20, 30)[group],
    var6 = dnorm(1:time_size, 5, 1) * c(40, 80, 120)[group],
    
    var7 = if (group == 1) dnorm(1:time_size, 4, 1) * 35 else if (group == 2) dnorm(1:time_size, 7, 1) * 45 else (0.4 * dnorm(1:time_size, 4, 1) + 0.6 * dnorm(1:time_size, 7, 1)) * 40
    
  )
}

df_traj_base <- bind_rows(
  generate_base_traj(1, 20),
  generate_base_traj(2, 20),
  generate_base_traj(3, 20)
)

df_traj_base %>% head

## Plotting: population and sample trajectories for all settings
vars <- c("var1", "var2", "var3", "var4", "var5", "var6", "var7", "var8", "var9")
fig_list <- list()

for (i in 1:length(vars)) {
  
  p <- ggplot(data = df_traj, aes(x = time, y = .data[[paste0("var", i)]], group = id)) +
    geom_line(alpha = 0.5) + 
    facet_grid(. ~ group) +
    xlab("Time") + 
    ylab(paste0("Var ", i)) +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          strip.text = element_text(size = 8))
  
  if (!(i %in% 8:9)) {
    p <- p + 
      geom_line(data = df_traj_base, aes(x = time, y = .data[[paste0("var", i)]], group = group), 
                color = "red", linewidth = 0.9)
  }
  
  fig_list[[i]] <- p
}

ggarrange(plotlist = fig_list[1:9], nrow = 3, ncol = 3, align = 'hv')


## -------------------------------------------------------------------------- ##
## Results of SFKmL
## -------------------------------------------------------------------------- ##

df_sim <- df_traj %>% select(-group)

## Compute scale values for each variable
vars_sim <- c("time", "var1", "var2", "var3", "var4", "var5", "var6", 
              "var7", "var8", "var9")  # time + 9 variables
scales_sim <- NULL
for (var in vars_sim) {
  scale_var = 1/diff(range(df_sim[[var]]))
  scales_sim = c(scales_sim, scale_var)
}

## Candidate time-scaling and l1bound
const_c_sim <- c(0.01, 0.1, 0.25, 0.5, 0.75, 1, seq(1.5, 4.5, by = 1))
lambda_sim <- scales_sim[1] * const_c_sim
l1b_sim <- seq(1.05, 2.95, length.out = 10)

## Loop over λ values to compute the distance array, apply SFclust(), 
## calculate the gap statistics for each
dist_ary_sim_list <- vector("list", length(const_c_sim))
res_sfperm_sim_list <- vector("list", length(const_c_sim))
df_gaps_sim_list <- vector("list", length(const_c_sim))

for (j in seq_along(const_c_sim)) {
  dist_ary_sim_list[[j]] <- dist.array(dt = df_sim, 
                                       time_scale = scales_sim[1] * const_c_sim[j], 
                                       var_scales = scales_sim[-1])
  
  res_sfperm_sim_list[[j]] <- SFclust.permute(dist.ary = dist_ary_sim_list[[j]], k = 3, 
                                              nperms = 50, seed = 10, 
                                              l1b = l1b_sim, 
                                              plot.gap = FALSE, plot.gap.l1b = FALSE)
  
  df_gaps_sim_list[[j]] <- data.frame(
    const_c_sim = rep(const_c_sim[j], length(res_sfperm_sim_list[[j]]$l1bounds)),
    lam = rep(lambda_sim[j], length(res_sfperm_sim_list[[j]]$l1bounds)),
    l1b = res_sfperm_sim_list[[j]]$l1bounds,
    gap = res_sfperm_sim_list[[j]]$gaps
  )
}

## Find the best (c, s) combination based on maximum gap
df_gaps_sim <- bind_rows(df_gaps_sim_list)
df_gaps_sim[which.max(df_gaps_sim$gap), ]
df_gaps_sim_max_by_lam <- df_gaps_sim %>%
  group_by(const_c_sim, lam) %>%
  summarise(gap_max = max(gap), .groups = "drop")

c_optimal_sim <- df_gaps_sim$const_c_sim[which.max(df_gaps_sim$gap)]
l1b_optimal_sim <- df_gaps_sim$l1b[which.max(df_gaps_sim$gap)]
c_opt_index_sim <- which(const_c_sim == c_optimal_sim)

## Apply SFKmL with optimal (λ, s)
SFclust_sim_optimal <- SFclust(k = 3, seed = 10, l1bound = l1b_optimal_sim, 
                               dist.ary = dist_ary_sim_list[[c_opt_index_sim]])

## Plotting: gap vs. const_c and variable weights vs. l1bound
weight_sim <- list()
for (j in seq_along(l1b_sim)) {
  dist_ary <- dist_ary_sim_list[[c_opt_index_sim]]
  res_sfclust <- SFclust(k = 3, seed = 10, l1bound = l1b_sim[j], dist.ary = dist_ary)
  weight_sim[[j]] <- res_sfclust$final.weight
}
df_weight_sim <- do.call(rbind, weight_sim) %>% data.frame
rownames(df_weight_sim) <- l1b_sim
colnames(df_weight_sim) <- colnames(df_sim)[-(1:2)]

df_weight_sim_long <- melt(data.frame(l1b_sim, df_weight_sim), id.vars = "l1b_sim")

p_gap_sim <- df_gaps_sim %>% ggplot(aes(x = const_c_sim, y = gap)) +
  geom_point() +
  geom_line(data = df_gaps_sim_max_by_lam, aes(x = const_c_sim, y = gap_max)) +
  geom_vline(xintercept = c_optimal_sim, size = 0.7, linetype = "dashed") +
  labs(x = "const_c (c)", y = "Gap") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray90"))

p_wght_sim <- ggplot(df_weight_sim_long, aes(x = l1b_sim, y = value, color = variable)) +
  geom_line() +
  geom_vline(xintercept = l1b_sim, linetype = "dashed", color = "gray90") +
  geom_vline(xintercept = l1b_optimal_sim, size = 0.7, linetype = "dashed", color = "gray10") +
  scale_x_continuous(breaks = sort(c(1, 1.5, 2, 2.5))) +
  labs(x = "l1bound (s)",
       y = "Weight",
       color = "Variable") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggarrange(p_gap_sim, p_wght_sim, nrow = 1, ncol = 2, align = 'hv')


## -------------------------------------------------------------------------- ##
## Clustering performance using ARI
## -------------------------------------------------------------------------- ##
## Compute the Adjusted Rand Index (ARI) for both clustering results (MFKmL and SFKmL)
clust_mf <- mfkml(dt = df_sim, clt_n = 3, scales = scales_sim, 
                  weight = 1, seed = 10, maxIter = 50)$Cluster[, 2]
clust_sf <- SFclust_sim_optimal$clust

unique_id_grp <- df_traj %>% distinct(id, group, .keep_all = FALSE)

ARI_res <- c(ARI_mf = randIndex(table(unique_id_grp$group, clust_mf)),
             ARI_sf = randIndex(table(unique_id_grp$group, clust_sf)))
ARI_res